/**************************************************************************/
/*  instance_data.cpp                                                     */
/**************************************************************************/
/*                         This file is part of:                          */
/*                             GODOT ENGINE                               */
/*                        https://godotengine.org                         */
/**************************************************************************/
/* Copyright (c) 2014-present Godot Engine contributors (see AUTHORS.md). */
/* Copyright (c) 2007-2014 Juan Linietsky, Ariel Manzur.                  */
/*                                                                        */
/* Permission is hereby granted, free of charge, to any person obtaining  */
/* a copy of this software and associated documentation files (the        */
/* "Software"), to deal in the Software without restriction, including    */
/* without limitation the rights to use, copy, modify, merge, publish,    */
/* distribute, sublicense, and/or sell copies of the Software, and to     */
/* permit persons to whom the Software is furnished to do so, subject to  */
/* the following conditions:                                              */
/*                                                                        */
/* The above copyright notice and this permission notice shall be         */
/* included in all copies or substantial portions of the Software.        */
/*                                                                        */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,        */
/* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF     */
/* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. */
/* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY   */
/* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,   */
/* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE      */
/* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                 */
/**************************************************************************/

#include "instance_data.h"
#include "../constants/voxel_constants.h"
#include "../util/math/funcs.h"
#include "../util/serialization.h"
#include <core/variant/variant.h>

namespace {
const uint32_t TRAILING_MAGIC = 0x900df00d;
}

const float VoxelInstanceBlockData::POSITION_RANGE_MINIMUM = 0.01f;

const float VoxelInstanceBlockData::SIMPLE_11B_V1_SCALE_RANGE_MINIMUM = 0.01f;

// TODO Unify with functions from VoxelBuffer?

inline uint8_t norm_to_u8(float x) {
	return clamp(static_cast<int>(128.f * x + 128.f), 0, 0xff);
}

inline float u8_to_norm(uint8_t v) {
	return (static_cast<real_t>(v) - 0x7f) * VoxelConstants::INV_0x7f;
}

struct CompressedQuaternion4b {
	uint8_t x;
	uint8_t y;
	uint8_t z;
	uint8_t w;

	static CompressedQuaternion4b from_quat(Quaternion q) {
		CompressedQuaternion4b c;
		c.x = norm_to_u8(q.x);
		c.y = norm_to_u8(q.y);
		c.z = norm_to_u8(q.z);
		c.w = norm_to_u8(q.w);
		return c;
	}

	Quaternion to_quat() const {
		Quaternion q;
		q.x = u8_to_norm(x);
		q.y = u8_to_norm(y);
		q.z = u8_to_norm(z);
		q.w = u8_to_norm(w);
		return q.normalized();
	}
};

bool serialize_instance_block_data(const VoxelInstanceBlockData &src,
		std::vector<uint8_t> &dst) {
	const uint8_t version = 0;
	const uint8_t instance_format = VoxelInstanceBlockData::FORMAT_SIMPLE_11B_V1;

	// TODO Apparently big-endian is dead
	// I chose it originally to match "network byte order",
	// but as I read comments about it there seem to be no reason to continue
	// using it. Needs a version increment.
	VoxelUtility::MemoryWriter w(dst, VoxelUtility::ENDIANESS_BIG_ENDIAN);

	ERR_FAIL_COND_V(src.position_range < 0.f, false);
	const float position_range =
			max(src.position_range, VoxelInstanceBlockData::POSITION_RANGE_MINIMUM);

	w.store_8(version);
	w.store_8(src.layers.size());
	w.store_float(position_range);

	// TODO Introduce a margin to position coordinates, stuff can spawn offset
	// from the ground. Or just compute the ranges
	const float pos_norm_scale = 1.f / position_range;

	for (size_t i = 0; i < src.layers.size(); ++i) {
		const VoxelInstanceBlockData::LayerData &layer = src.layers[i];

		ERR_FAIL_COND_V(layer.scale_max < layer.scale_min, false);

		float scale_min = layer.scale_min;
		float scale_max = layer.scale_max;
		if (scale_max - scale_min <
				VoxelInstanceBlockData::SIMPLE_11B_V1_SCALE_RANGE_MINIMUM) {
			scale_min = layer.scale_min;
			scale_max =
					scale_min + VoxelInstanceBlockData::SIMPLE_11B_V1_SCALE_RANGE_MINIMUM;
		}

		w.store_16(layer.id);
		w.store_16(layer.instances.size());
		w.store_float(scale_min);
		w.store_float(scale_max);
		w.store_8(instance_format);

		const float scale_norm_scale = 1.f / (scale_max - scale_min);

		for (size_t j = 0; j < layer.instances.size(); ++j) {
			const VoxelInstanceBlockData::InstanceData &instance = layer.instances[j];

			w.store_16(static_cast<uint16_t>(pos_norm_scale *
					instance.transform.origin.x * 0xffff));
			w.store_16(static_cast<uint16_t>(pos_norm_scale *
					instance.transform.origin.y * 0xffff));
			w.store_16(static_cast<uint16_t>(pos_norm_scale *
					instance.transform.origin.z * 0xffff));

			const float scale = instance.transform.get_basis().get_scale().y;
			w.store_8(
					static_cast<uint8_t>(scale_norm_scale * (scale - scale_min) * 0xff));

			const Quaternion q =
					instance.transform.get_basis().get_rotation_quaternion();
			const CompressedQuaternion4b cq = CompressedQuaternion4b::from_quat(q);
			w.store_8(cq.x);
			w.store_8(cq.y);
			w.store_8(cq.z);
			w.store_8(cq.w);
		}
	}

	w.store_32(TRAILING_MAGIC);

	return true;
}

bool deserialize_instance_block_data(VoxelInstanceBlockData &dst,
		Span<const uint8_t> src) {
	const uint8_t expected_version = 0;
	const uint8_t expected_instance_format =
			VoxelInstanceBlockData::FORMAT_SIMPLE_11B_V1;

	VoxelUtility::MemoryReader r(src, VoxelUtility::ENDIANESS_BIG_ENDIAN);

	const uint8_t version = r.get_8();
	ERR_FAIL_COND_V(version != expected_version, false);

	const uint8_t layers_count = r.get_8();
	dst.layers.resize(layers_count);

	dst.position_range = r.get_float();

	for (size_t i = 0; i < dst.layers.size(); ++i) {
		VoxelInstanceBlockData::LayerData &layer = dst.layers[i];

		layer.id = r.get_16();

		const uint16_t instance_count = r.get_16();
		layer.instances.resize(instance_count);

		layer.scale_min = r.get_float();
		layer.scale_max = r.get_float();
		ERR_FAIL_COND_V(layer.scale_max < layer.scale_min, false);
		const float scale_range = layer.scale_max - layer.scale_min;

		const uint8_t instance_format = r.get_8();
		ERR_FAIL_COND_V(instance_format != expected_instance_format, false);

		for (size_t j = 0; j < layer.instances.size(); ++j) {
			const float x =
					(static_cast<float>(r.get_16()) / 0xffff) * dst.position_range;
			const float y =
					(static_cast<float>(r.get_16()) / 0xffff) * dst.position_range;
			const float z =
					(static_cast<float>(r.get_16()) / 0xffff) * dst.position_range;

			const float s = (static_cast<float>(r.get_8()) / 0xff) * scale_range +
					layer.scale_min;

			CompressedQuaternion4b cq;
			cq.x = r.get_8();
			cq.y = r.get_8();
			cq.z = r.get_8();
			cq.w = r.get_8();
			const Quaternion q = cq.to_quat();

			VoxelInstanceBlockData::InstanceData &instance = layer.instances[j];
			instance.transform =
					Transform3D(Basis(q).scaled(Vector3(s, s, s)), Vector3(x, y, z));
		}
	}

	const uint32_t control_end = r.get_32();
	ERR_FAIL_COND_V_MSG(control_end != TRAILING_MAGIC, false,
			String("Expected {0}, found {1}")
					.format(varray(TRAILING_MAGIC, control_end)));

	return true;
}
