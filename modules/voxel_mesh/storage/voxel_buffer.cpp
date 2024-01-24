/**************************************************************************/
/*  voxel_buffer.cpp                                                      */
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

#define VOXEL_BUFFER_USE_MEMORY_POOL

#ifdef VOXEL_BUFFER_USE_MEMORY_POOL
#include "voxel_memory_pool.h"
#endif

#include "../edition/voxel_tool_buffer.h"
#include "../util/funcs.h"
#include "../util/profiling.h"
#include "voxel_buffer.h"

#include <core/io/image.h>
#include <core/io/marshalls.h>
#include <core/math/math_funcs.h>
#include <string.h>

namespace {
inline uint8_t *allocate_channel_data(uint32_t size) {
#ifdef VOXEL_BUFFER_USE_MEMORY_POOL
	return VoxelMemoryPool::get_singleton()->allocate(size);
#else
	return (uint8_t *)memalloc(size * sizeof(uint8_t));
#endif
}

inline void free_channel_data(uint8_t *data, uint32_t size) {
#ifdef VOXEL_BUFFER_USE_MEMORY_POOL
	VoxelMemoryPool::get_singleton()->recycle(data, size);
#else
	memfree(data);
#endif
}

uint64_t g_depth_max_values[] = {
	0xff, // 8
	0xffff, // 16
	0xffffffff, // 32
	0xffffffffffffffff // 64
};

inline uint32_t get_depth_bit_count(VoxelBuffer::Depth d) {
	CRASH_COND(d < 0 || d >= VoxelBuffer::DEPTH_COUNT);
	return VoxelBuffer::get_depth_byte_count(d) << 3;
}

inline uint64_t get_max_value_for_depth(VoxelBuffer::Depth d) {
	CRASH_COND(d < 0 || d >= VoxelBuffer::DEPTH_COUNT);
	return g_depth_max_values[d];
}

inline uint64_t clamp_value_for_depth(uint64_t value, VoxelBuffer::Depth d) {
	const uint64_t max_val = get_max_value_for_depth(d);
	if (value >= max_val) {
		return max_val;
	}
	return value;
}

static_assert(sizeof(uint32_t) == sizeof(float),
		"uint32_t and float cannot be marshalled back and forth");
static_assert(sizeof(uint64_t) == sizeof(double),
		"uint64_t and double cannot be marshalled back and forth");

inline uint64_t real_to_raw_voxel(real_t value, VoxelBuffer::Depth depth) {
	switch (depth) {
		case VoxelBuffer::DEPTH_8_BIT:
			return norm_to_u8(value);

		case VoxelBuffer::DEPTH_16_BIT:
			return norm_to_u16(value);

		case VoxelBuffer::DEPTH_32_BIT: {
			MarshallFloat m;
			m.f = value;
			return m.i;
		}
		case VoxelBuffer::DEPTH_64_BIT: {
			MarshallDouble m;
			m.d = value;
			return m.l;
		}
		default:
			CRASH_NOW();
			return 0;
	}
}

inline real_t raw_voxel_to_real(uint64_t value, VoxelBuffer::Depth depth) {
	// Depths below 32 are normalized between -1 and 1
	switch (depth) {
		case VoxelBuffer::DEPTH_8_BIT:
			return u8_to_norm(value);

		case VoxelBuffer::DEPTH_16_BIT:
			return u16_to_norm(value);

		case VoxelBuffer::DEPTH_32_BIT: {
			MarshallFloat m;
			m.i = value;
			return m.f;
		}

		case VoxelBuffer::DEPTH_64_BIT: {
			MarshallDouble m;
			m.l = value;
			return m.d;
		}

		default:
			CRASH_NOW();
			return 0;
	}
}
} // namespace

const char *VoxelBuffer::CHANNEL_ID_HINT_STRING =
		"Type,Sdf,Color,Indices,Weights,Data5,Data6,Data7";

VoxelBuffer::VoxelBuffer() {
	// Minecraft uses way more than 255 block types and there is room for eventual
	// metadata such as rotation
	_channels[CHANNEL_TYPE].depth = VoxelBuffer::DEFAULT_TYPE_CHANNEL_DEPTH;
	_channels[CHANNEL_TYPE].defval = 0;

	// 16-bit is better on average to handle large worlds
	_channels[CHANNEL_SDF].depth = VoxelBuffer::DEFAULT_SDF_CHANNEL_DEPTH;
	_channels[CHANNEL_SDF].defval = 0xffff;

	_channels[CHANNEL_INDICES].depth = VoxelBuffer::DEPTH_16_BIT;
	_channels[CHANNEL_INDICES].defval = encode_indices_to_packed_u16(0, 1, 2, 3);

	_channels[CHANNEL_WEIGHTS].depth = VoxelBuffer::DEPTH_16_BIT;
	_channels[CHANNEL_WEIGHTS].defval = encode_weights_to_packed_u16(15, 0, 0, 0);
}

VoxelBuffer::~VoxelBuffer() {
	clear();
}

void VoxelBuffer::create(unsigned int sx, unsigned int sy, unsigned int sz) {
	ERR_FAIL_COND(sx > MAX_SIZE || sy > MAX_SIZE || sz > MAX_SIZE);

	clear_voxel_metadata();

	VoxelVector3i new_size(sx, sy, sz);
	if (new_size != _size) {
		for (unsigned int i = 0; i < MAX_CHANNELS; ++i) {
			Channel &channel = _channels[i];
			if (channel.data) {
				// Channel already contained data
				delete_channel(i);
				create_channel(i, new_size, channel.defval);
			}
		}
		_size = new_size;
	}
}

void VoxelBuffer::create(VoxelVector3i size) {
	create(size.x, size.y, size.z);
}

void VoxelBuffer::clear() {
	for (unsigned int i = 0; i < MAX_CHANNELS; ++i) {
		Channel &channel = _channels[i];
		if (channel.data) {
			delete_channel(i);
		}
	}
	_size = VoxelVector3i();
	clear_voxel_metadata();
}

void VoxelBuffer::clear_channel(unsigned int channel_index,
		uint64_t clear_value) {
	ERR_FAIL_INDEX(channel_index, MAX_CHANNELS);
	Channel &channel = _channels[channel_index];
	if (channel.data != nullptr) {
		delete_channel(channel_index);
	}
	channel.defval = clamp_value_for_depth(clear_value, channel.depth);
}

void VoxelBuffer::clear_channel_f(unsigned int channel_index,
		real_t clear_value) {
	ERR_FAIL_INDEX(channel_index, MAX_CHANNELS);
	const Channel &channel = _channels[channel_index];
	clear_channel(channel_index, real_to_raw_voxel(clear_value, channel.depth));
}

void VoxelBuffer::set_default_values(
		FixedArray<uint64_t, VoxelBuffer::MAX_CHANNELS> values) {
	for (unsigned int i = 0; i < MAX_CHANNELS; ++i) {
		_channels[i].defval = clamp_value_for_depth(values[i], _channels[i].depth);
	}
}

uint64_t VoxelBuffer::get_voxel(int x, int y, int z,
		unsigned int channel_index) const {
	ERR_FAIL_INDEX_V(channel_index, MAX_CHANNELS, 0);
	ERR_FAIL_COND_V_MSG(
			!is_position_valid(x, y, z), 0,
			String("At position ({0}, {1}, {2})").format(varray(x, y, z)));

	const Channel &channel = _channels[channel_index];

	if (channel.data != nullptr) {
		const uint32_t i = get_index(x, y, z);

		switch (channel.depth) {
			case DEPTH_8_BIT:
				return channel.data[i];

			case DEPTH_16_BIT:
				return reinterpret_cast<uint16_t *>(channel.data)[i];

			case DEPTH_32_BIT:
				return reinterpret_cast<uint32_t *>(channel.data)[i];

			case DEPTH_64_BIT:
				return reinterpret_cast<uint64_t *>(channel.data)[i];

			default:
				CRASH_NOW();
				return 0;
		}

	} else {
		return channel.defval;
	}
}

void VoxelBuffer::set_voxel(uint64_t value, int x, int y, int z,
		unsigned int channel_index) {
	ERR_FAIL_INDEX(channel_index, MAX_CHANNELS);
	ERR_FAIL_COND_MSG(
			!is_position_valid(x, y, z),
			String("At position ({0}, {1}, {2})").format(varray(x, y, z)));

	Channel &channel = _channels[channel_index];

	value = clamp_value_for_depth(value, channel.depth);
	bool do_set = true;

	if (channel.data == nullptr) {
		if (channel.defval != value) {
			// Allocate channel with same initial values as defval
			create_channel(channel_index, _size, channel.defval);
		} else {
			do_set = false;
		}
	}

	if (do_set) {
		const uint32_t i = get_index(x, y, z);

		switch (channel.depth) {
			case DEPTH_8_BIT:
				channel.data[i] = value;
				break;

			case DEPTH_16_BIT:
				reinterpret_cast<uint16_t *>(channel.data)[i] = value;
				break;

			case DEPTH_32_BIT:
				reinterpret_cast<uint32_t *>(channel.data)[i] = value;
				break;

			case DEPTH_64_BIT:
				reinterpret_cast<uint64_t *>(channel.data)[i] = value;
				break;

			default:
				CRASH_NOW();
				break;
		}
	}
}

real_t VoxelBuffer::get_voxel_f(int x, int y, int z,
		unsigned int channel_index) const {
	ERR_FAIL_INDEX_V(channel_index, MAX_CHANNELS, 0);
	return raw_voxel_to_real(get_voxel(x, y, z, channel_index),
			_channels[channel_index].depth);
}

void VoxelBuffer::set_voxel_f(real_t value, int x, int y, int z,
		unsigned int channel_index) {
	ERR_FAIL_INDEX(channel_index, MAX_CHANNELS);
	set_voxel(real_to_raw_voxel(value, _channels[channel_index].depth), x, y, z,
			channel_index);
}

void VoxelBuffer::fill(uint64_t defval, unsigned int channel_index) {
	ERR_FAIL_INDEX(channel_index, MAX_CHANNELS);

	Channel &channel = _channels[channel_index];

	defval = clamp_value_for_depth(defval, channel.depth);

	if (channel.data == nullptr) {
		// Channel is already optimized and uniform
		if (channel.defval == defval) {
			// No change
			return;
		} else {
			// Just change default value
			channel.defval = defval;
			return;
		}
	}

	const unsigned int volume = get_volume();

	switch (channel.depth) {
		case DEPTH_8_BIT:
			memset(channel.data, defval, channel.size_in_bytes);
			break;

		case DEPTH_16_BIT:
			for (uint32_t i = 0; i < volume; ++i) {
				reinterpret_cast<uint16_t *>(channel.data)[i] = defval;
			}
			break;

		case DEPTH_32_BIT:
			for (uint32_t i = 0; i < volume; ++i) {
				reinterpret_cast<uint32_t *>(channel.data)[i] = defval;
			}
			break;

		case DEPTH_64_BIT:
			for (uint32_t i = 0; i < volume; ++i) {
				reinterpret_cast<uint64_t *>(channel.data)[i] = defval;
			}
			break;

		default:
			CRASH_NOW();
			break;
	}
}

void VoxelBuffer::fill_area(uint64_t defval, VoxelVector3i min,
		VoxelVector3i max, unsigned int channel_index) {
	ERR_FAIL_INDEX(channel_index, MAX_CHANNELS);

	VoxelVector3i::sort_min_max(min, max);
	min.clamp_to(VoxelVector3i(0, 0, 0), _size + VoxelVector3i(1, 1, 1));
	max.clamp_to(VoxelVector3i(0, 0, 0), _size + VoxelVector3i(1, 1, 1));
	const VoxelVector3i area_size = max - min;
	if (area_size.x == 0 || area_size.y == 0 || area_size.z == 0) {
		return;
	}

	Channel &channel = _channels[channel_index];
	defval = clamp_value_for_depth(defval, channel.depth);

	if (channel.data == nullptr) {
		if (channel.defval == defval) {
			return;
		} else {
			create_channel(channel_index, _size, channel.defval);
		}
	}

	VoxelVector3i pos;
	const unsigned int volume = get_volume();
	for (pos.z = min.z; pos.z < max.z; ++pos.z) {
		for (pos.x = min.x; pos.x < max.x; ++pos.x) {
			const unsigned int dst_ri = get_index(pos.x, pos.y + min.y, pos.z);
			CRASH_COND(dst_ri >= volume);

			switch (channel.depth) {
				case DEPTH_8_BIT:
					// Fill row by row
					memset(&channel.data[dst_ri], defval, area_size.y * sizeof(uint8_t));
					break;

				case DEPTH_16_BIT:
					for (int i = 0; i < area_size.y; ++i) {
						((uint16_t *)channel.data)[dst_ri + i] = defval;
					}
					break;

				case DEPTH_32_BIT:
					for (int i = 0; i < area_size.y; ++i) {
						((uint32_t *)channel.data)[dst_ri + i] = defval;
					}
					break;

				case DEPTH_64_BIT:
					for (int i = 0; i < area_size.y; ++i) {
						((uint64_t *)channel.data)[dst_ri + i] = defval;
					}
					break;

				default:
					CRASH_NOW();
					break;
			}
		}
	}
}

void VoxelBuffer::fill_area_f(float fvalue, VoxelVector3i min,
		VoxelVector3i max, unsigned int channel_index) {
	ERR_FAIL_INDEX(channel_index, MAX_CHANNELS);
	const Channel &channel = _channels[channel_index];
	fill_area(real_to_raw_voxel(fvalue, channel.depth), min, max, channel_index);
}

void VoxelBuffer::fill_f(real_t value, unsigned int channel) {
	ERR_FAIL_INDEX(channel, MAX_CHANNELS);
	fill(real_to_raw_voxel(value, _channels[channel].depth), channel);
}

template <typename T>
inline bool is_uniform_b(const uint8_t *data, unsigned int item_count) {
	return is_uniform<T>(reinterpret_cast<const T *>(data), item_count);
}

bool VoxelBuffer::is_uniform(unsigned int channel_index) const {
	ERR_FAIL_INDEX_V(channel_index, MAX_CHANNELS, true);

	const Channel &channel = _channels[channel_index];
	if (channel.data == nullptr) {
		// Channel has been optimized
		return true;
	}

	const unsigned int volume = get_volume();

	// Channel isn't optimized, so must look at each voxel
	switch (channel.depth) {
		case DEPTH_8_BIT:
			return ::is_uniform_b<uint8_t>(channel.data, volume);
		case DEPTH_16_BIT:
			return ::is_uniform_b<uint16_t>(channel.data, volume);
		case DEPTH_32_BIT:
			return ::is_uniform_b<uint32_t>(channel.data, volume);
		case DEPTH_64_BIT:
			return ::is_uniform_b<uint64_t>(channel.data, volume);
		default:
			CRASH_NOW();
			break;
	}

	return true;
}

void VoxelBuffer::compress_uniform_channels() {
	for (unsigned int i = 0; i < MAX_CHANNELS; ++i) {
		if (_channels[i].data != nullptr && is_uniform(i)) {
			// TODO More direct way
			const uint64_t v = get_voxel(0, 0, 0, i);
			clear_channel(i, v);
		}
	}
}

void VoxelBuffer::decompress_channel(unsigned int channel_index) {
	ERR_FAIL_INDEX(channel_index, MAX_CHANNELS);
	Channel &channel = _channels[channel_index];
	if (channel.data == nullptr) {
		create_channel(channel_index, _size, channel.defval);
	}
}

VoxelBuffer::Compression
VoxelBuffer::get_channel_compression(unsigned int channel_index) const {
	ERR_FAIL_INDEX_V(channel_index, MAX_CHANNELS, VoxelBuffer::COMPRESSION_NONE);
	const Channel &channel = _channels[channel_index];
	if (channel.data == nullptr) {
		return COMPRESSION_UNIFORM;
	}
	return COMPRESSION_NONE;
}

void VoxelBuffer::copy_format(const VoxelBuffer &other) {
	for (unsigned int i = 0; i < MAX_CHANNELS; ++i) {
		set_channel_depth(i, other.get_channel_depth(i));
	}
}

void VoxelBuffer::copy_from(const VoxelBuffer &other) {
	// Copy all channels, assuming sizes and formats match
	for (unsigned int i = 0; i < MAX_CHANNELS; ++i) {
		copy_from(other, i);
	}
}

void VoxelBuffer::copy_from(const VoxelBuffer &other,
		unsigned int channel_index) {
	ERR_FAIL_INDEX(channel_index, MAX_CHANNELS);
	ERR_FAIL_COND(other._size != _size);

	Channel &channel = _channels[channel_index];
	const Channel &other_channel = other._channels[channel_index];

	ERR_FAIL_COND(other_channel.depth != channel.depth);

	if (other_channel.data != nullptr) {
		if (channel.data == nullptr) {
			create_channel_noinit(channel_index, _size);
		}
		CRASH_COND(channel.size_in_bytes != other_channel.size_in_bytes);
		memcpy(channel.data, other_channel.data, channel.size_in_bytes);

	} else if (channel.data != nullptr) {
		delete_channel(channel_index);
	}

	channel.defval = other_channel.defval;
	channel.depth = other_channel.depth;
}

// TODO Disallow copying from overlapping areas of the same buffer
void VoxelBuffer::copy_from(const VoxelBuffer &other, VoxelVector3i src_min,
		VoxelVector3i src_max, VoxelVector3i dst_min,
		unsigned int channel_index) {
	ERR_FAIL_INDEX(channel_index, MAX_CHANNELS);

	Channel &channel = _channels[channel_index];
	const Channel &other_channel = other._channels[channel_index];

	ERR_FAIL_COND(other_channel.depth != channel.depth);

	if (channel.data == nullptr && other_channel.data == nullptr &&
			channel.defval == other_channel.defval) {
		// No action needed
		return;
	}

	if (other_channel.data != nullptr) {
		if (channel.data == nullptr) {
			// Note, we do this even if the pasted data happens to be all the same
			// value as our current channel. We assume that this case is not frequent
			// enough to bother, and compression can happen later
			create_channel(channel_index, _size, channel.defval);
		}
		const unsigned int item_size = get_depth_byte_count(channel.depth);
		Span<const uint8_t> src(other_channel.data, other_channel.size_in_bytes);
		Span<uint8_t> dst(channel.data, channel.size_in_bytes);
		copy_3d_region_zxy(dst, _size, dst_min, src, other._size, src_min, src_max,
				item_size);

	} else if (channel.defval != other_channel.defval) {
		// This logic is still required due to how source and destination regions
		// can be specified. The actual size of the destination area must be
		// determined from the source area, after it has been clipped.
		VoxelVector3i::sort_min_max(src_min, src_max);
		clip_copy_region(src_min, src_max, other._size, dst_min, _size);
		const VoxelVector3i area_size = src_max - src_min;
		if (area_size.x <= 0 || area_size.y <= 0 || area_size.z <= 0) {
			// Degenerate area, we'll not copy anything.
			return;
		}
		fill_area(other_channel.defval, dst_min, dst_min + area_size,
				channel_index);
	}
}

Ref<VoxelBuffer> VoxelBuffer::duplicate(bool include_metadata) const {
	VoxelBuffer *d = memnew(VoxelBuffer);
	d->create(_size);
	for (unsigned int i = 0; i < _channels.size(); ++i) {
		d->set_channel_depth(i, _channels[i].depth);
	}
	d->copy_from(*this);
	if (include_metadata) {
		d->copy_voxel_metadata(*this);
	}
	return Ref<VoxelBuffer>(d);
}

bool VoxelBuffer::get_channel_raw(unsigned int channel_index,
		Span<uint8_t> &slice) const {
	const Channel &channel = _channels[channel_index];
	if (channel.data != nullptr) {
		slice = Span<uint8_t>(channel.data, 0, channel.size_in_bytes);
		return true;
	}
	slice = Span<uint8_t>();
	return false;
}

void VoxelBuffer::create_channel(int i, VoxelVector3i size, uint64_t defval) {
	create_channel_noinit(i, size);
	fill(defval, i);
}

uint32_t VoxelBuffer::get_size_in_bytes_for_volume(VoxelVector3i size,
		Depth depth) {
	// Calculate appropriate size based on bit depth
	const unsigned int volume = size.x * size.y * size.z;
	const unsigned int bits = volume * ::get_depth_bit_count(depth);
	const unsigned int size_in_bytes = (bits >> 3);
	return size_in_bytes;
}

void VoxelBuffer::create_channel_noinit(int i, VoxelVector3i size) {
	Channel &channel = _channels[i];
	uint32_t size_in_bytes = get_size_in_bytes_for_volume(size, channel.depth);
	CRASH_COND(channel.data != nullptr);
	channel.data = allocate_channel_data(size_in_bytes);
	channel.size_in_bytes = size_in_bytes;
}

void VoxelBuffer::delete_channel(int i) {
	Channel &channel = _channels[i];
	ERR_FAIL_COND(channel.data == nullptr);
	free_channel_data(channel.data, channel.size_in_bytes);
	channel.data = nullptr;
	channel.size_in_bytes = 0;
}

void VoxelBuffer::downscale_to(VoxelBuffer &dst, VoxelVector3i src_min,
		VoxelVector3i src_max,
		VoxelVector3i dst_min) const {
	// TODO Align input to multiple of two

	src_min.clamp_to(VoxelVector3i(), _size);
	src_max.clamp_to(VoxelVector3i(), _size + VoxelVector3i(1));

	VoxelVector3i dst_max = dst_min + ((src_max - src_min) >> 1);

	// TODO This will be wrong if it overlaps the border?
	dst_min.clamp_to(VoxelVector3i(), dst._size);
	dst_max.clamp_to(VoxelVector3i(), dst._size + VoxelVector3i(1));

	for (int channel_index = 0; channel_index < MAX_CHANNELS; ++channel_index) {
		const Channel &src_channel = _channels[channel_index];
		const Channel &dst_channel = dst._channels[channel_index];

		if (src_channel.data == nullptr && dst_channel.data == nullptr &&
				src_channel.defval == dst_channel.defval) {
			// No action needed
			continue;
		}

		// Nearest-neighbor downscaling

		VoxelVector3i pos;
		for (pos.z = dst_min.z; pos.z < dst_max.z; ++pos.z) {
			for (pos.x = dst_min.x; pos.x < dst_max.x; ++pos.x) {
				for (pos.y = dst_min.y; pos.y < dst_max.y; ++pos.y) {
					const VoxelVector3i src_pos = src_min + ((pos - dst_min) << 1);

					// TODO Remove check once it works
					CRASH_COND(!is_position_valid(src_pos.x, src_pos.y, src_pos.z));

					uint64_t v;
					if (src_channel.data) {
						// TODO Optimized version?
						v = get_voxel(src_pos, channel_index);
					} else {
						v = src_channel.defval;
					}

					dst.set_voxel(v, pos, channel_index);
				}
			}
		}
	}
}

Ref<VoxelTool> VoxelBuffer::get_voxel_tool() {
	// I can't make this function `const`, because `Ref<T>` has no constructor
	// taking a `const T*`. The compiler would then choose Ref<T>(const Variant&),
	// which fumbles `this` into a null pointer
	Ref<VoxelBuffer> vb(this);
	return Ref<VoxelTool>(memnew(VoxelToolBuffer(vb)));
}

bool VoxelBuffer::equals(const VoxelBuffer &p_other) const {
	if (p_other._size != _size) {
		return false;
	}

	for (int channel_index = 0; channel_index < MAX_CHANNELS; ++channel_index) {
		const Channel &channel = _channels[channel_index];
		const Channel &other_channel = p_other._channels[channel_index];

		if ((channel.data == nullptr) != (other_channel.data == nullptr)) {
			// Note: they could still logically be equal if one channel contains
			// uniform voxel memory
			return false;
		}

		if (channel.depth != other_channel.depth) {
			return false;
		}

		if (channel.data == nullptr) {
			if (channel.defval != other_channel.defval) {
				return false;
			}

		} else {
			ERR_FAIL_COND_V(channel.size_in_bytes != other_channel.size_in_bytes,
					false);
			for (unsigned int i = 0; i < channel.size_in_bytes; ++i) {
				if (channel.data[i] != other_channel.data[i]) {
					return false;
				}
			}
		}
	}

	return true;
}

void VoxelBuffer::set_channel_depth(unsigned int channel_index,
		Depth new_depth) {
	ERR_FAIL_INDEX(channel_index, MAX_CHANNELS);
	ERR_FAIL_INDEX(new_depth, DEPTH_COUNT);
	Channel &channel = _channels[channel_index];
	if (channel.depth == new_depth) {
		return;
	}
	if (channel.data != nullptr) {
		// TODO Implement conversion and do it when specified
		WARN_PRINT("Changing VoxelBuffer depth with present data, this will reset "
				   "the channel");
		delete_channel(channel_index);
	}
	channel.defval = clamp_value_for_depth(channel.defval, new_depth);
	channel.depth = new_depth;
}

VoxelBuffer::Depth
VoxelBuffer::get_channel_depth(unsigned int channel_index) const {
	ERR_FAIL_INDEX_V(channel_index, MAX_CHANNELS, DEPTH_8_BIT);
	return _channels[channel_index].depth;
}

uint32_t VoxelBuffer::get_depth_bit_count(Depth d) {
	return ::get_depth_bit_count(d);
}

float VoxelBuffer::get_sdf_quantization_scale(Depth d) {
	switch (d) {
		// Normalized
		case DEPTH_8_BIT:
			return VoxelConstants::QUANTIZED_SDF_8_BITS_SCALE;
		case DEPTH_16_BIT:
			return VoxelConstants::QUANTIZED_SDF_16_BITS_SCALE;
		// Direct
		default:
			return 1.f;
	}
}

void VoxelBuffer::set_block_metadata(Variant meta) {
	_block_metadata = meta;
}

Variant VoxelBuffer::get_voxel_metadata(VoxelVector3i pos) const {
	ERR_FAIL_COND_V(!is_position_valid(pos), Variant());
	const HashMap<VoxelVector3i, Variant, Vector3iHasher>::ConstIterator elem = _voxel_metadata.find(pos);
	if (elem) {
		return elem->value;
	} else {
		return Variant();
	}
}

void VoxelBuffer::set_voxel_metadata(VoxelVector3i pos, Variant meta) {
	ERR_FAIL_COND(!is_position_valid(pos));
	if (meta.get_type() == Variant::NIL) {
		_voxel_metadata.erase(pos);
	} else {
		_voxel_metadata[pos] = meta;
	}
}

void VoxelBuffer::for_each_voxel_metadata(Callable callback) const {
	ERR_FAIL_COND(callback.is_null());
	for (const KeyValue<VoxelVector3i, Variant> &elem : _voxel_metadata) {
		const Variant key = elem.key.to_vec3();
		const Variant *args[2] = { &key, &elem.value };
		Callable::CallError err;
		Variant retval;
		callback.callp(args, 2, retval, err);
		ERR_FAIL_COND_MSG(err.error != Callable::CallError::CALL_OK,
				String("FuncRef call failed at {0}").format(varray(key)));
		// TODO Can't provide detailed error because FuncRef doesn't give us access
		// to the object ERR_FAIL_COND_MSG(err.error != Variant::CallError::CALL_OK,
		// false, 		Variant::get_call_error_text(callback->get_object(),
		// method_name, nullptr, 0, err));
	}
}

void VoxelBuffer::for_each_voxel_metadata_in_area(Callable callback,
		Box3i box) const {
	ERR_FAIL_COND(callback.is_null());
	for_each_voxel_metadata_in_area(box, [&callback](VoxelVector3i pos, Variant meta) {
		const Variant key = pos.to_vec3();
		const Variant *args[2] = { &key, &meta };
		Callable::CallError err;
		Variant retval;
		callback.callp(args, 2, retval, err);

		ERR_FAIL_COND_MSG(err.error != Callable::CallError::CALL_OK,
				String("FuncRef call failed at {0}").format(varray(key)));
		// TODO Can't provide detailed error because FuncRef doesn't give us access
		// to the object ERR_FAIL_COND_MSG(err.error != Variant::CallError::CALL_OK,
		// false, 		Variant::get_call_error_text(callback->get_object(),
		// method_name, nullptr, 0, err));
	});
}

void VoxelBuffer::clear_voxel_metadata() {
	_voxel_metadata.clear();
}

void VoxelBuffer::clear_voxel_metadata_in_area(Box3i box) {
	for (const KeyValue<VoxelVector3i, Variant> &elem : _voxel_metadata) {
		if (box.contains(elem.key)) {
			_voxel_metadata.erase(elem.key);
		}
	}
}

void VoxelBuffer::copy_voxel_metadata_in_area(Ref<VoxelBuffer> src_buffer,
		Box3i src_box,
		VoxelVector3i dst_origin) {
	ERR_FAIL_COND(src_buffer.is_null());
	ERR_FAIL_COND(!src_buffer->is_box_valid(src_box));

	const Box3i clipped_src_box =
			src_box.clipped(Box3i(src_box.pos - dst_origin, _size));
	const VoxelVector3i clipped_dst_offset =
			dst_origin + clipped_src_box.pos - src_box.pos;

	for (const KeyValue<VoxelVector3i, Variant> &elem : src_buffer->_voxel_metadata) {
		const VoxelVector3i src_pos = elem.key;
		if (src_box.contains(src_pos)) {
			const VoxelVector3i dst_pos = src_pos + clipped_dst_offset;
			CRASH_COND(!is_position_valid(dst_pos));
			_voxel_metadata[dst_pos] = elem.value.duplicate();
		}
	}
}

void VoxelBuffer::copy_voxel_metadata(const VoxelBuffer &src_buffer) {
	ERR_FAIL_COND(src_buffer.get_size() != _size);

	for (const KeyValue<VoxelVector3i, Variant> &elem : src_buffer._voxel_metadata) {
		const VoxelVector3i pos = elem.key;
		_voxel_metadata[pos] = elem.value.duplicate();
	}

	_block_metadata = src_buffer._block_metadata.duplicate();
}

Ref<Image> VoxelBuffer::debug_print_sdf_to_image_top_down() {
	Ref<Image> im = Image::create_empty(_size.x, _size.z, false, Image::FORMAT_RGB8);
	VoxelVector3i pos;
	for (pos.z = 0; pos.z < _size.z; ++pos.z) {
		for (pos.x = 0; pos.x < _size.x; ++pos.x) {
			for (pos.y = _size.y - 1; pos.y >= 0; --pos.y) {
				float v = get_voxel_f(pos.x, pos.y, pos.z, CHANNEL_SDF);
				if (v < 0.0) {
					break;
				}
			}
			float h = pos.y;
			float c = h / _size.y;
			im->set_pixel(pos.x, pos.z, Color(c, c, c));
		}
	}
	return Ref<Image>(im);
}

void VoxelBuffer::_bind_methods() {
	ClassDB::bind_method(D_METHOD("create", "sx", "sy", "sz"),
			&VoxelBuffer::_b_create);
	ClassDB::bind_method(D_METHOD("clear"), &VoxelBuffer::clear);

	ClassDB::bind_method(D_METHOD("get_size"), &VoxelBuffer::_b_get_size);
	ClassDB::bind_method(D_METHOD("get_size_x"), &VoxelBuffer::get_size_x);
	ClassDB::bind_method(D_METHOD("get_size_y"), &VoxelBuffer::get_size_y);
	ClassDB::bind_method(D_METHOD("get_size_z"), &VoxelBuffer::get_size_z);

	ClassDB::bind_method(D_METHOD("set_voxel", "value", "x", "y", "z", "channel"),
			&VoxelBuffer::_b_set_voxel, DEFVAL(0));
	ClassDB::bind_method(
			D_METHOD("set_voxel_f", "value", "x", "y", "z", "channel"),
			&VoxelBuffer::_b_set_voxel_f, DEFVAL(0));
	ClassDB::bind_method(D_METHOD("set_voxel_v", "value", "pos", "channel"),
			&VoxelBuffer::_b_set_voxel_v, DEFVAL(0));
	ClassDB::bind_method(D_METHOD("get_voxel", "x", "y", "z", "channel"),
			&VoxelBuffer::_b_get_voxel, DEFVAL(0));
	ClassDB::bind_method(D_METHOD("get_voxel_f", "x", "y", "z", "channel"),
			&VoxelBuffer::get_voxel_f, DEFVAL(0));
	ClassDB::bind_method(D_METHOD("get_voxel_tool"),
			&VoxelBuffer::get_voxel_tool);

	ClassDB::bind_method(D_METHOD("get_channel_depth", "channel"),
			&VoxelBuffer::get_channel_depth);
	ClassDB::bind_method(D_METHOD("set_channel_depth", "channel", "depth"),
			&VoxelBuffer::set_channel_depth);

	ClassDB::bind_method(D_METHOD("fill", "value", "channel"), &VoxelBuffer::fill,
			DEFVAL(0));
	ClassDB::bind_method(D_METHOD("fill_f", "value", "channel"),
			&VoxelBuffer::fill_f, DEFVAL(0));
	ClassDB::bind_method(D_METHOD("fill_area", "value", "min", "max", "channel"),
			&VoxelBuffer::_b_fill_area, DEFVAL(0));
	ClassDB::bind_method(D_METHOD("copy_channel_from", "other", "channel"),
			&VoxelBuffer::_b_copy_channel_from);
	ClassDB::bind_method(D_METHOD("copy_channel_from_area", "other", "src_min",
								 "src_max", "dst_min", "channel"),
			&VoxelBuffer::_b_copy_channel_from_area);
	ClassDB::bind_method(
			D_METHOD("downscale_to", "dst", "src_min", "src_max", "dst_min"),
			&VoxelBuffer::_b_downscale_to);

	ClassDB::bind_method(D_METHOD("is_uniform", "channel"),
			&VoxelBuffer::is_uniform);
	// TODO Rename `compress_uniform_channels`
	ClassDB::bind_method(D_METHOD("optimize"),
			&VoxelBuffer::compress_uniform_channels);
	ClassDB::bind_method(D_METHOD("get_channel_compression", "channel"),
			&VoxelBuffer::get_channel_compression);

	ClassDB::bind_method(D_METHOD("get_block_metadata"),
			&VoxelBuffer::get_block_metadata);
	ClassDB::bind_method(D_METHOD("set_block_metadata", "meta"),
			&VoxelBuffer::set_block_metadata);
	ClassDB::bind_method(D_METHOD("get_voxel_metadata", "pos"),
			&VoxelBuffer::_b_get_voxel_metadata);
	ClassDB::bind_method(D_METHOD("set_voxel_metadata", "pos", "value"),
			&VoxelBuffer::_b_set_voxel_metadata);
	ClassDB::bind_method(D_METHOD("for_each_voxel_metadata", "callback"),
			&VoxelBuffer::for_each_voxel_metadata);
	ClassDB::bind_method(D_METHOD("for_each_voxel_metadata_in_area", "callback",
								 "min_pos", "max_pos"),
			&VoxelBuffer::_b_for_each_voxel_metadata_in_area);
	ClassDB::bind_method(D_METHOD("clear_voxel_metadata"),
			&VoxelBuffer::clear_voxel_metadata);
	ClassDB::bind_method(
			D_METHOD("clear_voxel_metadata_in_area", "min_pos", "max_pos"),
			&VoxelBuffer::_b_clear_voxel_metadata_in_area);
	ClassDB::bind_method(D_METHOD("copy_voxel_metadata_in_area", "src_buffer",
								 "src_min_pos", "src_max_pos", "dst_min_pos"),
			&VoxelBuffer::_b_copy_voxel_metadata_in_area);

	BIND_ENUM_CONSTANT(CHANNEL_TYPE);
	BIND_ENUM_CONSTANT(CHANNEL_SDF);
	BIND_ENUM_CONSTANT(CHANNEL_COLOR);
	BIND_ENUM_CONSTANT(CHANNEL_INDICES);
	BIND_ENUM_CONSTANT(CHANNEL_WEIGHTS);
	BIND_ENUM_CONSTANT(CHANNEL_DATA5);
	BIND_ENUM_CONSTANT(CHANNEL_DATA6);
	BIND_ENUM_CONSTANT(CHANNEL_DATA7);
	BIND_ENUM_CONSTANT(MAX_CHANNELS);

	BIND_ENUM_CONSTANT(DEPTH_8_BIT);
	BIND_ENUM_CONSTANT(DEPTH_16_BIT);
	BIND_ENUM_CONSTANT(DEPTH_32_BIT);
	BIND_ENUM_CONSTANT(DEPTH_64_BIT);
	BIND_ENUM_CONSTANT(DEPTH_COUNT);

	BIND_ENUM_CONSTANT(COMPRESSION_NONE);
	BIND_ENUM_CONSTANT(COMPRESSION_UNIFORM);
	BIND_ENUM_CONSTANT(COMPRESSION_COUNT);

	BIND_CONSTANT(MAX_SIZE);
}

void VoxelBuffer::_b_copy_channel_from(Ref<VoxelBuffer> other,
		unsigned int channel) {
	ERR_FAIL_COND(other.is_null());
	copy_from(**other, channel);
}

void VoxelBuffer::_b_copy_channel_from_area(Ref<VoxelBuffer> other,
		Vector3 src_min, Vector3 src_max,
		Vector3 dst_min,
		unsigned int channel) {
	ERR_FAIL_COND(other.is_null());
	copy_from(**other, VoxelVector3i(src_min), VoxelVector3i(src_max),
			VoxelVector3i(dst_min), channel);
}

void VoxelBuffer::_b_downscale_to(Ref<VoxelBuffer> dst, Vector3 src_min,
		Vector3 src_max, Vector3 dst_min) const {
	ERR_FAIL_COND(dst.is_null());
	downscale_to(**dst, VoxelVector3i(src_min), VoxelVector3i(src_max),
			VoxelVector3i(dst_min));
}

void VoxelBuffer::_b_for_each_voxel_metadata_in_area(Callable callback,
		Vector3 min_pos,
		Vector3 max_pos) {
	for_each_voxel_metadata_in_area(
			callback,
			Box3i::from_min_max(VoxelVector3i(min_pos), VoxelVector3i(max_pos)));
}

void VoxelBuffer::_b_clear_voxel_metadata_in_area(Vector3 min_pos,
		Vector3 max_pos) {
	clear_voxel_metadata_in_area(
			Box3i::from_min_max(VoxelVector3i(min_pos), VoxelVector3i(max_pos)));
}

void VoxelBuffer::_b_copy_voxel_metadata_in_area(Ref<VoxelBuffer> src_buffer,
		Vector3 src_min_pos,
		Vector3 src_max_pos,
		Vector3 dst_pos) {
	copy_voxel_metadata_in_area(src_buffer,
			Box3i::from_min_max(VoxelVector3i(src_min_pos),
					VoxelVector3i(src_max_pos)),
			dst_pos);
}
