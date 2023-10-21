/**************************************************************************/
/*  csg_shape.h                                                           */
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

#ifndef CSG_SHAPE_H
#define CSG_SHAPE_H

#include "core/io/file_access.h"
#include "core/io/marshalls.h"
#include "core/math/math_defs.h"
#include "core/math/random_number_generator.h"
#include "core/object/class_db.h"
#include "core/object/object.h"
#include "core/os/memory.h"
#include "core/variant/variant.h"
#include "csg.h"

#include "scene/3d/path_3d.h"
#include "scene/3d/visual_instance_3d.h"
#include "scene/resources/concave_polygon_shape_3d.h"

#include "scene/resources/material.h"
#include "scene/resources/mesh.h"
#include "thirdparty/misc/mikktspace.h"
#include <cmath>

class CSGShape3D : public GeometryInstance3D {
	GDCLASS(CSGShape3D, GeometryInstance3D);

public:
	enum Operation {
		OPERATION_UNION,
		OPERATION_INTERSECTION,
		OPERATION_SUBTRACTION,

	};

private:
	Operation operation = OPERATION_UNION;
	CSGShape3D *parent_shape = nullptr;

	CSGBrush *brush = nullptr;

	AABB node_aabb;

	bool dirty = false;
	bool last_visible = false;
	float snap = 0.001;

	bool use_collision = false;
	uint32_t collision_layer = 1;
	uint32_t collision_mask = 1;
	real_t collision_priority = 1.0;
	Ref<ConcavePolygonShape3D> root_collision_shape;
	RID root_collision_instance;
	RID root_collision_debug_instance;
	Transform3D debug_shape_old_transform;

	bool calculate_tangents = true;

	Ref<ArrayMesh> root_mesh;

	struct Vector3Hasher {
		_ALWAYS_INLINE_ uint32_t hash(const Vector3 &p_vec3) const {
			uint32_t h = hash_murmur3_one_float(p_vec3.x);
			h = hash_murmur3_one_float(p_vec3.y, h);
			h = hash_murmur3_one_float(p_vec3.z, h);
			return h;
		}
	};

	struct ShapeUpdateSurface {
		Vector<Vector3> vertices;
		Vector<Vector3> normals;
		Vector<Vector2> uvs;
		Vector<real_t> tans;
		Ref<Material> material;
		int last_added = 0;

		Vector3 *verticesw = nullptr;
		Vector3 *normalsw = nullptr;
		Vector2 *uvsw = nullptr;
		real_t *tansw = nullptr;
	};

	//mikktspace callbacks
	static int mikktGetNumFaces(const SMikkTSpaceContext *pContext);
	static int mikktGetNumVerticesOfFace(const SMikkTSpaceContext *pContext, const int iFace);
	static void mikktGetPosition(const SMikkTSpaceContext *pContext, float fvPosOut[], const int iFace, const int iVert);
	static void mikktGetNormal(const SMikkTSpaceContext *pContext, float fvNormOut[], const int iFace, const int iVert);
	static void mikktGetTexCoord(const SMikkTSpaceContext *pContext, float fvTexcOut[], const int iFace, const int iVert);
	static void mikktSetTSpaceDefault(const SMikkTSpaceContext *pContext, const float fvTangent[], const float fvBiTangent[], const float fMagS, const float fMagT,
			const tbool bIsOrientationPreserving, const int iFace, const int iVert);

	void _update_shape();
	void _update_collision_faces();
	bool _is_debug_collision_shape_visible();
	void _update_debug_collision_shape();
	void _clear_debug_collision_shape();
	void _on_transform_changed();

protected:
	void _notification(int p_what);
	virtual CSGBrush *_build_brush() = 0;
	void _make_dirty(bool p_parent_removing = false);

	static void _bind_methods();

	friend class CSGCombiner3D;
	CSGBrush *_get_brush();

	void _validate_property(PropertyInfo &p_property) const;

public:
	Array get_meshes() const;

	void set_operation(Operation p_operation);
	Operation get_operation() const;

	virtual Vector<Vector3> get_brush_faces();

	virtual AABB get_aabb() const override;

	void set_use_collision(bool p_enable);
	bool is_using_collision() const;

	void set_collision_layer(uint32_t p_layer);
	uint32_t get_collision_layer() const;

	void set_collision_mask(uint32_t p_mask);
	uint32_t get_collision_mask() const;

	void set_collision_layer_value(int p_layer_number, bool p_value);
	bool get_collision_layer_value(int p_layer_number) const;

	void set_collision_mask_value(int p_layer_number, bool p_value);
	bool get_collision_mask_value(int p_layer_number) const;

	void set_collision_priority(real_t p_priority);
	real_t get_collision_priority() const;

	void set_snap(float p_snap);
	float get_snap() const;

	void set_calculate_tangents(bool p_calculate_tangents);
	bool is_calculating_tangents() const;

	bool is_root_shape() const;
	CSGShape3D();
	~CSGShape3D();
};

VARIANT_ENUM_CAST(CSGShape3D::Operation)

class CSGCombiner3D : public CSGShape3D {
	GDCLASS(CSGCombiner3D, CSGShape3D);

private:
	virtual CSGBrush *_build_brush() override;

public:
	CSGCombiner3D();
};

class CSGPrimitive3D : public CSGShape3D {
	GDCLASS(CSGPrimitive3D, CSGShape3D);

protected:
	bool flip_faces;
	CSGBrush *_create_brush_from_arrays(const Vector<Vector3> &p_vertices, const Vector<Vector2> &p_uv, const Vector<bool> &p_smooth, const Vector<Ref<Material>> &p_materials);
	static void _bind_methods();

public:
	void set_flip_faces(bool p_invert);
	bool get_flip_faces();

	CSGPrimitive3D();
};

class CSGDeform3D;
class CSGMesh3D : public CSGPrimitive3D {
	GDCLASS(CSGMesh3D, CSGPrimitive3D);

	virtual CSGBrush *_build_brush() override;

	Ref<Mesh> mesh;
	Ref<Material> material;

	void _mesh_changed();
	friend CSGDeform3D;

protected:
	static void _bind_methods();

public:
	void set_mesh(const Ref<Mesh> &p_mesh);
	Ref<Mesh> get_mesh();

	void set_material(const Ref<Material> &p_material);
	Ref<Material> get_material() const;
};

class CSGDeform3D : public CSGMesh3D {
	GDCLASS(CSGDeform3D, CSGMesh3D);
	PackedVector3Array lattice;
	PackedVector3Array start_lattice;
	Vector3 lattice_size = Vector3(2, 2, 2);
	Vector3i lattice_resolution = Vector3i(9, 9, 9);
	bool fix_normals = true;
	bool smooth = false;

public:
	CSGDeform3D() {}

	void set_lattice(const PackedVector3Array &p_value) {
		lattice = p_value;
		_mesh_changed();
	}

	PackedVector3Array get_lattice() const {
		return lattice;
	}

	void set_lattice_size(const Vector3 &p_value) {
		lattice_size = p_value;
		_mesh_changed();
	}

	Vector3 get_lattice_size() const {
		return lattice_size;
	}

	void set_lattice_resolution(const Vector3i &p_value) {
		lattice_resolution = p_value.clamp(Vector3i(1, 1, 1), Vector3i(64, 64, 64));

		_mesh_changed();
	}

	Vector3i get_lattice_resolution() const {
		return lattice_resolution;
	}

	void set_fix_normals(bool value) {
		fix_normals = value;
		_mesh_changed();
	}

	bool get_fix_normals() const {
		return fix_normals;
	}

	void set_smooth(bool value) {
		smooth = value;
		_mesh_changed();
	}

	bool get_smooth() const {
		return smooth;
	}

	Vector3 lattice_get_fast(const Vector3 &p_coordinates) {
		Vector3 resolution = lattice_resolution;
		Vector3 resm1 = lattice_resolution - Vector3(1, 1, 1);

		Vector3 coord_mod = (p_coordinates / lattice_size + Vector3(0.5, 0.5, 0.5)) * resm1;
		Vector3 cf = coord_mod.floor();

		Vector3i a = Vector3i(cf).clamp(Vector3i(), resm1) * Vector3i(1, resolution.x, resolution.x * resolution.y);
		Vector3i b = Vector3i(cf + Vector3(1, 1, 1)).clamp(Vector3i(), resm1) * Vector3i(1, resolution.x, resolution.x * resolution.y);

		Vector3 t = (coord_mod - cf).clamp(Vector3(), Vector3(1, 1, 1));

		Vector3 aaa = lattice[a.z + a.y + a.x];
		Vector3 baa = lattice[a.z + a.y + b.x];
		Vector3 aba = lattice[a.z + b.y + a.x];
		Vector3 bba = lattice[a.z + b.y + b.x];
		Vector3 aab = lattice[b.z + a.y + a.x];
		Vector3 bab = lattice[b.z + a.y + b.x];
		Vector3 abb = lattice[b.z + b.y + a.x];
		Vector3 bbb = lattice[b.z + b.y + b.x];

		Vector3 aa_ = aaa.lerp(baa, t.x);
		Vector3 ba_ = aba.lerp(bba, t.x);
		Vector3 ab_ = aab.lerp(bab, t.x);
		Vector3 bb_ = abb.lerp(bbb, t.x);

		Vector3 a__ = aa_.lerp(ba_, t.y);
		Vector3 b__ = ab_.lerp(bb_, t.y);

		return a__.lerp(b__, t.z);
	}

	Vector3 lattice_get_smooth(const Vector3 &p_coordinates) {
		Vector3 resm1 = lattice_resolution - Vector3(1, 1, 1);

		Vector3 coord_mod = (p_coordinates / lattice_size + Vector3(0.5, 0.5, 0.5)) * resm1;
		Vector3 cf = coord_mod.floor();

		Vector3i a = Vector3i(cf).clamp(Vector3i(), resm1);
		Vector3i b = Vector3i(cf + Vector3(1, 1, 1)).clamp(Vector3i(), resm1);

		Vector3 t = (coord_mod - cf).clamp(Vector3(), Vector3(1, 1, 1));

		return lattice_get_i_interp_3d(a, b, t);
	}

	void lattice_get_weights(const Vector3 &p_coordinates, float p_amount, HashMap<Vector3i, real_t> &p_weights, HashMap<Vector3i, int> &p_counts) {
		Vector3 resm1 = lattice_resolution - Vector3(1, 1, 1);

		Vector3 coord_mod = (p_coordinates / lattice_size) * resm1 + resm1 * 0.5;
		Vector3 cf = coord_mod.floor();

		Vector3i a = Vector3i(cf).clamp(Vector3i(), resm1);
		Vector3i b = Vector3i(cf + Vector3(1, 1, 1)).clamp(Vector3i(), resm1);

		Vector3 t = (coord_mod - cf).clamp(Vector3(), Vector3(1, 1, 1));

		for (int z = a.z; z <= b.z; ++z) {
			float w_z = Math::lerp(a.z - z + 1, z - a.z, t.z);
			for (int y = a.y; y <= b.y; ++y) {
				float w_y = Math::lerp(a.y - y + 1, y - a.y, t.y) * w_z;
				for (int x = a.x; x <= b.x; ++x) {
					float w = Math::lerp(a.x - x + 1, x - a.x, t.x) * w_y;
					Vector3i vec(x, y, z);

					if (!p_weights.has(vec)) {
						p_weights[vec] = 0.0;
						p_counts[vec] = 0;
					}

					p_weights[vec] = real_t(p_weights[vec]) + (w * p_amount);
					p_counts[vec] = real_t(p_counts[vec]) + 1;
				}
			}
		}
	}

	Vector3 lattice_get_i_interp_1d(const Vector3i &p_c1, const Vector3i &p_c2, float p_t) {
		Vector3 res = lattice_resolution;

		Vector3 a = lattice[p_c1.z * res.x * res.y + p_c1.y * res.x + p_c1.x];
		Vector3 b = lattice[p_c2.z * res.x * res.y + p_c2.y * res.x + p_c2.x];

		Vector3 d = p_c2 - p_c1;

		Vector3i c0 = (p_c1 - d).clamp(Vector3i(), res - Vector3i(1, 1, 1));
		Vector3i c3 = (p_c2 + d).clamp(Vector3i(), res - Vector3i(1, 1, 1));

		Vector3 pre = lattice[c0.z * res.x * res.y + c0.y * res.x + c0.x];
		Vector3 post = lattice[c3.z * res.x * res.y + c3.y * res.x + c3.x];

		return a.cubic_interpolate(b, pre, post, p_t);
	}

	Vector3 lattice_get_i_interp_2d(const Vector3i &p_c1, const Vector3i &p_c2, float p_tx, float p_ty) {
		Vector3 res = lattice_resolution;

		Vector3i dy = (p_c2 - p_c1) * Vector3i(0, 1, 0);

		Vector3 a = lattice_get_i_interp_1d(p_c1, p_c2 - dy, p_tx);
		Vector3 b = lattice_get_i_interp_1d(p_c1 + dy, p_c2, p_tx);

		Vector3i c0_a = (p_c1 - dy).clamp(Vector3i(), res - Vector3i(1, 1, 1));
		Vector3i c0_b = (p_c2 - dy - dy).clamp(Vector3i(), res - Vector3i(1, 1, 1));

		Vector3i c1_a = (p_c1 + dy + dy).clamp(Vector3i(), res - Vector3i(1, 1, 1));
		Vector3i c1_b = (p_c2 + dy).clamp(Vector3i(), res - Vector3i(1, 1, 1));

		Vector3 pre = lattice_get_i_interp_1d(c0_a, c0_b, p_tx);
		Vector3 post = lattice_get_i_interp_1d(c1_a, c1_b, p_tx);

		return a.cubic_interpolate(b, pre, post, p_ty);
	}

	Vector3 lattice_get_i_interp_3d(const Vector3i &p_c1, const Vector3i &p_c2, const Vector3 &p_tv) {
		Vector3 res = lattice_resolution;

		Vector3i dz = (p_c2 - p_c1) * Vector3i(0, 0, 1);

		Vector3 a = lattice_get_i_interp_2d(p_c1, p_c2 - dz, p_tv.x, p_tv.y);
		Vector3 b = lattice_get_i_interp_2d(p_c1 + dz, p_c2, p_tv.x, p_tv.y);

		Vector3i c0_a = (p_c1 - dz).clamp(Vector3i(), res - Vector3i(1, 1, 1));
		Vector3i c0_b = (p_c2 - dz - dz).clamp(Vector3i(), res - Vector3i(1, 1, 1));

		Vector3i c1_a = (p_c1 + dz + dz).clamp(Vector3i(), res - Vector3i(1, 1, 1));
		Vector3i c1_b = (p_c2 + dz).clamp(Vector3i(), res - Vector3i(1, 1, 1));

		Vector3 pre = lattice_get_i_interp_2d(c0_a, c0_b, p_tv.x, p_tv.y);
		Vector3 post = lattice_get_i_interp_2d(c1_a, c1_b, p_tv.x, p_tv.y);

		return a.cubic_interpolate(b, pre, post, p_tv.z);
	}

	void build_lattice() {
		int size = lattice_resolution.x * lattice_resolution.y * lattice_resolution.z;

		lattice.clear();
		lattice.resize(size);

		for (int i = 0; i < size; ++i) {
			lattice.write[i] = Vector3();
		}
	}

	void begin_operation() {
		start_lattice = lattice.duplicate();
	}

	void end_operation() {
		start_lattice.resize(0);
	}

	PackedVector3Array get_lattice_difference() {
		PackedVector3Array ret = start_lattice.duplicate();
		for (int lattice_i = 0; lattice_i < ret.size(); lattice_i++) {
			ret.set(lattice_i, lattice[lattice_i] - ret[lattice_i]);
		}
		return ret;
	}

	void add_lattice(const PackedVector3Array &p_lattice, float p_sign) {
		ERR_FAIL_COND(p_lattice.size() > lattice.size());

		for (int i = 0; i < p_lattice.size(); i++) {
			Vector3 new_lattice = p_lattice[i];
			lattice.set(i, p_lattice[i] + new_lattice * p_sign);
		}
		_mesh_changed();
	}

	void affect_lattice(const Vector3 &p_where, float p_radius, const Vector3 &p_normal, float p_delta, const String &p_mode) {
		CSGBrush *mesh_brush = _get_brush();
		if (!mesh_brush) {
			return;
		}
		HashMap<Vector3i, real_t> hit_lattice_weights;
		HashMap<Vector3i, int> hit_lattice_counts;
		int closest_i = -1;
		float closest_dist = 10000000000.0;
		CSGBrush::Face *closest_face = nullptr; // Store the closest face
		int closest_v = -1; // Store the vertex index of the closest vertex
		for (int f = 0; f < mesh_brush->faces.size(); ++f) {
			CSGBrush::Face &face = mesh_brush->faces.write[f];
			for (int v = 0; v < 3; ++v) {
				Vector3 vert = face.vertices[v];
				Vector3 diff = (vert - p_where) / p_radius;
				float dist = diff.length_squared();
				if (dist < closest_dist) {
					closest_dist = dist;
					closest_i = f * 3 + v;
					closest_face = &face; // Update the closest face
					closest_v = v; // Update the vertex index
				}
				float l = 1.0 - dist;
				l = MAX(0, l);
				l *= l;
				if (l > 0.0) {
					Vector3 original_vert = face.vertices[v];
					lattice_get_weights(original_vert, l, hit_lattice_weights, hit_lattice_counts);
				}
			}
		}

		// ensure the minimum weight is at least 1.0
		float max_weight = 0.0;
		for (KeyValue<Vector3i, real_t> coord : hit_lattice_weights) {
			max_weight = MAX(max_weight, (hit_lattice_weights[coord.key] / hit_lattice_counts[coord.key]));
		}
		if (max_weight > 0.0) {
			max_weight = 1.0 / MIN(max_weight, 1.0);
		} else {
			if (closest_i >= 0 && closest_face != nullptr) { // Check if a closest vertex was found
				Vector3 original_vert = closest_face->vertices[closest_v]; // Use the stored face and vertex index
				lattice_get_weights(original_vert, 1.0, hit_lattice_weights, hit_lattice_counts);

				max_weight = 1.0;
			}
		}

		Vector3 res = lattice_resolution;
		Vector3 avg_deform;
		float avg_deform_weight = 0.0;

		for (KeyValue<Vector3i, real_t> coord : hit_lattice_weights) {
			hit_lattice_weights[coord.key] *= max_weight;
			if (p_mode == "Smooth" || p_mode == "Relax" || p_mode == "Average") {
				float weight = hit_lattice_weights[coord.key];
				avg_deform_weight += weight;
				int index = coord.key.z * res.x * res.y + coord.key.y * res.x + coord.key.x;
				avg_deform += lattice[index] * weight;
			}
		}

		if (avg_deform_weight > 0.0) {
			avg_deform /= avg_deform_weight;
		}

		for (KeyValue<Vector3i, real_t> coord : hit_lattice_weights) {
			float weight = hit_lattice_weights[coord.key] / hit_lattice_counts[coord.key];
			int index = coord.key.z * res.x * res.y + coord.key.y * res.x + coord.key.x;
			if (p_mode == "Grow") {
				lattice.write[index] += p_normal * weight * p_delta;
			} else if (p_mode == "Erase") {
				Vector3 erased = lattice[index].lerp(Vector3(), 1.0 - pow(0.5, abs(p_delta) * weight * 10.0));
				if (p_delta > 0.0) {
					lattice.write[index] = erased;
				} else {
					lattice.write[index] = lattice[index].lerp(erased, -1.0);
				}
			} else if (p_mode == "Smooth" || p_mode == "Relax" || p_mode == "Average") {
				Vector3 smoothed = lattice[index].lerp(avg_deform, 1.0 - pow(0.5, abs(p_delta) * weight * 10.0));
				Vector3 diff = smoothed - lattice[index];
				if (p_mode == "Smooth") {
					diff = diff.project(p_normal);
				} else if (p_mode == "Relax") {
					diff = diff.slide(p_normal);
				}
				if (p_delta > 0.0) {
					lattice.write[index] += diff;
				} else {
					lattice.write[index] -= diff;
				}
			}
		}
		_mesh_changed();
	}

	Vector3 translate_coord(const Vector3 &p_coordinate, bool p_force_fast = false) {
		if (smooth && !p_force_fast) {
			return p_coordinate + lattice_get_smooth(p_coordinate);
		} else {
			return p_coordinate + lattice_get_fast(p_coordinate);
		}
	}

	Vector3 get_normal_at_coordinate(const Vector3 &p_c, const Vector3 &p_n) {
		Vector3 y = p_n.cross(Vector3(p_n.y, p_n.z, -p_n.x));
		Vector3 x = p_n.cross(y);
		Vector3 s = lattice_size / Vector3(lattice_resolution) / 2.0; // * (0.5 if smooth else 1.0)
		Vector3 tan = (translate_coord(p_c + x * s, true) - translate_coord(p_c - x * s, true));
		Vector3 bitan = (translate_coord(p_c + y * s, true) - translate_coord(p_c - y * s, true));
		return -tan.cross(bitan).normalized();
	}

	void randomize_lattice() {
		float s = 0.1;
		Ref<RandomNumberGenerator> random_generator;
		random_generator.instantiate();

		for (int i = 0; i < lattice.size(); ++i) {
			random_generator->set_seed(HashMapHasherDefault::hash(i) + OS::get_singleton()->get_ticks_usec());
			lattice.write[i] = Vector3(random_generator->randf_range(-s, s), random_generator->randf_range(-s, s), random_generator->randf_range(-s, s));
		}
	}

protected:
	static void _bind_methods() {
		ClassDB::bind_method(D_METHOD("set_lattice", "value"), &CSGDeform3D::set_lattice);
		ClassDB::bind_method(D_METHOD("get_lattice"), &CSGDeform3D::get_lattice);

		ADD_PROPERTY(PropertyInfo(Variant::PACKED_VECTOR3_ARRAY, "lattice",PROPERTY_HINT_NONE, "", PROPERTY_USAGE_NO_EDITOR), "set_lattice", "get_lattice");

		ClassDB::bind_method(D_METHOD("set_lattice_size", "value"), &CSGDeform3D::set_lattice_size);
		ClassDB::bind_method(D_METHOD("get_lattice_size"), &CSGDeform3D::get_lattice_size);
		ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "lattice_size"), "set_lattice_size", "get_lattice_size");

		ClassDB::bind_method(D_METHOD("set_lattice_resolution", "value"), &CSGDeform3D::set_lattice_resolution);
		ClassDB::bind_method(D_METHOD("get_lattice_resolution"), &CSGDeform3D::get_lattice_resolution);
		ADD_PROPERTY(PropertyInfo(Variant::VECTOR3I, "lattice_resolution"), "set_lattice_resolution", "get_lattice_resolution");

		ClassDB::bind_method(D_METHOD("set_fix_normals", "value"), &CSGDeform3D::set_fix_normals);
		ClassDB::bind_method(D_METHOD("get_fix_normals"), &CSGDeform3D::get_fix_normals);
		ADD_PROPERTY(PropertyInfo(Variant::BOOL, "fix_normals"), "set_fix_normals", "get_fix_normals");

		ClassDB::bind_method(D_METHOD("set_smooth", "value"), &CSGDeform3D::set_smooth);
		ClassDB::bind_method(D_METHOD("get_smooth"), &CSGDeform3D::get_smooth);
		ADD_PROPERTY(PropertyInfo(Variant::BOOL, "smooth"), "set_smooth", "get_smooth");

		ClassDB::bind_method(D_METHOD("lattice_get_fast", "coord"), &CSGDeform3D::lattice_get_fast);
		ClassDB::bind_method(D_METHOD("lattice_get_smooth", "coord"), &CSGDeform3D::lattice_get_smooth);
		// ClassDB::bind_method(D_METHOD("lattice_get_weights", "coord", "amount", "weights", "counts"), &CSGDeform3D::lattice_get_weights);
		ClassDB::bind_method(D_METHOD("lattice_get_i_interp_1d", "c1", "c2", "t"), &CSGDeform3D::lattice_get_i_interp_1d);
		ClassDB::bind_method(D_METHOD("lattice_get_i_interp_2d", "c1", "c2", "tx", "ty"), &CSGDeform3D::lattice_get_i_interp_2d);
		ClassDB::bind_method(D_METHOD("lattice_get_i_interp_3d", "c1", "c2", "tv"), &CSGDeform3D::lattice_get_i_interp_3d);
		ClassDB::bind_method(D_METHOD("build_lattice"), &CSGDeform3D::build_lattice);
		ClassDB::bind_method(D_METHOD("begin_operation"), &CSGDeform3D::begin_operation);
		ClassDB::bind_method(D_METHOD("end_operation"), &CSGDeform3D::end_operation);
		ClassDB::bind_method(D_METHOD("get_lattice_difference"), &CSGDeform3D::get_lattice_difference);
		ClassDB::bind_method(D_METHOD("add_lattice", "lattice", "sign"), &CSGDeform3D::add_lattice);
		ClassDB::bind_method(D_METHOD("affect_lattice", "where", "radius", "normal", "delta", "mode"), &CSGDeform3D::affect_lattice, DEFVAL(Vector3(0, 0, 0)));
		ClassDB::bind_method(D_METHOD("translate_coord", "coord", "force_fast"), &CSGDeform3D::translate_coord);
		ClassDB::bind_method(D_METHOD("get_normal_at_coordinate", "coordinate", "normal"), &CSGDeform3D::get_normal_at_coordinate);
		ClassDB::bind_method(D_METHOD("randomize_lattice"), &CSGDeform3D::randomize_lattice);
	}

	virtual CSGBrush *_build_brush() override {
		set_use_collision(true);
		int capacity = lattice_resolution.x * lattice_resolution.y * lattice_resolution.z;
		if (lattice.size() != capacity) {
			lattice.resize(capacity);
			build_lattice();
		}
		CSGBrush *mesh_brush = CSGMesh3D::_build_brush();
		for (int f = 0; f < mesh_brush->faces.size(); ++f) {
			CSGBrush::Face &face = mesh_brush->faces.write[f];
			for (int v = 0; v < 3; ++v) {
				face.vertices[v] = translate_coord(face.vertices[v]);
			}
		}
		return mesh_brush;
	}
};

class CSGSphere3D : public CSGPrimitive3D {
	GDCLASS(CSGSphere3D, CSGPrimitive3D);
	virtual CSGBrush *_build_brush() override;

	Ref<Material> material;
	bool smooth_faces;
	float radius;
	int radial_segments;
	int rings;

protected:
	static void _bind_methods();

public:
	void set_radius(const float p_radius);
	float get_radius() const;

	void set_radial_segments(const int p_radial_segments);
	int get_radial_segments() const;

	void set_rings(const int p_rings);
	int get_rings() const;

	void set_material(const Ref<Material> &p_material);
	Ref<Material> get_material() const;

	void set_smooth_faces(bool p_smooth_faces);
	bool get_smooth_faces() const;

	CSGSphere3D();
};

class CSGBox3D : public CSGPrimitive3D {
	GDCLASS(CSGBox3D, CSGPrimitive3D);
	virtual CSGBrush *_build_brush() override;

	Ref<Material> material;
	Vector3 size = Vector3(1, 1, 1);

protected:
	static void _bind_methods();
#ifndef DISABLE_DEPRECATED
	// Kept for compatibility from 3.x to 4.0.
	bool _set(const StringName &p_name, const Variant &p_value);
#endif

public:
	void set_size(const Vector3 &p_size);
	Vector3 get_size() const;

	void set_material(const Ref<Material> &p_material);
	Ref<Material> get_material() const;

	CSGBox3D() {}
};

class CSGCylinder3D : public CSGPrimitive3D {
	GDCLASS(CSGCylinder3D, CSGPrimitive3D);
	virtual CSGBrush *_build_brush() override;

	Ref<Material> material;
	float radius;
	float height;
	int sides;
	bool cone;
	bool smooth_faces;

protected:
	static void _bind_methods();

public:
	void set_radius(const float p_radius);
	float get_radius() const;

	void set_height(const float p_height);
	float get_height() const;

	void set_sides(const int p_sides);
	int get_sides() const;

	void set_cone(const bool p_cone);
	bool is_cone() const;

	void set_smooth_faces(bool p_smooth_faces);
	bool get_smooth_faces() const;

	void set_material(const Ref<Material> &p_material);
	Ref<Material> get_material() const;

	CSGCylinder3D();
};

class CSGTorus3D : public CSGPrimitive3D {
	GDCLASS(CSGTorus3D, CSGPrimitive3D);
	virtual CSGBrush *_build_brush() override;

	Ref<Material> material;
	float inner_radius;
	float outer_radius;
	int sides;
	int ring_sides;
	bool smooth_faces;

protected:
	static void _bind_methods();

public:
	void set_inner_radius(const float p_inner_radius);
	float get_inner_radius() const;

	void set_outer_radius(const float p_outer_radius);
	float get_outer_radius() const;

	void set_sides(const int p_sides);
	int get_sides() const;

	void set_ring_sides(const int p_ring_sides);
	int get_ring_sides() const;

	void set_smooth_faces(bool p_smooth_faces);
	bool get_smooth_faces() const;

	void set_material(const Ref<Material> &p_material);
	Ref<Material> get_material() const;

	CSGTorus3D();
};

class CSGPolygon3D : public CSGPrimitive3D {
	GDCLASS(CSGPolygon3D, CSGPrimitive3D);

public:
	enum Mode {
		MODE_DEPTH,
		MODE_SPIN,
		MODE_PATH
	};

	enum PathIntervalType {
		PATH_INTERVAL_DISTANCE,
		PATH_INTERVAL_SUBDIVIDE
	};

	enum PathRotation {
		PATH_ROTATION_POLYGON,
		PATH_ROTATION_PATH,
		PATH_ROTATION_PATH_FOLLOW,
	};

private:
	virtual CSGBrush *_build_brush() override;

	Vector<Vector2> polygon;
	Ref<Material> material;

	Mode mode;

	float depth;

	float spin_degrees;
	int spin_sides;

	NodePath path_node;
	PathIntervalType path_interval_type;
	float path_interval;
	float path_simplify_angle;
	PathRotation path_rotation;
	bool path_local;

	Path3D *path = nullptr;

	bool smooth_faces;
	bool path_continuous_u;
	real_t path_u_distance;
	bool path_joined;

	bool _is_editable_3d_polygon() const;
	bool _has_editable_3d_polygon_no_depth() const;

	void _path_changed();
	void _path_exited();

protected:
	static void _bind_methods();
	void _validate_property(PropertyInfo &p_property) const;
	void _notification(int p_what);

public:
	void set_polygon(const Vector<Vector2> &p_polygon);
	Vector<Vector2> get_polygon() const;

	void set_mode(Mode p_mode);
	Mode get_mode() const;

	void set_depth(float p_depth);
	float get_depth() const;

	void set_spin_degrees(float p_spin_degrees);
	float get_spin_degrees() const;

	void set_spin_sides(int p_spin_sides);
	int get_spin_sides() const;

	void set_path_node(const NodePath &p_path);
	NodePath get_path_node() const;

	void set_path_interval_type(PathIntervalType p_interval_type);
	PathIntervalType get_path_interval_type() const;

	void set_path_interval(float p_interval);
	float get_path_interval() const;

	void set_path_simplify_angle(float p_angle);
	float get_path_simplify_angle() const;

	void set_path_rotation(PathRotation p_rotation);
	PathRotation get_path_rotation() const;

	void set_path_local(bool p_enable);
	bool is_path_local() const;

	void set_path_continuous_u(bool p_enable);
	bool is_path_continuous_u() const;

	void set_path_u_distance(real_t p_path_u_distance);
	real_t get_path_u_distance() const;

	void set_path_joined(bool p_enable);
	bool is_path_joined() const;

	void set_smooth_faces(bool p_smooth_faces);
	bool get_smooth_faces() const;

	void set_material(const Ref<Material> &p_material);
	Ref<Material> get_material() const;

	CSGPolygon3D();
};

VARIANT_ENUM_CAST(CSGPolygon3D::Mode)
VARIANT_ENUM_CAST(CSGPolygon3D::PathRotation)
VARIANT_ENUM_CAST(CSGPolygon3D::PathIntervalType)

#endif // CSG_SHAPE_H
