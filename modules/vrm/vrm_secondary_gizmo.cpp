/**************************************************************************/
/*  vrm_secondary_gizmo.cpp                                               */
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

#include "vrm_secondary_gizmo.h"
#include "scene/resources/immediate_mesh.h"

void SecondaryGizmo::_bind_methods() {
	ClassDB::bind_method(D_METHOD("ready", "secondary"), &SecondaryGizmo::ready);
	ClassDB::bind_method(D_METHOD("draw"), &SecondaryGizmo::draw);
	ClassDB::bind_method(D_METHOD("draw_spring_bones", "color"), &SecondaryGizmo::draw_spring_bones);
	ClassDB::bind_method(D_METHOD("draw_collider_groups"), &SecondaryGizmo::draw_collider_groups);
	ClassDB::bind_method(D_METHOD("draw_line", "begin_pos", "end_pos", "color"), &SecondaryGizmo::draw_line);
	ClassDB::bind_method(D_METHOD("draw_sphere", "bas", "center", "radius", "color"), &SecondaryGizmo::draw_sphere);
}

void SecondaryGizmo::draw_spring_bones(const Color &color) {
	if (!secondary_node) {
		return;
	}
	set_material_override(m);

	for (int i = 0; i < secondary_node->spring_bones_internal.size(); ++i) {
		Ref<VRMSpringBone> spring_bone = secondary_node->spring_bones_internal[i];
		if (spring_bone->verlets.size() == 0 && spring_bone->verlets.size() == 0) {
			continue;
		}
		Ref<ImmediateMesh>(mesh)->surface_begin(Mesh::PRIMITIVE_LINES);

		for (int j = 0; j < spring_bone->verlets.size(); ++j) {
			Ref<VRMSpringBoneLogic> v = spring_bone->verlets[j];
			Transform3D s_tr = Transform3D();
			Skeleton3D *s_sk = spring_bone->skeleton;

			if (Engine::get_singleton()->is_editor_hint()) {
				if (v->bone_idx != -1) {
					s_tr = s_sk->get_bone_global_pose(v->bone_idx);
				}
			} else {
				s_tr = spring_bone->skeleton->get_bone_global_pose_no_override(v->bone_idx);
			}

			draw_line(s_tr.origin, s_sk->get_global_transform().xform_inv(v->get_current_tail()), color);
		}

		for (int j = 0; j < spring_bone->verlets.size(); ++j) {
			Ref<VRMSpringBoneLogic> v = spring_bone->verlets[j];
			Transform3D s_tr = Transform3D();
			Skeleton3D *s_sk = spring_bone->skeleton;

			if (Engine::get_singleton()->is_editor_hint()) {
				if (v->bone_idx != -1) {
					s_tr = s_sk->get_bone_global_pose(v->bone_idx);
				}
			} else {
				s_tr = spring_bone->skeleton->get_bone_global_pose_no_override(v->bone_idx);
			}

			draw_sphere(s_tr.basis, s_sk->get_global_transform().xform_inv(v->get_current_tail()), spring_bone->get_hit_radius(), color);
		}

		Ref<ImmediateMesh>(mesh)->surface_end();
	}
}

void SecondaryGizmo::draw_collider_groups() {
	if (!secondary_node) {
		return;
	}
	set_material_override(m);

	TypedArray<VRMColliderGroup> collider_groups = Engine::get_singleton()->is_editor_hint() ? secondary_node->get_collider_groups() : secondary_node->get_collider_groups_internal();

	for (int collider_group_i = 0; collider_group_i < collider_groups.size(); ++collider_group_i) {
		Ref<VRMColliderGroup> collider_group = collider_groups[collider_group_i];
		Ref<ImmediateMesh>(mesh)->surface_begin(Mesh::PRIMITIVE_LINE_STRIP);

		Transform3D c_tr = Transform3D();

		if (Engine::get_singleton()->is_editor_hint()) {
			Node *c_node = secondary_node->get_node_or_null(collider_group->skeleton_or_node);
			Skeleton3D *c_sk = Object::cast_to<Skeleton3D>(c_node);
			if (c_sk) {
				if (collider_group->get_bone_idx() == -1) {
					collider_group->set_bone_idx(c_sk->find_bone(collider_group->bone));
				}
				c_tr = c_sk->get_bone_global_pose(collider_group->get_bone_idx());
			}
		} else if (Object::cast_to<Skeleton3D>(collider_group->get_skel())) {
			if (collider_group->get_skel()) {
				c_tr = collider_group->get_skel()->get_bone_global_pose_no_override(collider_group->get_skel()->find_bone(collider_group->get_bone()));
			}
		}
		Array sphere_colliders = collider_group->get_sphere_colliders();
		for (int j = 0; j < sphere_colliders.size(); ++j) {
			Plane collider = sphere_colliders[j];
			Vector3 c_ps = Vector3(collider.normal.x, collider.normal.y, collider.normal.z);
			draw_sphere(c_tr.basis, c_tr.xform(c_ps), collider.d, collider_group->get_gizmo_color());
		}

		Ref<ImmediateMesh>(mesh)->surface_end();
	}
}

void SecondaryGizmo::draw_line(Vector3 begin_pos, Vector3 end_pos, Color color) {
	Ref<ImmediateMesh> immediate_mesh = Ref<ImmediateMesh>(mesh);
	if (immediate_mesh.is_null()) {
		return;
	}
	immediate_mesh->surface_set_color(color);
	immediate_mesh->surface_add_vertex(begin_pos);
	immediate_mesh->surface_set_color(color);
	immediate_mesh->surface_add_vertex(end_pos);
}

void SecondaryGizmo::draw_sphere(Basis bas, Vector3 center, float radius, const Color color) {
	Ref<ImmediateMesh> immediate_mesh = Ref<ImmediateMesh>(mesh);
	if (immediate_mesh.is_null()) {
		return;
	}
	int step = 15;
	float sppi = 2 * Math_PI / step;
	for (int i = 0; i <= step; ++i) {
		immediate_mesh->surface_set_color(color);
		immediate_mesh->surface_add_vertex(center + ((bas.xform(Vector3(0, radius, 0))).rotated(bas.xform(Vector3(1, 0, 0)), sppi * (i % step))));
	}

	for (int i = 0; i <= step; ++i) {
		immediate_mesh->surface_set_color(color);
		immediate_mesh->surface_add_vertex(center + ((bas.xform(Vector3(radius, 0, 0))).rotated(bas.xform(Vector3(0, 0, 1)), sppi * (i % step))));
	}

	for (int i = 0; i <= step; ++i) {
		immediate_mesh->surface_set_color(color);
		immediate_mesh->surface_add_vertex(center + ((bas.xform(Vector3(0, 0, radius))).rotated(bas.xform(Vector3(0, 1, 0)), sppi * (i % step))));
	}
}

void SecondaryGizmo::draw() {
	Ref<ImmediateMesh> immediate_mesh = Ref<ImmediateMesh>(mesh);
	if (immediate_mesh.is_null()) {
		return;
	}
	immediate_mesh->clear_surfaces();
	if (secondary_node && Object::cast_to<VRMTopLevel>(secondary_node->get_parent())) {
		draw_spring_bones(Object::cast_to<VRMTopLevel>(secondary_node->get_parent())->get_gizmo_spring_bone_color());
		draw_collider_groups();
	}
}

void SecondaryGizmo::ready(Node *p_secondary_node) {
	Ref<ImmediateMesh> immediate_mesh;
	immediate_mesh.instantiate();
	set_mesh(immediate_mesh);
	secondary_node = cast_to<VRMSecondary>(p_secondary_node);
	m.instantiate();
	m->set_depth_draw_mode(BaseMaterial3D::DEPTH_DRAW_DISABLED);
	m->set_shading_mode(BaseMaterial3D::SHADING_MODE_UNSHADED);
	m->set_flag(StandardMaterial3D::FLAG_ALBEDO_FROM_VERTEX_COLOR, true);
	m->set_transparency(BaseMaterial3D::TRANSPARENCY_ALPHA);
}
