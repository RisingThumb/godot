/**************************************************************************/
/*  many_bone_ik_3d.cpp                                                   */
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

#include "many_bone_ik_3d.h"
#include "core/core_string_names.h"
#include "core/error/error_macros.h"
#include "core/io/json.h"
#include "core/object/class_db.h"
#include "core/string/string_name.h"
#include "core/variant/typed_array.h"
#include "ik_bone_3d.h"
#include "ik_kusudama_3d.h"
#include "ik_limit_cone_3d.h"
#include "scene/3d/marker_3d.h"
#include "scene/3d/physics_body_3d.h"
#include "scene/3d/skeleton_3d.h"
#include "scene/resources/skeleton_profile.h"

#ifdef TOOLS_ENABLED
#include "editor/editor_node.h"
#endif

void ManyBoneIK3D::set_pin_count(int32_t p_value) {
	int32_t old_count = pins.size();
	pin_count = p_value;
	pins.resize(p_value);
	for (int32_t pin_i = p_value; pin_i-- > old_count;) {
		pins.write[pin_i].instantiate();
	}
	set_dirty();
}

int32_t ManyBoneIK3D::get_pin_count() const {
	return pin_count;
}

void ManyBoneIK3D::set_pin_bone(int32_t p_pin_index, const String &p_bone) {
	ERR_FAIL_INDEX(p_pin_index, pins.size());
	Ref<IKEffectorTemplate3D> effector_template = pins[p_pin_index];
	if (effector_template.is_null()) {
		effector_template.instantiate();
		pins.write[p_pin_index] = effector_template;
	}
	effector_template->set_name(p_bone);
	set_dirty();
}

void ManyBoneIK3D::set_pin_target_nodepath(int32_t p_pin_index, const NodePath &p_target_node) {
	ERR_FAIL_INDEX(p_pin_index, pins.size());
	Ref<IKEffectorTemplate3D> effector_template = pins[p_pin_index];
	if (effector_template.is_null()) {
		effector_template.instantiate();
		pins.write[p_pin_index] = effector_template;
	}
	effector_template->set_target_node(p_target_node);
	set_dirty();
}

NodePath ManyBoneIK3D::get_pin_target_nodepath(int32_t p_pin_index) {
	ERR_FAIL_INDEX_V(p_pin_index, pins.size(), NodePath());
	const Ref<IKEffectorTemplate3D> effector_template = pins[p_pin_index];
	return effector_template->get_target_node();
}

Vector<Ref<IKEffectorTemplate3D>> ManyBoneIK3D::get_bone_effectors() const {
	return pins;
}

void ManyBoneIK3D::_remove_pin(int32_t p_index) {
	ERR_FAIL_INDEX(p_index, pins.size());
	pins.remove_at(p_index);
	pin_count--;
	pins.resize(pin_count);
	set_dirty();
}

void ManyBoneIK3D::update_ik_bones_transform() {
	for (int32_t bone_i = bone_list.size(); bone_i-- > 0;) {
		Ref<IKBone3D> bone = bone_list[bone_i];
		if (bone.is_null()) {
			continue;
		}
		bone->set_initial_pose(get_skeleton());
		if (bone->is_pinned()) {
			bone->get_pin()->update_target_global_transform(get_skeleton(), this);
		}
	}
}

void ManyBoneIK3D::update_skeleton_bones_transform() {
	for (int32_t bone_i = bone_list.size(); bone_i-- > 0;) {
		Ref<IKBone3D> bone = bone_list[bone_i];
		if (bone.is_null()) {
			continue;
		}
		if (bone->get_bone_id() == -1) {
			continue;
		}
		bone->set_skeleton_bone_pose(get_skeleton());
	}
}

void ManyBoneIK3D::_get_property_list(List<PropertyInfo> *p_list) const {
	RBSet<StringName> existing_pins;
	for (int32_t pin_i = 0; pin_i < get_pin_count(); pin_i++) {
		const String name = get_pin_bone_name(pin_i);
		existing_pins.insert(name);
	}
	p_list->push_back(
			PropertyInfo(Variant::INT, "pin_count",
					PROPERTY_HINT_RANGE, "0,65536,or_greater", PROPERTY_USAGE_DEFAULT | PROPERTY_USAGE_ARRAY | PROPERTY_USAGE_READ_ONLY,
					"Pins,pins/"));
	for (int pin_i = 0; pin_i < get_pin_count(); pin_i++) {
		PropertyInfo effector_name;
		effector_name.type = Variant::STRING_NAME;
		effector_name.name = "pins/" + itos(pin_i) + "/bone_name";
		const uint32_t pin_usage = PROPERTY_USAGE_DEFAULT;
		effector_name.usage = pin_usage | PROPERTY_USAGE_READ_ONLY;
		if (get_skeleton()) {
			String names;
			for (int bone_i = 0; bone_i < get_skeleton()->get_bone_count(); bone_i++) {
				String name = get_skeleton()->get_bone_name(bone_i);
				StringName string_name = StringName(name);
				if (existing_pins.has(string_name)) {
					continue;
				}
			}
			effector_name.hint = PROPERTY_HINT_ENUM_SUGGESTION;
			effector_name.hint_string = names;
		} else {
			effector_name.hint = PROPERTY_HINT_NONE;
			effector_name.hint_string = "";
		}
		p_list->push_back(effector_name);
		p_list->push_back(
				PropertyInfo(Variant::NODE_PATH, "pins/" + itos(pin_i) + "/target_node", PROPERTY_HINT_NODE_PATH_VALID_TYPES, "Node3D", pin_usage));
		p_list->push_back(
				PropertyInfo(Variant::FLOAT, "pins/" + itos(pin_i) + "/passthrough_factor", PROPERTY_HINT_RANGE, "0,1,0.1,or_greater", pin_usage));
		p_list->push_back(
				PropertyInfo(Variant::FLOAT, "pins/" + itos(pin_i) + "/weight", PROPERTY_HINT_RANGE, "0,1,0.1,or_greater", pin_usage));
		p_list->push_back(
				PropertyInfo(Variant::VECTOR3, "pins/" + itos(pin_i) + "/direction_priorities", PROPERTY_HINT_RANGE, "0,1,0.1,or_greater", pin_usage));
	}

	const Vector<Ref<IKBone3D>> ik_bones = get_bone_list();

	RBSet<String> existing_constraints;
	for (int32_t constraint_i = 0; constraint_i < ik_bones.size(); constraint_i++) {
		const String name = ik_bones[constraint_i]->get_name();
		existing_constraints.insert(name);
	}
	p_list->push_back(
			PropertyInfo(Variant::INT, "constraint_count",
					PROPERTY_HINT_RANGE, "0,256,or_greater", PROPERTY_USAGE_DEFAULT | PROPERTY_USAGE_ARRAY | PROPERTY_USAGE_READ_ONLY,
					"Kusudama Constraints,constraints/"));
	for (int constraint_i = 0; constraint_i < get_constraint_count(); constraint_i++) {
		PropertyInfo bone_name;
		bone_name.type = Variant::STRING_NAME;
		const uint32_t constraint_usage = PROPERTY_USAGE_DEFAULT;
		bone_name.usage = constraint_usage;
		bone_name.name = "constraints/" + itos(constraint_i) + "/bone_name";
		if (get_skeleton()) {
			String names;
			for (int bone_i = 0; bone_i < get_constraint_count(); bone_i++) {
				String name = get_constraint_name(bone_i);
				if (existing_constraints.has(name)) {
					continue;
				}
				name += ",";
				names += name;
				existing_constraints.insert(name);
			}
			bone_name.hint = PROPERTY_HINT_ENUM_SUGGESTION;
			bone_name.hint_string = names;
		} else {
			bone_name.hint = PROPERTY_HINT_NONE;
			bone_name.hint_string = "";
		}
		p_list->push_back(bone_name);
		p_list->push_back(
				PropertyInfo(Variant::FLOAT, "constraints/" + itos(constraint_i) + "/painfulness", PROPERTY_HINT_RANGE, "0,1,0.01,exp", constraint_usage));
		p_list->push_back(
				PropertyInfo(Variant::FLOAT, "constraints/" + itos(constraint_i) + "/twist_from", PROPERTY_HINT_RANGE, "-359.9,359.9,0.1,radians,exp", constraint_usage));
		p_list->push_back(
				PropertyInfo(Variant::FLOAT, "constraints/" + itos(constraint_i) + "/twist_range", PROPERTY_HINT_RANGE, "-359.9,359.9,0.1,radians,exp", constraint_usage));
		p_list->push_back(
				PropertyInfo(Variant::FLOAT, "constraints/" + itos(constraint_i) + "/twist_current", PROPERTY_HINT_RANGE, "0,1,0.1,or_greater", constraint_usage));
		p_list->push_back(
				PropertyInfo(Variant::INT, "constraints/" + itos(constraint_i) + "/kusudama_limit_cone_count", PROPERTY_HINT_RANGE, "0,10,1", constraint_usage | PROPERTY_USAGE_ARRAY | PROPERTY_USAGE_READ_ONLY,
						"Limit Cones,constraints/" + itos(constraint_i) + "/kusudama_limit_cone/"));
		for (int cone_i = 0; cone_i < get_kusudama_limit_cone_count(constraint_i); cone_i++) {
			p_list->push_back(
					PropertyInfo(Variant::VECTOR3, "constraints/" + itos(constraint_i) + "/kusudama_limit_cone/" + itos(cone_i) + "/center", PROPERTY_HINT_RANGE, "-1.0,1.0,0.01,or_greater,exp", constraint_usage));
			p_list->push_back(
					PropertyInfo(Variant::FLOAT, "constraints/" + itos(constraint_i) + "/kusudama_limit_cone/" + itos(cone_i) + "/radius", PROPERTY_HINT_RANGE, "0,180,0.1,radians,exp", constraint_usage));
		}
		p_list->push_back(
				PropertyInfo(Variant::TRANSFORM3D, "constraints/" + itos(constraint_i) + "/kusudama_twist", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_NO_EDITOR));
		p_list->push_back(
				PropertyInfo(Variant::TRANSFORM3D, "constraints/" + itos(constraint_i) + "/kusudama_orientation", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_NO_EDITOR));
		p_list->push_back(
				PropertyInfo(Variant::TRANSFORM3D, "constraints/" + itos(constraint_i) + "/bone_direction", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_NO_EDITOR));
	}
}

bool ManyBoneIK3D::_get(const StringName &p_name, Variant &r_ret) const {
	String name = p_name;
	if (name == "constraint_count") {
		r_ret = get_constraint_count();
		return true;
	} else if (name == "pin_count") {
		r_ret = get_pin_count();
		return true;
	} else if (name == "bone_count") {
		r_ret = get_bone_count();
		return true;
	} else if (name.begins_with("pins/")) {
		int index = name.get_slicec('/', 1).to_int();
		String what = name.get_slicec('/', 2);
		ERR_FAIL_INDEX_V(index, pins.size(), false);
		Ref<IKEffectorTemplate3D> effector_template = pins[index];
		ERR_FAIL_NULL_V(effector_template, false);
		if (what == "bone_name") {
			r_ret = effector_template->get_name();
			return true;
		} else if (what == "target_node") {
			r_ret = effector_template->get_target_node();
			return true;
		} else if (what == "passthrough_factor") {
			r_ret = get_pin_passthrough_factor(index);
			return true;
		} else if (what == "weight") {
			r_ret = get_pin_weight(index);
			return true;
		} else if (what == "direction_priorities") {
			r_ret = get_pin_direction_priorities(index);
			return true;
		}
	} else if (name.begins_with("constraints/")) {
		int index = name.get_slicec('/', 1).to_int();
		String what = name.get_slicec('/', 2);
		ERR_FAIL_INDEX_V(index, constraint_count, false);
		String begins = "constraints/" + itos(index) + "/kusudama_limit_cone";
		if (what == "bone_name") {
			ERR_FAIL_INDEX_V(index, constraint_names.size(), false);
			r_ret = constraint_names[index];
			return true;
		} else if (what == "painfulness") {
			r_ret = get_kusudama_painfulness(index);
			return true;
		} else if (what == "twist_current") {
			r_ret = get_kusudama_twist_current(index);
			return true;
		} else if (what == "twist_from") {
			r_ret = get_kusudama_twist(index).x;
			return true;
		} else if (what == "twist_range") {
			r_ret = get_kusudama_twist(index).y;
			return true;
		} else if (what == "kusudama_limit_cone_count") {
			r_ret = get_kusudama_limit_cone_count(index);
			return true;
		} else if (name.begins_with(begins)) {
			int32_t cone_index = name.get_slicec('/', 3).to_int();
			String cone_what = name.get_slicec('/', 4);
			if (cone_what == "center") {
				r_ret = get_kusudama_limit_cone_center(index, cone_index);
				return true;
			} else if (cone_what == "radius") {
				r_ret = get_kusudama_limit_cone_radius(index, cone_index);
				return true;
			}
		} else if (what == "bone_direction") {
			r_ret = get_bone_direction_transform(index);
			return true;
		} else if (what == "kusudama_orientation") {
			r_ret = get_constraint_orientation_transform(index);
			return true;
		} else if (what == "kusudama_twist") {
			r_ret = get_constraint_twist_transform(index);
			return true;
		}
	}
	return false;
}

bool ManyBoneIK3D::_set(const StringName &p_name, const Variant &p_value) {
	String name = p_name;
	if (name == "constraint_count") {
		set_constraint_count(p_value);
		return true;
	} else if (name == "pin_count") {
		set_pin_count(p_value);
		return true;
	} else if (name.begins_with("pins/")) {
		int index = name.get_slicec('/', 1).to_int();
		String what = name.get_slicec('/', 2);
		if (index >= pins.size()) {
			set_pin_count(constraint_count);
		}
		if (what == "bone_name") {
			set_pin_bone(index, p_value);
			return true;
		} else if (what == "target_node") {
			set_pin_target_nodepath(index, p_value);
			String existing_bone = get_pin_bone_name(index);
			if (existing_bone.is_empty()) {
				return false;
			}
			return true;
		} else if (what == "passthrough_factor") {
			set_pin_passthrough_factor(index, p_value);
			return true;
		} else if (what == "weight") {
			set_pin_weight(index, p_value);
			return true;
		} else if (what == "direction_priorities") {
			set_pin_direction_priorities(index, p_value);
			return true;
		}
	} else if (name.begins_with("constraints/")) {
		int index = name.get_slicec('/', 1).to_int();
		String what = name.get_slicec('/', 2);
		String begins = "constraints/" + itos(index) + "/kusudama_limit_cone/";
		if (index >= constraint_names.size()) {
			set_constraint_count(constraint_count);
		}
		if (what == "bone_name") {
			set_constraint_name(index, p_value);
			return true;
		} else if (what == "painfulness") {
			set_kusudama_painfulness(index, p_value);
			return true;
		} else if (what == "twist_current") {
			set_kusudama_twist_current(index, p_value);
			return true;
		} else if (what == "twist_from") {
			Vector2 twist_from = get_kusudama_twist(index);
			set_kusudama_twist(index, Vector2(p_value, twist_from.y));
			return true;
		} else if (what == "twist_range") {
			Vector2 twist_range = get_kusudama_twist(index);
			set_kusudama_twist(index, Vector2(twist_range.x, p_value));
			return true;
		} else if (what == "kusudama_limit_cone_count") {
			set_kusudama_limit_cone_count(index, p_value);
			return true;
		} else if (name.begins_with(begins)) {
			int cone_index = name.get_slicec('/', 3).to_int();
			String cone_what = name.get_slicec('/', 4);
			if (cone_what == "center") {
				Vector3 center = p_value;
				if (Math::is_zero_approx(center.length_squared())) {
					center = Vector3(0.0, 1.0, 0.0);
				}
				set_kusudama_limit_cone_center(index, cone_index, center);
				return true;
			} else if (cone_what == "radius") {
				set_kusudama_limit_cone_radius(index, cone_index, p_value);
				return true;
			}
		} else if (what == "bone_direction") {
			set_bone_direction_transform(index, p_value);
			return true;
		} else if (what == "kusudama_orientation") {
			set_constraint_orientation_transform(index, p_value);
			return true;
		} else if (what == "kusudama_twist") {
			set_constraint_twist_transform(index, p_value);
			return true;
		}
	}

	return false;
}

void ManyBoneIK3D::_bind_methods() {
	ClassDB::bind_method(D_METHOD("get_constraint_twist_transform", "index"), &ManyBoneIK3D::get_constraint_twist_transform);
	ClassDB::bind_method(D_METHOD("set_constraint_twist_transform", "index", "transform"), &ManyBoneIK3D::set_constraint_twist_transform);
	ClassDB::bind_method(D_METHOD("get_constraint_orientation_transform", "index"), &ManyBoneIK3D::get_constraint_orientation_transform);
	ClassDB::bind_method(D_METHOD("set_constraint_orientation_transform", "index", "transform"), &ManyBoneIK3D::set_constraint_orientation_transform);
	ClassDB::bind_method(D_METHOD("get_bone_direction_transform", "index"), &ManyBoneIK3D::get_bone_direction_transform);
	ClassDB::bind_method(D_METHOD("set_bone_direction_transform", "index", "transform"), &ManyBoneIK3D::set_bone_direction_transform);
	ClassDB::bind_method(D_METHOD("get_pin_enabled", "index"), &ManyBoneIK3D::get_pin_enabled);
	ClassDB::bind_method(D_METHOD("remove_constraint", "index"), &ManyBoneIK3D::remove_constraint);
	ClassDB::bind_method(D_METHOD("set_skeleton_node_path", "path"), &ManyBoneIK3D::set_skeleton_node_path);
	ClassDB::bind_method(D_METHOD("get_skeleton_node_path"), &ManyBoneIK3D::get_skeleton_node_path);
	ClassDB::bind_method(D_METHOD("register_skeleton"), &ManyBoneIK3D::register_skeleton);
	ClassDB::bind_method(D_METHOD("reset_constraints"), &ManyBoneIK3D::register_skeleton);
	ClassDB::bind_method(D_METHOD("set_pin_weight", "index", "weight"), &ManyBoneIK3D::set_pin_weight);
	ClassDB::bind_method(D_METHOD("get_pin_weight", "index"), &ManyBoneIK3D::get_pin_weight);
	ClassDB::bind_method(D_METHOD("set_dirty"), &ManyBoneIK3D::set_dirty);
	ClassDB::bind_method(D_METHOD("set_kusudama_limit_cone_radius", "index", "cone_index", "radius"), &ManyBoneIK3D::set_kusudama_limit_cone_radius);
	ClassDB::bind_method(D_METHOD("get_kusudama_limit_cone_radius", "index", "cone_index"), &ManyBoneIK3D::get_kusudama_limit_cone_radius);
	ClassDB::bind_method(D_METHOD("set_kusudama_limit_cone_center", "index", "cone_index", "center"), &ManyBoneIK3D::set_kusudama_limit_cone_center);
	ClassDB::bind_method(D_METHOD("get_kusudama_limit_cone_center", "index", "cone_index"), &ManyBoneIK3D::get_kusudama_limit_cone_center);
	ClassDB::bind_method(D_METHOD("set_kusudama_limit_cone_count", "index", "count"), &ManyBoneIK3D::set_kusudama_limit_cone_count);
	ClassDB::bind_method(D_METHOD("get_kusudama_limit_cone_count", "index"), &ManyBoneIK3D::get_kusudama_limit_cone_count);
	ClassDB::bind_method(D_METHOD("set_kusudama_twist", "index", "limit"), &ManyBoneIK3D::set_kusudama_twist);
	ClassDB::bind_method(D_METHOD("get_kusudama_twist", "index"), &ManyBoneIK3D::get_kusudama_twist);
	ClassDB::bind_method(D_METHOD("set_pin_passthrough_factor", "index", "falloff"), &ManyBoneIK3D::set_pin_passthrough_factor);
	ClassDB::bind_method(D_METHOD("get_pin_passthrough_factor", "index"), &ManyBoneIK3D::get_pin_passthrough_factor);
	ClassDB::bind_method(D_METHOD("get_constraint_name", "index"), &ManyBoneIK3D::get_constraint_name);
	ClassDB::bind_method(D_METHOD("get_iterations_per_frame"), &ManyBoneIK3D::get_iterations_per_frame);
	ClassDB::bind_method(D_METHOD("set_iterations_per_frame", "count"), &ManyBoneIK3D::set_iterations_per_frame);
	ClassDB::bind_method(D_METHOD("find_constraint", "name"), &ManyBoneIK3D::find_constraint);
	ClassDB::bind_method(D_METHOD("get_constraint_count"), &ManyBoneIK3D::get_constraint_count);
	ClassDB::bind_method(D_METHOD("get_pin_count"), &ManyBoneIK3D::get_pin_count);
	ClassDB::bind_method(D_METHOD("get_pin_bone_name", "index"), &ManyBoneIK3D::get_pin_bone_name);
	ClassDB::bind_method(D_METHOD("get_pin_direction_priorities", "index"), &ManyBoneIK3D::get_pin_direction_priorities);
	ClassDB::bind_method(D_METHOD("set_pin_direction_priorities", "index", "priority"), &ManyBoneIK3D::set_pin_direction_priorities);
	ClassDB::bind_method(D_METHOD("queue_print_skeleton"), &ManyBoneIK3D::queue_print_skeleton);
	ClassDB::bind_method(D_METHOD("get_default_damp"), &ManyBoneIK3D::get_default_damp);
	ClassDB::bind_method(D_METHOD("set_default_damp", "damp"), &ManyBoneIK3D::set_default_damp);
	ClassDB::bind_method(D_METHOD("get_pin_nodepath", "index"), &ManyBoneIK3D::get_pin_nodepath);
	ClassDB::bind_method(D_METHOD("set_pin_nodepath", "index", "nodepath"), &ManyBoneIK3D::set_pin_nodepath);
	ClassDB::bind_method(D_METHOD("get_bone_count"), &ManyBoneIK3D::get_bone_count);
	ClassDB::bind_method(D_METHOD("set_constraint_mode", "enabled"), &ManyBoneIK3D::set_constraint_mode);
	ClassDB::bind_method(D_METHOD("get_constraint_mode"), &ManyBoneIK3D::get_constraint_mode);
	ClassDB::bind_method(D_METHOD("set_ui_selected_bone", "bone"), &ManyBoneIK3D::set_ui_selected_bone);
	ClassDB::bind_method(D_METHOD("get_ui_selected_bone"), &ManyBoneIK3D::get_ui_selected_bone);
	ClassDB::bind_method(D_METHOD("set_twist_constraint_defaults", "defaults"), &ManyBoneIK3D::set_twist_constraint_defaults);
	ClassDB::bind_method(D_METHOD("get_twist_constraint_defaults"), &ManyBoneIK3D::get_twist_constraint_defaults);
	ClassDB::bind_method(D_METHOD("set_orientation_constraint_defaults", "defaults"), &ManyBoneIK3D::set_orientation_constraint_defaults);
	ClassDB::bind_method(D_METHOD("get_orientation_constraint_defaults"), &ManyBoneIK3D::get_orientation_constraint_defaults);
	ClassDB::bind_method(D_METHOD("set_bone_direction_constraint_defaults", "defaults"), &ManyBoneIK3D::set_bone_direction_constraint_defaults);
	ClassDB::bind_method(D_METHOD("get_bone_direction_constraint_defaults"), &ManyBoneIK3D::get_bone_direction_constraint_defaults);
	ClassDB::bind_method(D_METHOD("set_stabilization_passes", "passes"), &ManyBoneIK3D::set_stabilization_passes);
	ClassDB::bind_method(D_METHOD("get_stabilization_passes"), &ManyBoneIK3D::get_stabilization_passes);

	ClassDB::bind_method(D_METHOD("setup_humanoid_bones", "enable"), &ManyBoneIK3D::setup_humanoid_bones);

	ClassDB::bind_method(D_METHOD("set_setup_humanoid_bones", "set_targets"), &ManyBoneIK3D::set_setup_humanoid_bones);
	ClassDB::bind_method(D_METHOD("get_setup_humanoid_bones"), &ManyBoneIK3D::get_setup_humanoid_bones);

	ADD_PROPERTY(PropertyInfo(Variant::BOOL, "initialize_humanoid_bones"), "set_setup_humanoid_bones", "get_setup_humanoid_bones");

	ADD_PROPERTY(PropertyInfo(Variant::INT, "humanoid_mode", PROPERTY_HINT_ENUM, "All,Humanoid,Body"), "set_humanoid_mode", "get_humanoid_mode");
	ADD_PROPERTY(PropertyInfo(Variant::NODE_PATH, "skeleton_node_path"), "set_skeleton_node_path", "get_skeleton_node_path");
	ADD_PROPERTY(PropertyInfo(Variant::INT, "iterations_per_frame", PROPERTY_HINT_RANGE, "1,150,1,or_greater"), "set_iterations_per_frame", "get_iterations_per_frame");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "default_damp", PROPERTY_HINT_RANGE, "0.01,180.0,0.1,radians,exp", PROPERTY_USAGE_DEFAULT | PROPERTY_USAGE_UPDATE_ALL_IF_MODIFIED), "set_default_damp", "get_default_damp");
	ADD_PROPERTY(PropertyInfo(Variant::BOOL, "constraint_mode"), "set_constraint_mode", "get_constraint_mode");
	ADD_PROPERTY(PropertyInfo(Variant::INT, "ui_selected_bone", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_NO_EDITOR), "set_ui_selected_bone", "get_ui_selected_bone");
	ADD_PROPERTY(PropertyInfo(Variant::INT, "stabilization_passes"), "set_stabilization_passes", "get_stabilization_passes");
	ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "twist_constraint_defaults", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_NO_EDITOR), "set_twist_constraint_defaults", "get_twist_constraint_defaults");
	ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "orientation_constraint_defaults", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_NO_EDITOR), "set_orientation_constraint_defaults", "get_orientation_constraint_defaults");
	ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "bone_direction_constraint_defaults", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_NO_EDITOR), "set_bone_direction_constraint_defaults", "get_bone_direction_constraint_defaults");
}

ManyBoneIK3D::ManyBoneIK3D() {
}

ManyBoneIK3D::~ManyBoneIK3D() {
}

void ManyBoneIK3D::queue_print_skeleton() {
	queue_debug_skeleton = true;
}

float ManyBoneIK3D::get_pin_passthrough_factor(int32_t p_effector_index) const {
	ERR_FAIL_INDEX_V(p_effector_index, pins.size(), 0.0f);
	const Ref<IKEffectorTemplate3D> effector_template = pins[p_effector_index];
	return effector_template->get_passthrough_factor();
}

void ManyBoneIK3D::set_pin_passthrough_factor(int32_t p_effector_index, const float p_passthrough_factor) {
	ERR_FAIL_INDEX(p_effector_index, pins.size());
	Ref<IKEffectorTemplate3D> effector_template = pins[p_effector_index];
	ERR_FAIL_NULL(effector_template);
	effector_template->set_passthrough_factor(p_passthrough_factor);
	set_dirty();
}

void ManyBoneIK3D::set_constraint_count(int32_t p_count) {
	int32_t old_count = constraint_names.size();
	constraint_count = p_count;
	constraint_names.resize(p_count);
	kusudama_twist.resize(p_count);
	kusudama_limit_cone_count.resize(p_count);
	kusudama_limit_cones.resize(p_count);
	bone_painfulness.resize(p_count);
	for (int32_t constraint_i = p_count; constraint_i-- > old_count;) {
		constraint_names.write[constraint_i] = String();
		kusudama_limit_cone_count.write[constraint_i] = 0;
		kusudama_limit_cones.write[constraint_i].resize(1);
		kusudama_limit_cones.write[constraint_i].write[0] = Vector4(0, 1, 0, Math_PI);
		kusudama_twist.write[constraint_i] = Vector2(0, Math_TAU - CMP_EPSILON);
		bone_painfulness.write[constraint_i] = 0.0f;
	}
	set_dirty();
}

int32_t ManyBoneIK3D::get_constraint_count() const {
	return constraint_count;
}

inline StringName ManyBoneIK3D::get_constraint_name(int32_t p_index) const {
	ERR_FAIL_INDEX_V(p_index, constraint_names.size(), StringName());
	return constraint_names[p_index];
}

void ManyBoneIK3D::set_kusudama_twist(int32_t p_index, Vector2 p_to) {
	ERR_FAIL_INDEX(p_index, constraint_count);
	kusudama_twist.write[p_index] = p_to;
	set_dirty();
}

int32_t ManyBoneIK3D::find_effector_id(StringName p_bone_name) {
	for (int32_t constraint_i = 0; constraint_i < constraint_count; constraint_i++) {
		if (constraint_names[constraint_i] == p_bone_name) {
			return constraint_i;
		}
	}
	return -1;
}

void ManyBoneIK3D::set_kusudama_limit_cone(int32_t p_constraint_index, int32_t p_index,
		Vector3 p_center, float p_radius) {
	ERR_FAIL_INDEX(p_constraint_index, kusudama_limit_cones.size());
	Vector<Vector4> cones = kusudama_limit_cones.write[p_constraint_index];
	if (Math::is_zero_approx(p_center.length_squared())) {
		p_center = Vector3(0.0f, 1.0f, 0.0f);
	}
	Vector3 center = p_center.normalized();
	Vector4 cone;
	cone.x = center.x;
	cone.y = center.y;
	cone.z = center.z;
	cone.w = p_radius;
	cones.write[p_index] = cone;
	kusudama_limit_cones.write[p_constraint_index] = cones;
	set_dirty();
}

Vector3 ManyBoneIK3D::get_kusudama_limit_cone_center(int32_t p_constraint_index, int32_t p_index) const {
	if (unlikely((p_constraint_index) < 0 || (p_constraint_index) >= (kusudama_limit_cone_count.size()))) {
		ERR_PRINT_ONCE("Can't get limit cone center.");
		return Vector3(0.0, 1.0, 0.0);
	}
	if (unlikely((p_constraint_index) < 0 || (p_constraint_index) >= (kusudama_limit_cones.size()))) {
		ERR_PRINT_ONCE("Can't get limit cone center.");
		return Vector3(0.0, 1.0, 0.0);
	}
	if (unlikely((p_index) < 0 || (p_index) >= (kusudama_limit_cones[p_constraint_index].size()))) {
		ERR_PRINT_ONCE("Can't get limit cone center.");
		return Vector3(0.0, 1.0, 0.0);
	}
	const Vector4 &cone = kusudama_limit_cones[p_constraint_index][p_index];
	Vector3 ret;
	ret.x = cone.x;
	ret.y = cone.y;
	ret.z = cone.z;
	return ret;
}

float ManyBoneIK3D::get_kusudama_limit_cone_radius(int32_t p_constraint_index, int32_t p_index) const {
	ERR_FAIL_INDEX_V(p_constraint_index, kusudama_limit_cone_count.size(), Math_TAU);
	ERR_FAIL_INDEX_V(p_constraint_index, kusudama_limit_cones.size(), Math_TAU);
	ERR_FAIL_INDEX_V(p_index, kusudama_limit_cones[p_constraint_index].size(), Math_TAU);
	return kusudama_limit_cones[p_constraint_index][p_index].w;
}

int32_t ManyBoneIK3D::get_kusudama_limit_cone_count(int32_t p_constraint_index) const {
	ERR_FAIL_INDEX_V(p_constraint_index, kusudama_limit_cone_count.size(), 0);
	return kusudama_limit_cone_count[p_constraint_index];
}

void ManyBoneIK3D::set_kusudama_limit_cone_count(int32_t p_constraint_index, int32_t p_count) {
	ERR_FAIL_INDEX(p_constraint_index, kusudama_limit_cone_count.size());
	ERR_FAIL_INDEX(p_constraint_index, kusudama_limit_cones.size());
	int32_t old_cone_count = kusudama_limit_cones[p_constraint_index].size();
	kusudama_limit_cone_count.write[p_constraint_index] = p_count;
	Vector<Vector4> &cones = kusudama_limit_cones.write[p_constraint_index];
	cones.resize(p_count);
	String bone_name = get_constraint_name(p_constraint_index);
	Transform3D bone_transform = get_bone_direction_transform(p_constraint_index);
	Vector3 forward_axis = -bone_transform.basis.get_column(Vector3::AXIS_Y).normalized();
	for (int32_t cone_i = p_count; cone_i-- > old_cone_count;) {
		Vector4 &cone = cones.write[cone_i];
		cone.x = forward_axis.x;
		cone.y = forward_axis.y;
		cone.z = forward_axis.z;
		cone.w = Math::deg_to_rad(0.0f);
	}
}

real_t ManyBoneIK3D::get_default_damp() const {
	return default_damp;
}

void ManyBoneIK3D::set_default_damp(float p_default_damp) {
	default_damp = p_default_damp;
	set_dirty();
}

StringName ManyBoneIK3D::get_pin_bone_name(int32_t p_effector_index) const {
	ERR_FAIL_INDEX_V(p_effector_index, pins.size(), "");
	Ref<IKEffectorTemplate3D> effector_template = pins[p_effector_index];
	return effector_template->get_name();
}

void ManyBoneIK3D::set_kusudama_limit_cone_radius(int32_t p_effector_index, int32_t p_index, float p_radius) {
	ERR_FAIL_INDEX(p_effector_index, kusudama_limit_cone_count.size());
	ERR_FAIL_INDEX(p_effector_index, kusudama_limit_cones.size());
	ERR_FAIL_INDEX(p_index, kusudama_limit_cone_count[p_effector_index]);
	ERR_FAIL_INDEX(p_index, kusudama_limit_cones[p_effector_index].size());
	Vector4 &cone = kusudama_limit_cones.write[p_effector_index].write[p_index];
	cone.w = p_radius;
	set_dirty();
}

void ManyBoneIK3D::set_kusudama_limit_cone_center(int32_t p_effector_index, int32_t p_index, Vector3 p_center) {
	ERR_FAIL_INDEX(p_effector_index, kusudama_limit_cone_count.size());
	ERR_FAIL_INDEX(p_effector_index, kusudama_limit_cones.size());
	ERR_FAIL_INDEX(p_index, kusudama_limit_cones[p_effector_index].size());
	Vector4 &cone = kusudama_limit_cones.write[p_effector_index].write[p_index];
	if (Math::is_zero_approx(p_center.length_squared())) {
		cone.x = 0;
		cone.y = 1;
		cone.z = 0;
	} else {
		cone.x = p_center.x;
		cone.y = p_center.y;
		cone.z = p_center.z;
	}
	set_dirty();
}

Vector2 ManyBoneIK3D::get_kusudama_twist(int32_t p_index) const {
	ERR_FAIL_INDEX_V(p_index, kusudama_twist.size(), Vector2());
	return kusudama_twist[p_index];
}

void ManyBoneIK3D::set_constraint_name(int32_t p_index, String p_name) {
	ERR_FAIL_INDEX(p_index, constraint_names.size());
	constraint_names.write[p_index] = p_name;
	set_dirty();
}

Vector<Ref<IKBoneSegment3D>> ManyBoneIK3D::get_segmented_skeletons() {
	return segmented_skeletons;
}
float ManyBoneIK3D::get_iterations_per_frame() const {
	return iterations_per_frame;
}

void ManyBoneIK3D::set_iterations_per_frame(const float &p_iterations_per_frame) {
	iterations_per_frame = p_iterations_per_frame;
}

void ManyBoneIK3D::set_pin_bone_name(int32_t p_effector_index, StringName p_name) const {
	ERR_FAIL_INDEX(p_effector_index, pins.size());
	Ref<IKEffectorTemplate3D> effector_template = pins[p_effector_index];
	effector_template->set_name(p_name);
}

void ManyBoneIK3D::set_pin_nodepath(int32_t p_effector_index, NodePath p_node_path) {
	ERR_FAIL_INDEX(p_effector_index, pins.size());
	Node *node = get_node_or_null(p_node_path);
	if (!node) {
		return;
	}
	Ref<IKEffectorTemplate3D> effector_template = pins[p_effector_index];
	effector_template->set_target_node(p_node_path);
}

NodePath ManyBoneIK3D::get_pin_nodepath(int32_t p_effector_index) const {
	ERR_FAIL_INDEX_V(p_effector_index, pins.size(), NodePath());
	Ref<IKEffectorTemplate3D> effector_template = pins[p_effector_index];
	return effector_template->get_target_node();
}

void ManyBoneIK3D::execute(real_t delta) {
	if (!get_skeleton()) {
		return;
	}
	if (get_pin_count() == 0) {
		return;
	}
	if (!segmented_skeletons.size()) {
		set_dirty();
	}
	if (is_dirty) {
		skeleton_changed(get_skeleton());
		is_dirty = false;
		for (int32_t constraint_i = 0; constraint_i < get_constraint_count(); constraint_i++) {
			String constraint_name = get_constraint_name(constraint_i);
			twist_constraint_defaults[constraint_name] = get_constraint_twist_transform(constraint_i);
			orientation_constraint_defaults[constraint_name] = get_constraint_orientation_transform(constraint_i);
			bone_direction_constraint_defaults[constraint_name] = get_bone_direction_transform(constraint_i);
		}
	}
	if (bone_list.size()) {
		Ref<IKNode3D> root_ik_bone = bone_list.write[0]->get_ik_transform();
		if (root_ik_bone.is_null()) {
			return;
		}
		Skeleton3D *skeleton = get_skeleton();
		godot_skeleton_transform->set_transform(skeleton->get_transform());
		godot_skeleton_transform_inverse = skeleton->get_transform().affine_inverse();
	}
	bool has_pins = false;
	for (Ref<IKEffectorTemplate3D> pin : pins) {
		if (pin.is_valid() && !pin->get_name().is_empty()) {
			has_pins = true;
			break;
		}
	}
	if (!has_pins) {
		return;
	}
	update_ik_bones_transform();
	for (Ref<IKBoneSegment3D> segmented_skeleton : segmented_skeletons) {
		if (segmented_skeleton.is_null()) {
			continue;
		}
		segmented_skeleton->update_returnfulness_damp(get_iterations_per_frame());
	}
	for (int32_t i = 0; i < get_iterations_per_frame(); i++) {
		for (Ref<IKBoneSegment3D> segmented_skeleton : segmented_skeletons) {
			if (segmented_skeleton.is_null()) {
				continue;
			}
			segmented_skeleton->segment_solver(bone_damp, get_default_damp(), get_constraint_mode(), i, get_iterations_per_frame());
		}
	}
	update_skeleton_bones_transform();
}

void ManyBoneIK3D::skeleton_changed(Skeleton3D *p_skeleton) {
	if (!p_skeleton) {
		return;
	}
	Vector<int32_t> roots = p_skeleton->get_parentless_bones();
	if (roots.is_empty()) {
		return;
	}
	bone_list.clear();
	segmented_skeletons.clear();
	for (BoneId root_bone_index : roots) {
		StringName parentless_bone = p_skeleton->get_bone_name(root_bone_index);
		Ref<IKBoneSegment3D> segmented_skeleton = Ref<IKBoneSegment3D>(memnew(IKBoneSegment3D(p_skeleton, parentless_bone, pins, this, nullptr, root_bone_index, -1, stabilize_passes)));
		segmented_skeleton->get_root()->get_ik_transform()->set_parent(ik_origin);
		segmented_skeleton->generate_default_segments(pins, root_bone_index, -1, this);
		Vector<Ref<IKBone3D>> new_bone_list;
		segmented_skeleton->create_bone_list(new_bone_list, true, queue_debug_skeleton);
		bone_list.append_array(new_bone_list);
		Vector<Vector<double>> weight_array;
		segmented_skeleton->update_pinned_list(weight_array);
		segmented_skeleton->recursive_create_headings_arrays_for(segmented_skeleton);
		segmented_skeletons.push_back(segmented_skeleton);
	}
	update_ik_bones_transform();
	for (Ref<IKBone3D> &ik_bone_3d : bone_list) {
		ik_bone_3d->update_default_bone_direction_transform(p_skeleton);
	}
	for (int constraint_i = 0; constraint_i < constraint_count; ++constraint_i) {
		String bone = constraint_names[constraint_i];
		BoneId bone_id = p_skeleton->find_bone(bone);
		for (Ref<IKBone3D> &ik_bone_3d : bone_list) {
			if (ik_bone_3d->get_bone_id() != bone_id) {
				continue;
			}
			Ref<IKKusudama3D> constraint = Ref<IKKusudama3D>(memnew(IKKusudama3D()));
			constraint->enable_orientational_limits();

			int32_t cone_count = kusudama_limit_cone_count[constraint_i];
			const Vector<Vector4> &cones = kusudama_limit_cones[constraint_i];
			for (int32_t cone_i = 0; cone_i < cone_count; ++cone_i) {
				const Vector4 &cone = cones[cone_i];
				constraint->add_limit_cone(Vector3(cone.x, cone.y, cone.z), cone.w);
			}

			const Vector2 axial_limit = get_kusudama_twist(constraint_i);
			constraint->enable_axial_limits();
			constraint->set_axial_limits(axial_limit.x, axial_limit.y);
			constraint->set_painfulness(get_kusudama_painfulness(constraint_i));
			ik_bone_3d->add_constraint(constraint);
			constraint->_update_constraint();
			break;
		}
	}
	if (!twist_constraint_defaults.size() && !orientation_constraint_defaults.size() && !bone_direction_constraint_defaults.size()) {
		for (Ref<IKBone3D> &ik_bone_3d : bone_list) {
			ik_bone_3d->update_default_constraint_transform();
		}
		for (int32_t constraint_i = 0; constraint_i < get_constraint_count(); ++constraint_i) {
			String constraint_name = get_constraint_name(constraint_i);
			twist_constraint_defaults[constraint_name] = get_constraint_twist_transform(constraint_i);
			orientation_constraint_defaults[constraint_name] = get_constraint_orientation_transform(constraint_i);
			bone_direction_constraint_defaults[constraint_name] = get_bone_direction_transform(constraint_i);
		}
	}
	for (int32_t constraint_i = 0; constraint_i < get_constraint_count(); ++constraint_i) {
		String constraint_name = get_constraint_name(constraint_i);
		set_constraint_twist_transform(constraint_i, twist_constraint_defaults[constraint_name]);
		set_constraint_orientation_transform(constraint_i, orientation_constraint_defaults[constraint_name]);
		set_bone_direction_transform(constraint_i, bone_direction_constraint_defaults[constraint_name]);
	}
	if (queue_debug_skeleton) {
		queue_debug_skeleton = false;
	}
}

real_t ManyBoneIK3D::get_pin_weight(int32_t p_pin_index) const {
	ERR_FAIL_INDEX_V(p_pin_index, pins.size(), 0.0);
	const Ref<IKEffectorTemplate3D> effector_template = pins[p_pin_index];
	return effector_template->get_weight();
}

void ManyBoneIK3D::set_pin_weight(int32_t p_pin_index, const real_t &p_weight) {
	ERR_FAIL_INDEX(p_pin_index, pins.size());
	Ref<IKEffectorTemplate3D> effector_template = pins[p_pin_index];
	if (effector_template.is_null()) {
		effector_template.instantiate();
		pins.write[p_pin_index] = effector_template;
	}
	effector_template->set_weight(p_weight);
	set_dirty();
}

Vector3 ManyBoneIK3D::get_pin_direction_priorities(int32_t p_pin_index) const {
	ERR_FAIL_INDEX_V(p_pin_index, pins.size(), Vector3(0, 0, 0));
	const Ref<IKEffectorTemplate3D> effector_template = pins[p_pin_index];
	return effector_template->get_direction_priorities();
}

void ManyBoneIK3D::set_pin_direction_priorities(int32_t p_pin_index, const Vector3 &p_priority_direction) {
	ERR_FAIL_INDEX(p_pin_index, pins.size());
	Ref<IKEffectorTemplate3D> effector_template = pins[p_pin_index];
	if (effector_template.is_null()) {
		effector_template.instantiate();
		pins.write[p_pin_index] = effector_template;
	}
	effector_template->set_direction_priorities(p_priority_direction);
	set_dirty();
}

void ManyBoneIK3D::set_dirty() {
	is_dirty = true;
	is_gizmo_dirty = true;
	if (timer->is_inside_tree()) {
		timer->start();
	}
}

void ManyBoneIK3D::_on_timer_timeout() {
	notify_property_list_changed();
}

int32_t ManyBoneIK3D::find_constraint(String p_string) const {
	for (int32_t constraint_i = 0; constraint_i < constraint_count; constraint_i++) {
		if (get_constraint_name(constraint_i) == p_string) {
			return constraint_i;
		}
	}
	return -1;
}

Skeleton3D *ManyBoneIK3D::get_skeleton() const {
	Node *node = get_node_or_null(skeleton_node_path);
	if (!node) {
		return nullptr;
	}
	return cast_to<Skeleton3D>(node);
}

NodePath ManyBoneIK3D::get_skeleton_node_path() {
	return skeleton_node_path;
}

void ManyBoneIK3D::set_skeleton_node_path(NodePath p_skeleton_node_path) {
	skeleton_node_path = p_skeleton_node_path;
	register_skeleton();
	set_dirty(); // Duplicated for ease of verification.
}

void ManyBoneIK3D::_notification(int p_what) {
	switch (p_what) {
		case NOTIFICATION_READY: {
			add_child(timer, true, INTERNAL_MODE_BACK);
			timer->set_owner(get_owner());
			timer->start(true);
			timer->set_wait_time(0.5);
			timer->set_one_shot(true);
			timer->connect("timeout", callable_mp(this, &ManyBoneIK3D::_on_timer_timeout));
			set_process_priority(1);
		} break;
		case NOTIFICATION_ENTER_TREE: {
			set_physics_process_internal(true);
		} break;
		case NOTIFICATION_EXIT_TREE: {
			set_physics_process_internal(false);
			if (timer) {
				timer->queue_free();
			}
		} break;
		case NOTIFICATION_INTERNAL_PHYSICS_PROCESS: {
			execute(get_process_delta_time());
		} break;
	}
}

void ManyBoneIK3D::remove_constraint(int32_t p_index) {
	ERR_FAIL_INDEX(p_index, constraint_count);

	constraint_names.remove_at(p_index);
	kusudama_limit_cone_count.remove_at(p_index);
	kusudama_limit_cones.remove_at(p_index);
	kusudama_twist.remove_at(p_index);
	bone_painfulness.remove_at(p_index);

	constraint_count--;

	set_dirty();
}

void ManyBoneIK3D::_set_bone_count(int32_t p_count) {
	bone_damp.resize(p_count);
	for (int32_t bone_i = p_count; bone_i-- > bone_count;) {
		bone_damp.write[bone_i] = get_default_damp();
	}
	bone_count = p_count;
}

int32_t ManyBoneIK3D::get_bone_count() const {
	return bone_count;
}

Vector<Ref<IKBone3D>> ManyBoneIK3D::get_bone_list() const {
	return bone_list;
}

void ManyBoneIK3D::set_bone_direction_transform(int32_t p_index, Transform3D p_transform) {
	ERR_FAIL_INDEX(p_index, constraint_names.size());
	if (!get_skeleton()) {
		return;
	}
	String bone_name = constraint_names[p_index];
	int32_t bone_index = get_skeleton()->find_bone(bone_name);
	for (Ref<IKBoneSegment3D> segmented_skeleton : segmented_skeletons) {
		if (segmented_skeleton.is_null()) {
			continue;
		}
		Ref<IKBone3D> ik_bone = segmented_skeleton->get_ik_bone(bone_index);
		if (ik_bone.is_null() || ik_bone->get_constraint().is_null()) {
			continue;
		}
		if (ik_bone->get_bone_direction_transform().is_null()) {
			continue;
		}
		ik_bone->get_bone_direction_transform()->set_transform(p_transform);
		break;
	}
}

Transform3D ManyBoneIK3D::get_bone_direction_transform(int32_t p_index) const {
	if (p_index < 0 || p_index >= constraint_names.size() || get_skeleton() == nullptr) {
		return Transform3D();
	}

	String bone_name = constraint_names[p_index];
	int32_t bone_index = get_skeleton()->find_bone(bone_name);
	for (const Ref<IKBoneSegment3D> &segmented_skeleton : segmented_skeletons) {
		if (segmented_skeleton.is_null()) {
			continue;
		}
		Ref<IKBone3D> ik_bone = segmented_skeleton->get_ik_bone(bone_index);
		if (ik_bone.is_null() || ik_bone->get_constraint().is_null()) {
			continue;
		}
		return ik_bone->get_bone_direction_transform()->get_transform();
	}
	return Transform3D();
}

Transform3D ManyBoneIK3D::get_constraint_orientation_transform(int32_t p_index) const {
	ERR_FAIL_INDEX_V(p_index, constraint_names.size(), Transform3D());
	String bone_name = constraint_names[p_index];
	if (!segmented_skeletons.size()) {
		return Transform3D();
	}
	if (!get_skeleton()) {
		return Transform3D();
	}
	for (Ref<IKBoneSegment3D> segmented_skeleton : segmented_skeletons) {
		if (segmented_skeleton.is_null()) {
			continue;
		}
		Ref<IKBone3D> ik_bone = segmented_skeleton->get_ik_bone(get_skeleton()->find_bone(bone_name));
		if (ik_bone.is_null()) {
			continue;
		}
		if (ik_bone->get_constraint().is_null()) {
			continue;
		}
		return ik_bone->get_constraint_orientation_transform()->get_transform();
	}
	return Transform3D();
}

void ManyBoneIK3D::set_constraint_orientation_transform(int32_t p_index, Transform3D p_transform) {
	ERR_FAIL_INDEX(p_index, constraint_names.size());
	String bone_name = constraint_names[p_index];
	if (!get_skeleton()) {
		return;
	}
	for (Ref<IKBoneSegment3D> segmented_skeleton : segmented_skeletons) {
		if (segmented_skeleton.is_null()) {
			continue;
		}
		Ref<IKBone3D> ik_bone = segmented_skeleton->get_ik_bone(get_skeleton()->find_bone(bone_name));
		if (ik_bone.is_null()) {
			continue;
		}
		if (ik_bone->get_constraint().is_null()) {
			continue;
		}
		ik_bone->get_constraint_orientation_transform()->set_transform(p_transform);
		break;
	}
}

Transform3D ManyBoneIK3D::get_constraint_twist_transform(int32_t p_index) const {
	ERR_FAIL_INDEX_V(p_index, constraint_names.size(), Transform3D());
	String bone_name = constraint_names[p_index];
	if (!segmented_skeletons.size()) {
		return Transform3D();
	}
	if (!get_skeleton()) {
		return Transform3D();
	}
	for (Ref<IKBoneSegment3D> segmented_skeleton : segmented_skeletons) {
		if (segmented_skeleton.is_null()) {
			continue;
		}
		Ref<IKBone3D> ik_bone = segmented_skeleton->get_ik_bone(get_skeleton()->find_bone(bone_name));
		if (ik_bone.is_null()) {
			continue;
		}
		if (ik_bone->get_constraint().is_null()) {
			continue;
		}
		return ik_bone->get_constraint_twist_transform()->get_transform();
	}
	return Transform3D();
}

void ManyBoneIK3D::set_constraint_twist_transform(int32_t p_index, Transform3D p_transform) {
	ERR_FAIL_INDEX(p_index, constraint_names.size());
	String bone_name = constraint_names[p_index];
	if (!get_skeleton()) {
		return;
	}
	for (Ref<IKBoneSegment3D> segmented_skeleton : segmented_skeletons) {
		if (segmented_skeleton.is_null()) {
			continue;
		}
		Ref<IKBone3D> ik_bone = segmented_skeleton->get_ik_bone(get_skeleton()->find_bone(bone_name));
		if (ik_bone.is_null()) {
			continue;
		}
		if (ik_bone->get_constraint().is_null()) {
			continue;
		}
		ik_bone->get_constraint_twist_transform()->set_transform(p_transform);
		break;
	}
}

bool ManyBoneIK3D::get_pin_enabled(int32_t p_effector_index) const {
	ERR_FAIL_INDEX_V(p_effector_index, pins.size(), false);
	Ref<IKEffectorTemplate3D> effector_template = pins[p_effector_index];
	return !effector_template->get_target_node().is_empty();
}

void ManyBoneIK3D::register_skeleton() {
	if (!get_pin_count() && !get_constraint_count()) {
		reset_constraints();
	}
	set_dirty();
}

void ManyBoneIK3D::reset_constraints() {
	Skeleton3D *skeleton = get_skeleton();
	if (skeleton) {
		int32_t pin_count = get_pin_count();
		set_pin_count(0);
		set_pin_count(pin_count);
		int32_t get_constraint_count = constraint_names.size();
		set_constraint_count(0);
		set_constraint_count(get_constraint_count);
		_set_bone_count(0);
		_set_bone_count(get_constraint_count);
	}
	set_dirty();
}

bool ManyBoneIK3D::get_constraint_mode() const {
	return is_constraint_mode;
}

void ManyBoneIK3D::set_constraint_mode(bool p_enabled) {
	is_constraint_mode = p_enabled;
}

int32_t ManyBoneIK3D::get_ui_selected_bone() const {
	return ui_selected_bone;
}

void ManyBoneIK3D::set_ui_selected_bone(int32_t p_ui_selected_bone) {
	ui_selected_bone = p_ui_selected_bone;
}

void ManyBoneIK3D::set_kusudama_twist_current(int32_t p_index, real_t p_rotation) {
	ERR_FAIL_INDEX(p_index, constraint_names.size());
	String bone_name = constraint_names[p_index];
	for (Ref<IKBoneSegment3D> segmented_skeleton : segmented_skeletons) {
		if (segmented_skeleton.is_null()) {
			continue;
		}
		Ref<IKBone3D> ik_bone = segmented_skeleton->get_ik_bone(get_skeleton()->find_bone(bone_name));
		if (ik_bone.is_null()) {
			continue;
		}
		if (ik_bone->get_constraint().is_null()) {
			continue;
		}
		ik_bone->get_constraint()->set_current_twist_rotation(ik_bone, p_rotation);
		ik_bone->set_skeleton_bone_pose(get_skeleton());
	}
	set_dirty();
}

real_t ManyBoneIK3D::get_kusudama_twist_current(int32_t p_index) const {
	ERR_FAIL_INDEX_V(p_index, constraint_names.size(), 0.0f);
	String bone_name = constraint_names[p_index];
	if (!segmented_skeletons.size()) {
		return 0;
	}
	for (Ref<IKBoneSegment3D> segmented_skeleton : segmented_skeletons) {
		if (segmented_skeleton.is_null()) {
			continue;
		}
		Ref<IKBone3D> ik_bone = segmented_skeleton->get_ik_bone(get_skeleton()->find_bone(bone_name));
		if (ik_bone.is_null()) {
			continue;
		}
		if (ik_bone->get_constraint().is_null()) {
			continue;
		}
		return CLAMP(ik_bone->get_constraint()->get_current_twist_rotation(ik_bone), 0, 1);
	}
	return 0;
}

void ManyBoneIK3D::set_stabilization_passes(int32_t p_passes) {
	stabilize_passes = p_passes;
	set_dirty();
}

int32_t ManyBoneIK3D::get_stabilization_passes() {
	return stabilize_passes;
}

void ManyBoneIK3D::set_twist_constraint_defaults(Dictionary p_defaults) {
	twist_constraint_defaults = p_defaults;
	set_dirty();
}

Dictionary ManyBoneIK3D::get_twist_constraint_defaults() {
	return twist_constraint_defaults;
}

void ManyBoneIK3D::set_orientation_constraint_defaults(Dictionary p_defaults) {
	orientation_constraint_defaults = p_defaults;
}

Dictionary ManyBoneIK3D::get_orientation_constraint_defaults() {
	return orientation_constraint_defaults;
}

void ManyBoneIK3D::set_bone_direction_constraint_defaults(Dictionary p_defaults) {
	bone_direction_constraint_defaults = p_defaults;
	set_dirty();
}

Dictionary ManyBoneIK3D::get_bone_direction_constraint_defaults() {
	return bone_direction_constraint_defaults;
}

Transform3D ManyBoneIK3D::get_godot_skeleton_transform_inverse() {
	return godot_skeleton_transform_inverse;
}

Ref<IKNode3D> ManyBoneIK3D::get_godot_skeleton_transform() {
	return godot_skeleton_transform;
}

void ManyBoneIK3D::set_setup_humanoid_bones(bool set_targets) {
	is_setup_humanoid_bones = set_targets;
	setup_humanoid_bones(is_setup_humanoid_bones);
	set_dirty();
}

bool ManyBoneIK3D::get_setup_humanoid_bones() const {
	return is_setup_humanoid_bones;
}

void ManyBoneIK3D::setup_humanoid_bones(bool p_set_targets) {
	// **Rotation Twist**
	// | Body Part       | Description                                                                                                                                                                                                                   |
	// |-----------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
	// | Head            | The head can rotate side-to-side up to 60-70 degrees, enabling the character to look left and right.                                                                                                   |
	// | Neck            | The neck can rotate side-to-side up to 50-60 degrees for looking left and right.                                                                                                       |
	// | [Side]UpperLeg  | The upper leg can rotate slightly up to 5-10 degrees for sitting.                                                                                                  |
	// | [Side]Foot      | The foot can also rotate slightly inward or outward (inversion and eversion) up to 10-20 degrees for balance.         |
	// | [Side]UpperArm  | The upper arm can also rotate slightly up to 30-40 degrees for more natural arm movement.                                                                             |
	// | [Side]Hand      | The wrist can also rotate slightly up to 20-30 degrees, enabling the hand to twist inward or outward for grasping and gesturing.                             |

	// **Rotation Swing**
	// | Body Part       | Description                                                                                                                                                                                                                   |
	// |-----------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
	// | Hips            | The hips can tilt forward and backward up to 20-30 degrees, allowing the legs to swing in a wide arc during walking or running. They can also move side-to-side up to 10-20 degrees, enabling the legs to spread apart or come together.                               |
	// | UpperChest      | The upper chest can tilt forward and backward up to 10-20 degrees, allowing for natural breathing and posture adjustments.                                                                                                                         |
	// | Chest           | The chest can tilt forward and backward up to 10-20 degrees, allowing for natural breathing and posture adjustments.                                                                                                                               |
	// | Spine           | The spine can tilt forward and backward up to 35-45 degrees, allowing for bending and straightening of the torso.                                                                                                                                  |
	// | [Side]UpperLeg  | The upper leg can swing forward and backward up to 80-90 degrees, allowing for steps during walking and running.                                                                                                  |
	// | [Side]LowerLeg  | The knee can bend and straighten up to 110-120 degrees, allowing the lower leg to move towards or away from the upper leg during walking, running, and stepping.                                                                                     |
	// | [Side]Foot      | The ankle can tilt up (dorsiflexion) up to 10-20 degrees and down (plantarflexion) up to 35-40 degrees, allowing the foot to step and adjust during walking and running.          |
	// | [Side]Shoulder  | The shoulder can tilt forward and backward up to 160 degrees, allowing the arms to swing in a wide arc. They can also move side-to-side up to 40-50 degrees, enabling the arms to extend outwards or cross over the chest.                                       |
	// | [Side]UpperArm  | The upper arm can swing forward and backward up to 110-120 degrees, allowing for reaching and swinging motions.                                                                             |
	// | [Side]LowerArm  | The elbow can bend and straighten up to 120-130 degrees, allowing the forearm to move towards or away from the upper arm during reaching and swinging motions.                                                                                       |
	// | [Side]Hand      | The wrist can tilt up and down up to 50-60 degrees, allowing the hand to move towards or away from the forearm.
	set_process_thread_group(PROCESS_THREAD_GROUP_SUB_THREAD);
	set_process_thread_group_order(100);
	Skeleton3D *skeleton = cast_to<Skeleton3D>(get_node_or_null(get_skeleton_node_path()));
	ERR_FAIL_NULL(skeleton);
	skeleton->reset_bone_poses();
	Ref<SkeletonProfileHumanoid> humanoid_profile = memnew(SkeletonProfileHumanoid);
	PackedStringArray humanoid_bones;
	if (!p_set_targets) {
		return;
	}
	Vector<String> bones = {
		"Head",
		"LeftHand",
		"RightHand",
		"Hips",
		"LeftFoot",
		"RightFoot",
	};
	set_pin_count(0);
	set_pin_count(bones.size());
	TypedArray<Node> children = find_children("*", "Marker3D");
	for (int i = 0; i < children.size(); ++i) {
		Node *node = cast_to<Node>(children[i]);
		node->queue_free();
	}
	for (int pin_i = 0; pin_i < bones.size(); pin_i++) {
		String bone_name = bones[pin_i];
		Marker3D *marker_3d = memnew(Marker3D);
		marker_3d->set_name(bone_name);
		add_child(marker_3d, true);
		marker_3d->set_owner(get_owner());
		int32_t bone_i = skeleton->find_bone(bone_name);
		Transform3D pose = skeleton->get_global_transform().affine_inverse() * skeleton->get_bone_global_pose_no_override(bone_i);
		marker_3d->set_global_transform(pose);
		set_pin_nodepath(pin_i, get_path_to(marker_3d));
		set_pin_bone_name(pin_i, bone_name);
		set_pin_passthrough_factor(pin_i, 1.0f);
		if (bone_name.ends_with("Foot")) {
			set_pin_passthrough_factor(pin_i, 0.0f);
		} else if (bone_name.ends_with("Hand")) {
			set_pin_passthrough_factor(pin_i, 0.0f);
		} else if (bone_name.find("Head") != -1) {
			set_pin_passthrough_factor(pin_i, 0.0f);
		}
	}
	set_constraint_count(0);
	for (int human_bone_i = 0; human_bone_i < humanoid_profile->get_bone_size(); human_bone_i++) {
		String bone_name = humanoid_profile->get_bone_name(human_bone_i);
		add_constraint();
		set_constraint_name(human_bone_i, bone_name);
	}
	skeleton_changed(get_skeleton());
	for (int constraint_i = 0; constraint_i < get_constraint_count(); constraint_i++) {
		set_kusudama_limit_cone_count(constraint_i, 1);
		const int FIRST_CONE = 0;
		const int SECOND_CONE = 1;
		Transform3D bone_transform = get_bone_direction_transform(constraint_i);
		Vector3 forward = bone_transform.basis.get_column(Vector3::AXIS_Y).normalized();
		Quaternion twist_rotation, swing_rotation;
		IKKusudama3D::get_swing_twist(bone_transform.basis, forward, swing_rotation, twist_rotation);
		Vector3 backwards = -forward;
		set_kusudama_twist(constraint_i, Vector2(Math::deg_to_rad(0.0f), Math::deg_to_rad(180.0f)));
		String bone_name = get_constraint_name(constraint_i);
		if (bone_name == "Spine" || bone_name == "Chest") {
			set_kusudama_twist(constraint_i, Vector2(Math::deg_to_rad(10.0f), Math::deg_to_rad(5.0f)));
			set_kusudama_painfulness(constraint_i, 0.9);
			set_kusudama_limit_cone_center(constraint_i, FIRST_CONE, forward);
			set_kusudama_limit_cone_radius(constraint_i, FIRST_CONE, Math::deg_to_rad(2.5f));
		} else if (bone_name == "UpperChest") {
			set_kusudama_twist(constraint_i, Vector2(Math::deg_to_rad(-20.0f), Math::deg_to_rad(5.0f)));
			set_kusudama_painfulness(constraint_i, 0.9);
			set_kusudama_limit_cone_center(constraint_i, FIRST_CONE, forward);
			set_kusudama_limit_cone_radius(constraint_i, FIRST_CONE, Math::deg_to_rad(10.0f));
		} else if (bone_name == "Hips") {
			set_kusudama_twist(constraint_i, Vector2(Math::deg_to_rad(-90.0f), Math::deg_to_rad(5.0f)));
			set_kusudama_painfulness(constraint_i, 0.8);
			set_kusudama_limit_cone_center(constraint_i, FIRST_CONE, Vector3(0, -1, 0));
			set_kusudama_limit_cone_radius(constraint_i, FIRST_CONE, Math::deg_to_rad(10.0f));
		} else if (bone_name == "Neck") {
			set_kusudama_painfulness(constraint_i, 0.5);
			set_kusudama_limit_cone_center(constraint_i, FIRST_CONE, forward);
			set_kusudama_limit_cone_radius(constraint_i, FIRST_CONE, Math::deg_to_rad(2.5f));
		} else if (bone_name == "Head") {
			set_kusudama_limit_cone_center(constraint_i, FIRST_CONE, forward);
			set_kusudama_limit_cone_radius(constraint_i, FIRST_CONE, Math::deg_to_rad(15.0f));
		} else if (bone_name.ends_with("UpperLeg")) {
			set_kusudama_twist(constraint_i, Vector2(Math::deg_to_rad(0.0f), Math::deg_to_rad(180.0f)));
			set_kusudama_limit_cone_center(constraint_i, FIRST_CONE, backwards);
			set_kusudama_limit_cone_radius(constraint_i, FIRST_CONE, Math::deg_to_rad(25.0f));
		} else if (bone_name.ends_with("LowerLeg")) {
			set_kusudama_twist(constraint_i, Vector2(Math::deg_to_rad(0.0f), Math::deg_to_rad(180.0f)));
			set_kusudama_painfulness(constraint_i, 0.7);
			set_kusudama_limit_cone_count(constraint_i, 2);
			set_kusudama_limit_cone_center(constraint_i, FIRST_CONE, forward);
			set_kusudama_limit_cone_radius(constraint_i, FIRST_CONE, Math::deg_to_rad(2.5f));
			backwards.z += -1;
			backwards.normalize();
			set_kusudama_limit_cone_center(constraint_i, SECOND_CONE, backwards);
			set_kusudama_limit_cone_radius(constraint_i, SECOND_CONE, Math::deg_to_rad(2.5f));
		} else if (bone_name.ends_with("Foot")) {
			set_kusudama_twist(constraint_i, Vector2(Math::deg_to_rad(0.0f), Math::deg_to_rad(180.0f)));
			set_kusudama_painfulness(constraint_i, 0.3);
			backwards.y += -1;
			backwards.z += -1;
			backwards.normalize();
			set_kusudama_limit_cone_center(constraint_i, FIRST_CONE, Vector3(0, -1, 0));
			set_kusudama_limit_cone_radius(constraint_i, FIRST_CONE, Math::deg_to_rad(45.0f));
		} else if (bone_name.ends_with("Shoulder")) {
			set_kusudama_painfulness(constraint_i, 0.6);
			set_kusudama_limit_cone_center(constraint_i, FIRST_CONE, forward);
			set_kusudama_limit_cone_radius(constraint_i, FIRST_CONE, Math::deg_to_rad(30.0f));
		} else if (bone_name.ends_with("UpperArm")) {
			set_kusudama_limit_cone_center(constraint_i, FIRST_CONE, forward);
			set_kusudama_limit_cone_radius(constraint_i, FIRST_CONE, Math::deg_to_rad(90.0f));
		} else if (bone_name.ends_with("LowerArm")) {
			set_kusudama_limit_cone_count(constraint_i, 2);
			set_kusudama_limit_cone_center(constraint_i, FIRST_CONE, forward);
			set_kusudama_limit_cone_radius(constraint_i, FIRST_CONE, Math::deg_to_rad(2.5f));
			if (bone_name.begins_with("Left")) {
				backwards.x += -1;
			} else {
				backwards.x += 1;
			}
			backwards.z += -1;
			backwards.normalize();
			set_kusudama_limit_cone_center(constraint_i, SECOND_CONE, backwards);
			set_kusudama_limit_cone_radius(constraint_i, SECOND_CONE, Math::deg_to_rad(2.5f));
		} else if (bone_name.ends_with("Hand")) {
			set_kusudama_painfulness(constraint_i, 0.4);
			set_kusudama_limit_cone_center(constraint_i, FIRST_CONE, forward);
			set_kusudama_limit_cone_radius(constraint_i, FIRST_CONE, Math::deg_to_rad(60.0f));
		} else if (bone_name.ends_with("Thumb")) {
			set_kusudama_limit_cone_center(constraint_i, FIRST_CONE, forward);
			set_kusudama_limit_cone_radius(constraint_i, FIRST_CONE, Math::deg_to_rad(90.0f));
		} else if (bone_name.ends_with("Eye")) {
			set_kusudama_limit_cone_center(constraint_i, FIRST_CONE, forward);
			set_kusudama_limit_cone_radius(constraint_i, FIRST_CONE, Math::deg_to_rad(10.0f));
		} else {
			set_kusudama_limit_cone_center(constraint_i, FIRST_CONE, forward);
			set_kusudama_limit_cone_radius(constraint_i, FIRST_CONE, Math::deg_to_rad(2.5f));
		}
	}
	is_setup_humanoid_bones = false;
}

void ManyBoneIK3D::set_kusudama_painfulness(int32_t p_index, real_t p_painfulness) {
	ERR_FAIL_INDEX(p_index, constraint_names.size());
	String bone_name = constraint_names[p_index];
	bone_painfulness.write[p_index] = p_painfulness;
	for (Ref<IKBoneSegment3D> segmented_skeleton : segmented_skeletons) {
		if (segmented_skeleton.is_null()) {
			continue;
		}
		Ref<IKBone3D> ik_bone = segmented_skeleton->get_ik_bone(get_skeleton()->find_bone(bone_name));
		if (ik_bone.is_null()) {
			continue;
		}
		if (ik_bone->get_constraint().is_null()) {
			continue;
		}
		ik_bone->get_constraint()->set_painfulness(p_painfulness);
		ik_bone->set_skeleton_bone_pose(get_skeleton());
		break;
	}
	set_dirty();
}

real_t ManyBoneIK3D::get_kusudama_painfulness(int32_t p_index) const {
	ERR_FAIL_INDEX_V(p_index, constraint_names.size(), 0.0f);
	return bone_painfulness[p_index];
}

void ManyBoneIK3D::add_constraint() {
	int32_t old_count = constraint_count;
	set_constraint_count(constraint_count + 1);
	constraint_names.write[old_count] = String();
	kusudama_limit_cone_count.write[old_count] = 0;
	kusudama_limit_cones.write[old_count].resize(1);
	kusudama_limit_cones.write[old_count].write[0] = Vector4(0, 1, 0, Math_PI);
	kusudama_twist.write[old_count] = Vector2(0, Math_TAU - CMP_EPSILON);
	bone_painfulness.write[old_count] = 0.0f;
	set_dirty();
}
