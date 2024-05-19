/**************************************************************************/
/*  domain.cpp                                                            */
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

#include "domain.h"

void Domain::_bind_methods() {
	ClassDB::bind_static_method("Domain", D_METHOD("method_verify_goal", "state", "method", "state_var", "arguments", "desired_values", "depth", "verbose"), &Domain::method_verify_goal);
	ClassDB::bind_static_method("Domain", D_METHOD("method_goals_not_achieved", "state", "multigoal"), &Domain::method_goals_not_achieved);
	ClassDB::bind_static_method("Domain", D_METHOD("method_verify_multigoal", "state", "method", "multigoal", "depth", "verbose"), &Domain::method_verify_multigoal);

	ClassDB::bind_method(D_METHOD("print_domain"), &Domain::print_domain);
	ClassDB::bind_method(D_METHOD("print_actions"), &Domain::print_actions);
	ClassDB::bind_method(D_METHOD("print_task_methods"), &Domain::print_task_methods);
	ClassDB::bind_method(D_METHOD("print_unigoal_methods"), &Domain::print_unigoal_methods);
	ClassDB::bind_method(D_METHOD("print_multigoal_methods"), &Domain::print_multigoal_methods);
	ClassDB::bind_method(D_METHOD("print_methods"), &Domain::print_methods);

	ClassDB::bind_method(D_METHOD("set_verbose", "value"), &Domain::set_verbose);
	ClassDB::bind_method(D_METHOD("get_verbose"), &Domain::get_verbose);

	ClassDB::bind_method(D_METHOD("set_actions", "value"), &Domain::set_actions);
	ClassDB::bind_method(D_METHOD("get_actions"), &Domain::get_actions);

	ClassDB::bind_method(D_METHOD("set_task_methods", "value"), &Domain::set_task_methods);
	ClassDB::bind_method(D_METHOD("get_task_methods"), &Domain::get_task_methods);

	ClassDB::bind_method(D_METHOD("set_unigoal_methods", "value"), &Domain::set_unigoal_methods);
	ClassDB::bind_method(D_METHOD("get_unigoal_methods"), &Domain::get_unigoal_methods);

	ClassDB::bind_method(D_METHOD("set_multigoal_methods", "value"), &Domain::set_multigoal_methods);
	ClassDB::bind_method(D_METHOD("get_multigoal_methods"), &Domain::get_multigoal_methods);

	ADD_PROPERTY(PropertyInfo(Variant::INT, "verbose"), "set_verbose", "get_verbose");
	ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "actions"), "set_actions", "get_actions");
	ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "tasks"), "set_task_methods", "get_task_methods");
	ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "unigoal_methods"), "set_unigoal_methods", "get_unigoal_methods");
	ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "multigoal_methods"), "set_multigoal_methods", "get_multigoal_methods");
}

Variant Domain::method_verify_goal(Dictionary p_state, String p_method, String p_state_variable, String p_arguments, Variant p_desired_value, int p_depth, int p_verbose) {
	Dictionary state_dict = p_state[p_state_variable];

	if (state_dict[p_arguments] != p_desired_value) {
		if (p_verbose >= 3) {
			print_line(vformat("Depth %d: method %s didn't achieve\nGoal %s[%s] = %s", p_depth, p_method, p_state_variable, p_arguments, p_desired_value));
		}
		return false;
	}

	if (p_verbose >= 3) {
		print_line(vformat("Depth %d: method %s achieved\nGoal %s[%s] = %s", p_depth, p_method, p_state_variable, p_arguments, p_desired_value));
	}
	return Array();
}

Dictionary Domain::method_goals_not_achieved(Dictionary p_state, Ref<Multigoal> p_multigoal) {
	Dictionary incomplete;
	Array keys = p_multigoal->get_state().keys();
	for (int i = 0; i < keys.size(); ++i) {
		String n = keys[i];
		Dictionary sub_dict = p_multigoal->get_state()[n];
		Array sub_keys = sub_dict.keys();
		for (int j = 0; j < sub_keys.size(); ++j) {
			String arg = sub_keys[j];
			Variant val = sub_dict[arg];
			if (p_state[n].get_type() == Variant::DICTIONARY && Dictionary(p_state[n]).has(arg) && val != Dictionary(p_state[n])[arg]) {
				if (!incomplete.has(n)) {
					incomplete[n] = Dictionary();
				}
				Dictionary temp = incomplete[n];
				temp[arg] = val;
				incomplete[n] = temp;
			}
		}
	}
	return incomplete;
}

Variant Domain::method_verify_multigoal(Dictionary p_state, String p_method, Ref<Multigoal> p_multigoal, int p_depth, int p_verbose) {
	Dictionary goal_dict = method_goals_not_achieved(p_state, p_multigoal);
	if (!goal_dict.is_empty()) {
		if (p_verbose >= 3) {
			print_line(vformat("Depth %d: method %s didn't achieve %s", p_depth, p_method, p_multigoal));
		}
		return false;
	}

	if (p_verbose >= 3) {
		print_line(vformat("Depth %d: method %s achieved %s", p_depth, p_method, p_multigoal));
	}
	return Array();
}

Domain::Domain(String p_name) {
	set_name(p_name);
	task_method_dictionary["_verify_g"] = varray(callable_mp_static(&Domain::method_verify_goal));
	task_method_dictionary["_verify_mg"] = varray(callable_mp_static(&Domain::method_verify_multigoal));
}

void Domain::print_domain() const {
	print_line(vformat("Domain name: %s", get_name()));
	print_actions();
	print_methods();
}

void Domain::print_actions() const {
	Dictionary action_dict = get_actions();
	if (action_dict.is_empty()) {
		print_line("-- There are no actions --");
		return;
	}
	String actions = "-- Actions: ";

	Array keys = action_dict.keys();
	for (int i = 0; i < keys.size(); ++i) {
		if (i != 0) {
			actions += ", ";
		}
		String key = keys[i];
		actions += key;
	}
	print_line(actions);
}

void Domain::print_task_methods() const {
	if (get_task_methods().is_empty()) {
		print_line("-- There are no task methods --");
		return;
	}
	print_line("Task name:         Relevant task methods:");
	print_line("---------------    ----------------------");

	String string_array;
	Array keys = get_task_methods().keys();
	for (int i = 0; i < keys.size(); ++i) {
		string_array += String(keys[i]) + ", ";
	}
	print_line(string_array.substr(0, string_array.length() - 2)); // Remove last comma and space

	print_line("");
}

void Domain::print_unigoal_methods() const {
	if (get_unigoal_methods().is_empty()) {
		print_line("-- There are no unigoal methods --");
		return;
	}
	print_line("Blackboard var name:    Relevant unigoal methods:");
	print_line("---------------    -------------------------");
	Array keys = get_unigoal_methods().keys();
	for (int j = 0; j < keys.size(); ++j) {
		String string_array;
		Array methods = get_unigoal_methods()[keys[j]];
		for (int i = 0; i < methods.size(); ++i) {
			String m = methods[i];
			string_array += m + ", ";
		}
		string_array = string_array.substr(0, string_array.length() - 2);
		print_line(String(keys[j]) + ":    " + string_array);
	}

	print_line("");
}

void Domain::print_multigoal_methods() const {
	if (get_multigoal_methods().is_empty()) {
		print_line("-- There are no multigoal methods --");
		return;
	}
	String string_array;
	Array methods = get_multigoal_methods();
	for (int i = 0; i < methods.size(); ++i) {
		Variant m = methods[i];
		string_array += String(m) + ", ";
	}
	print_line("-- Multigoal methods: " + string_array.substr(0, string_array.length() - 2)); // Remove last comma and space
}

void Domain::print_methods() const {
	print_task_methods();
	print_unigoal_methods();
	print_multigoal_methods();
}
