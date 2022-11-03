/*************************************************************************/
/*  ik_effector_template.h                                               */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                      https://godotengine.org                          */
/*************************************************************************/
/* Copyright (c) 2007-2019 Juan Linietsky, Ariel Manzur.                 */
/* Copyright (c) 2014-2019 Godot Engine contributors (cf. AUTHORS.md)    */
/*                                                                       */
/* Permission is hereby granted, free of charge, to any person obtaining */
/* a copy of this software and associated documentation files (the       */
/* "Software"), to deal in the Software without restriction, including   */
/* without limitation the rights to use, copy, modify, merge, publish,   */
/* distribute, sublicense, and/or sell copies of the Software, and to    */
/* permit persons to whom the Software is furnished to do so, subject to */
/* the following conditions:                                             */
/*                                                                       */
/* The above copyright notice and this permission notice shall be        */
/* included in all copies or substantial portions of the Software.       */
/*                                                                       */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,       */
/* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF    */
/* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.*/
/* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY  */
/* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,  */
/* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE     */
/* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                */
/*************************************************************************/

#pragma once

#include "core/io/resource.h"
#include "core/object/ref_counted.h"
#include "core/string/node_path.h"

class IKEffectorTemplate : public Resource {
	GDCLASS(IKEffectorTemplate, Resource);

	NodePath target_node;
	real_t depth_falloff = 0.0f;
	real_t weight = 1.0f;
	Vector3 priority_direction = Vector3(0.5f, 0.0f, 0.5f);

protected:
	static void _bind_methods();

public:
	NodePath get_target_node() const;
	void set_target_node(NodePath p_node_path);
	float get_depth_falloff() const;
	void set_depth_falloff(float p_depth_falloff);
	real_t get_weight() const { return weight; }
	void set_weight(real_t p_weight) { weight = p_weight; }
	Vector3 get_direction_priorities() const { return priority_direction; }
	void set_direction_priorities(Vector3 p_priority_direction) { priority_direction = p_priority_direction; }

	IKEffectorTemplate();
};
