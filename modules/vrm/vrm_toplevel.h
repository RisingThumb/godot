/**************************************************************************/
/*  vrm_toplevel.h                                                        */
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

#ifndef VRM_VRMTOPLEVEL_H
#define VRM_VRMTOPLEVEL_H

#include "core/math/quaternion.h"
#include "core/math/transform_3d.h"
#include "scene/3d/node_3d.h"

class VRMTopLevel : public Node3D {
	GDCLASS(VRMTopLevel, Node3D);

private:
	NodePath vrm_skeleton;
	NodePath vrm_animplayer;

	Ref<Resource> vrm_meta;

protected:
	static void _bind_methods();

public:
	VRMTopLevel() {}
	~VRMTopLevel() {}
	// Getters and setters for the exported variables
	void set_vrm_skeleton(const NodePath &path);
	NodePath get_vrm_skeleton() const;

	void set_vrm_animplayer(const NodePath &path);
	NodePath get_vrm_animplayer() const;

	void set_vrm_meta(const Ref<Resource> &meta);
	Ref<Resource> get_vrm_meta() const;
};

#endif // VRM_VRMTOPLEVEL_H