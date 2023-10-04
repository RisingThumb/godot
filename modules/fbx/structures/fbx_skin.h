/**************************************************************************/
/*  fbx_skin.h                                                            */
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

#ifndef FBX_SKIN_H
#define FBX_SKIN_H

#include "../fbx_defines.h"

#include "core/io/resource.h"

template <typename T>
class TypedArray;

class FBXSkin : public Resource {
	GDCLASS(FBXSkin, Resource);
	friend class FBXDocument;

private:
	// The "skeleton" property defined in the gltf spec. -1 = Scene Root
	FBXNodeIndex skin_root = -1;

	Vector<FBXNodeIndex> joints_original;
	Vector<Transform3D> inverse_binds;

	// Note: joints + non_joints should form a complete subtree, or subtrees
	// with a common parent

	// All nodes that are skins that are caught in-between the original joints
	// (inclusive of joints_original)
	Vector<FBXNodeIndex> joints;

	// All Nodes that are caught in-between skin joint nodes, and are not
	// defined as joints by any skin
	Vector<FBXNodeIndex> non_joints;

	// The roots of the skin. In the case of multiple roots, their parent *must*
	// be the same (the roots must be siblings)
	Vector<FBXNodeIndex> roots;

	// The GLTF Skeleton this Skin points to (after we determine skeletons)
	FBXSkeletonIndex skeleton = -1;

	// A mapping from the joint indices (in the order of joints_original) to the
	// Godot Skeleton's bone_indices
	HashMap<int, int> joint_i_to_bone_i;
	HashMap<int, StringName> joint_i_to_name;

	// The Actual Skin that will be created as a mapping between the IBM's of
	// this skin to the generated skeleton for the mesh instances.
	Ref<Skin> godot_skin;

protected:
	static void _bind_methods();

public:
	FBXNodeIndex get_skin_root();
	void set_skin_root(FBXNodeIndex p_skin_root);

	Vector<FBXNodeIndex> get_joints_original();
	void set_joints_original(Vector<FBXNodeIndex> p_joints_original);

	TypedArray<Transform3D> get_inverse_binds();
	void set_inverse_binds(TypedArray<Transform3D> p_inverse_binds);

	Vector<FBXNodeIndex> get_joints();
	void set_joints(Vector<FBXNodeIndex> p_joints);

	Vector<FBXNodeIndex> get_non_joints();
	void set_non_joints(Vector<FBXNodeIndex> p_non_joints);

	Vector<FBXNodeIndex> get_roots();
	void set_roots(Vector<FBXNodeIndex> p_roots);

	int get_skeleton();
	void set_skeleton(int p_skeleton);

	Dictionary get_joint_i_to_bone_i();
	void set_joint_i_to_bone_i(Dictionary p_joint_i_to_bone_i);

	Dictionary get_joint_i_to_name();
	void set_joint_i_to_name(Dictionary p_joint_i_to_name);

	Ref<Skin> get_godot_skin();
	void set_godot_skin(Ref<Skin> p_godot_skin);
};

#endif // FBX_SKIN_H
