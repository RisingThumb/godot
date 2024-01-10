#include "fbx_skin_utility.h"

#include "modules/fbx/fbx_defines.h"
#include "scene/3d/skeleton_3d.h"
#include "structures/fbx_node.h"
#include "structures/fbx_skin.h"
#include "fbx_document.h"
#include "fbx_state.h"

FBXNodeIndex FBXSkinUtility::_find_highest_node(Vector<Ref<FBXNode>> &r_nodes, const Vector<FBXNodeIndex> &p_subset) {
	int highest = -1;
	FBXNodeIndex best_node = -1;

	for (int i = 0; i < p_subset.size(); ++i) {
		const FBXNodeIndex node_i = p_subset[i];
		const Ref<FBXNode> node = r_nodes[node_i];

		if (highest == -1 || node->height < highest) {
			highest = node->height;
			best_node = node_i;
		}
	}

	return best_node;
}

bool FBXSkinUtility::_capture_nodes_in_skin(const Vector<Ref<FBXNode>> &nodes, Ref<FBXSkin> p_skin, const FBXNodeIndex p_node_index) {
	bool found_joint = false;
	Ref<FBXNode> current_node = nodes[p_node_index];

	for (int i = 0; i < current_node->children.size(); ++i) {
		found_joint |= _capture_nodes_in_skin(nodes, p_skin, current_node->children[i]);
	}

	if (found_joint) {
		// Mark it if we happen to find another skins joint...
		if (current_node->joint && p_skin->joints.find(p_node_index) < 0) {
			p_skin->joints.push_back(p_node_index);
		} else if (p_skin->non_joints.find(p_node_index) < 0) {
			p_skin->non_joints.push_back(p_node_index);
		}
	}

	if (p_skin->joints.find(p_node_index) > 0) {
		return true;
	}

	return false;
}
void FBXSkinUtility::_capture_nodes_for_multirooted_skin(Vector<Ref<FBXNode>> &r_nodes, Ref<FBXSkin> p_skin) {
	DisjointSet<FBXNodeIndex> disjoint_set;

	for (int i = 0; i < p_skin->joints.size(); ++i) {
		const FBXNodeIndex node_index = p_skin->joints[i];
		const FBXNodeIndex parent = r_nodes[node_index]->parent;
		disjoint_set.insert(node_index);

		if (p_skin->joints.find(parent) >= 0) {
			disjoint_set.create_union(parent, node_index);
		}
	}

	Vector<FBXNodeIndex> roots;
	disjoint_set.get_representatives(roots);

	if (roots.size() <= 1) {
		return;
	}

	int maxHeight = -1;

	// Determine the max height rooted tree
	for (int i = 0; i < roots.size(); ++i) {
		const FBXNodeIndex root = roots[i];

		if (maxHeight == -1 || r_nodes[root]->height < maxHeight) {
			maxHeight = r_nodes[root]->height;
		}
	}

	// Go up the tree till all of the multiple roots of the skin are at the same hierarchy level.
	// This sucks, but 99% of all game engines (not just Godot) would have this same issue.
	for (int i = 0; i < roots.size(); ++i) {
		FBXNodeIndex current_node = roots[i];
		while (r_nodes[current_node]->height > maxHeight) {
			FBXNodeIndex parent = r_nodes[current_node]->parent;

			if (r_nodes[parent]->joint && p_skin->joints.find(parent) < 0) {
				p_skin->joints.push_back(parent);
			} else if (p_skin->non_joints.find(parent) < 0) {
				p_skin->non_joints.push_back(parent);
			}

			current_node = parent;
		}

		// replace the roots
		roots.write[i] = current_node;
	}

	// Climb up the tree until they all have the same parent
	bool all_same;

	do {
		all_same = true;
		const FBXNodeIndex first_parent = r_nodes[roots[0]]->parent;

		for (int i = 1; i < roots.size(); ++i) {
			all_same &= (first_parent == r_nodes[roots[i]]->parent);
		}

		if (!all_same) {
			for (int i = 0; i < roots.size(); ++i) {
				const FBXNodeIndex current_node = roots[i];
				const FBXNodeIndex parent = r_nodes[current_node]->parent;

				if (r_nodes[parent]->joint && p_skin->joints.find(parent) < 0) {
					p_skin->joints.push_back(parent);
				} else if (p_skin->non_joints.find(parent) < 0) {
					p_skin->non_joints.push_back(parent);
				}

				roots.write[i] = parent;
			}
		}

	} while (!all_same);
}
Error FBXSkinUtility::_expand_skin(Vector<Ref<FBXNode>> &r_nodes, Ref<FBXSkin> p_skin) {
	_capture_nodes_for_multirooted_skin(r_nodes, p_skin);

	// Grab all nodes that lay in between skin joints/nodes
	DisjointSet<FBXNodeIndex> disjoint_set;

	Vector<FBXNodeIndex> all_skin_nodes;
	all_skin_nodes.append_array(p_skin->joints);
	all_skin_nodes.append_array(p_skin->non_joints);

	for (int i = 0; i < all_skin_nodes.size(); ++i) {
		const FBXNodeIndex node_index = all_skin_nodes[i];
		const FBXNodeIndex parent = r_nodes[node_index]->parent;
		disjoint_set.insert(node_index);

		if (all_skin_nodes.find(parent) >= 0) {
			disjoint_set.create_union(parent, node_index);
		}
	}

	Vector<FBXNodeIndex> out_owners;
	disjoint_set.get_representatives(out_owners);

	Vector<FBXNodeIndex> out_roots;

	for (int i = 0; i < out_owners.size(); ++i) {
		Vector<FBXNodeIndex> set;
		disjoint_set.get_members(set, out_owners[i]);

		const FBXNodeIndex root = _find_highest_node(r_nodes, set);
		ERR_FAIL_COND_V(root < 0, FAILED);
		out_roots.push_back(root);
	}

	out_roots.sort();

	for (int i = 0; i < out_roots.size(); ++i) {
		_capture_nodes_in_skin(r_nodes, p_skin, out_roots[i]);
	}

	p_skin->roots = out_roots;

	return OK;
}
Error FBXSkinUtility::_verify_skin(Vector<Ref<FBXNode>> &r_nodes, Ref<FBXSkin> p_skin) {
	// This may seem duplicated from expand_skins, but this is really a sanity check! (so it kinda is)
	// In case additional interpolating logic is added to the skins, this will help ensure that you
	// do not cause it to self implode into a fiery blaze

	// We are going to re-calculate the root nodes and compare them to the ones saved in the skin,
	// then ensure the multiple trees (if they exist) are on the same sublevel

	// Grab all nodes that lay in between skin joints/nodes
	DisjointSet<FBXNodeIndex> disjoint_set;

	Vector<FBXNodeIndex> all_skin_nodes;
	all_skin_nodes.append_array(p_skin->joints);
	all_skin_nodes.append_array(p_skin->non_joints);

	for (int i = 0; i < all_skin_nodes.size(); ++i) {
		const FBXNodeIndex node_index = all_skin_nodes[i];
		const FBXNodeIndex parent = r_nodes[node_index]->parent;
		disjoint_set.insert(node_index);

		if (all_skin_nodes.find(parent) >= 0) {
			disjoint_set.create_union(parent, node_index);
		}
	}

	Vector<FBXNodeIndex> out_owners;
	disjoint_set.get_representatives(out_owners);

	Vector<FBXNodeIndex> out_roots;

	for (int i = 0; i < out_owners.size(); ++i) {
		Vector<FBXNodeIndex> set;
		disjoint_set.get_members(set, out_owners[i]);

		const FBXNodeIndex root = _find_highest_node(r_nodes, set);
		ERR_FAIL_COND_V(root < 0, FAILED);
		out_roots.push_back(root);
	}

	out_roots.sort();

	ERR_FAIL_COND_V(out_roots.size() == 0, FAILED);

	// Make sure the roots are the exact same (they better be)
	ERR_FAIL_COND_V(out_roots.size() != p_skin->roots.size(), FAILED);
	for (int i = 0; i < out_roots.size(); ++i) {
		ERR_FAIL_COND_V(out_roots[i] != p_skin->roots[i], FAILED);
	}

	// Single rooted skin? Perfectly ok!
	if (out_roots.size() == 1) {
		return OK;
	}

	// Make sure all parents of a multi-rooted skin are the SAME
	const FBXNodeIndex parent = r_nodes[out_roots[0]]->parent;
	for (int i = 1; i < out_roots.size(); ++i) {
		if (r_nodes[out_roots[i]]->parent != parent) {
			return FAILED;
		}
	}

	return OK;
}
void FBXSkinUtility::_recurse_children(
		Vector<Ref<FBXNode>> &nodes,
		const FBXNodeIndex p_node_index,
		RBSet<FBXNodeIndex> &p_all_skin_nodes,
		HashSet<FBXNodeIndex> &p_child_visited_set) {
	if (p_child_visited_set.has(p_node_index)) {
		return;
	}
	p_child_visited_set.insert(p_node_index);

	Ref<FBXNode> current_node = nodes[p_node_index];
	for (int i = 0; i < current_node->children.size(); ++i) {
		_recurse_children(nodes, current_node->children[i], p_all_skin_nodes, p_child_visited_set);
	}

	// Continue to use 'current_node' for clarity and direct access.
	if (current_node->skin < 0 || current_node->mesh < 0 || !current_node->children.is_empty()) {
		p_all_skin_nodes.insert(p_node_index);
	}
}
Error FBXSkinUtility::_determine_skeletons(
		Vector<Ref<FBXSkin>> &skins,
		Vector<Ref<FBXNode>> &nodes,
		Vector<Ref<FBXSkeleton>> &skeletons) {
	// Using a disjoint set, we are going to potentially combine all skins that are actually branches
	// of a main skeleton, or treat skins defining the same set of nodes as ONE skeleton.
	// This is another unclear issue caused by the current glTF specification.

	DisjointSet<FBXNodeIndex> skeleton_sets;

	for (FBXSkinIndex skin_i = 0; skin_i < skins.size(); ++skin_i) {
		const Ref<FBXSkin> skin = skins[skin_i];

		HashSet<FBXNodeIndex> child_visited_set;
		RBSet<FBXNodeIndex> all_skin_nodes;
		for (int i = 0; i < skin->joints.size(); ++i) {
			all_skin_nodes.insert(skin->joints[i]);
			FBXSkinUtility::_recurse_children(nodes, skin->joints[i], all_skin_nodes, child_visited_set);
		}
		for (int i = 0; i < skin->non_joints.size(); ++i) {
			all_skin_nodes.insert(skin->non_joints[i]);
			FBXSkinUtility::_recurse_children(nodes, skin->non_joints[i], all_skin_nodes, child_visited_set);
		}
		for (FBXNodeIndex node_index : all_skin_nodes) {
			const FBXNodeIndex parent = nodes[node_index]->parent;
			skeleton_sets.insert(node_index);

			if (all_skin_nodes.has(parent)) {
				skeleton_sets.create_union(parent, node_index);
			}
		}

		// We are going to connect the separate skin subtrees in each skin together
		// so that the final roots are entire sets of valid skin trees
		for (int i = 1; i < skin->roots.size(); ++i) {
			skeleton_sets.create_union(skin->roots[0], skin->roots[i]);
		}
	}

	{ // attempt to joint all touching subsets (siblings/parent are part of another skin)
		Vector<FBXNodeIndex> groups_representatives;
		skeleton_sets.get_representatives(groups_representatives);

		Vector<FBXNodeIndex> highest_group_members;
		Vector<Vector<FBXNodeIndex>> groups;
		for (int i = 0; i < groups_representatives.size(); ++i) {
			Vector<FBXNodeIndex> group;
			skeleton_sets.get_members(group, groups_representatives[i]);
			highest_group_members.push_back(FBXSkinUtility::_find_highest_node(nodes, group));
			groups.push_back(group);
		}

		for (int i = 0; i < highest_group_members.size(); ++i) {
			const FBXNodeIndex node_i = highest_group_members[i];

			// Attach any siblings together (this needs to be done n^2/2 times)
			for (int j = i + 1; j < highest_group_members.size(); ++j) {
				const FBXNodeIndex node_j = highest_group_members[j];

				// Even if they are siblings under the root! :)
				if (nodes[node_i]->parent == nodes[node_j]->parent) {
					skeleton_sets.create_union(node_i, node_j);
				}
			}

			// Attach any parenting going on together (we need to do this n^2 times)
			const FBXNodeIndex node_i_parent = nodes[node_i]->parent;
			if (node_i_parent >= 0) {
				for (int j = 0; j < groups.size() && i != j; ++j) {
					const Vector<FBXNodeIndex> &group = groups[j];

					if (group.find(node_i_parent) >= 0) {
						const FBXNodeIndex node_j = highest_group_members[j];
						skeleton_sets.create_union(node_i, node_j);
					}
				}
			}
		}
	}

	// At this point, the skeleton groups should be finalized
	Vector<FBXNodeIndex> skeleton_owners;
	skeleton_sets.get_representatives(skeleton_owners);

	// Mark all the skins actual skeletons, after we have merged them
	for (FBXSkeletonIndex skel_i = 0; skel_i < skeleton_owners.size(); ++skel_i) {
		const FBXNodeIndex skeleton_owner = skeleton_owners[skel_i];
		Ref<FBXSkeleton> skeleton;
		skeleton.instantiate();

		Vector<FBXNodeIndex> skeleton_nodes;
		skeleton_sets.get_members(skeleton_nodes, skeleton_owner);

		for (FBXSkinIndex skin_i = 0; skin_i < skins.size(); ++skin_i) {
			Ref<FBXSkin> skin = skins.write[skin_i];

			// If any of the the skeletons nodes exist in a skin, that skin now maps to the skeleton
			for (int i = 0; i < skeleton_nodes.size(); ++i) {
				FBXNodeIndex skel_node_i = skeleton_nodes[i];
				if (skin->joints.find(skel_node_i) >= 0 || skin->non_joints.find(skel_node_i) >= 0) {
					skin->skeleton = skel_i;
					continue;
				}
			}
		}

		Vector<FBXNodeIndex> non_joints;
		for (int i = 0; i < skeleton_nodes.size(); ++i) {
			const FBXNodeIndex node_i = skeleton_nodes[i];

			if (nodes[node_i]->joint) {
				skeleton->joints.push_back(node_i);
			} else {
				non_joints.push_back(node_i);
			}
		}

		skeletons.push_back(skeleton);

		FBXSkinUtility::_reparent_non_joint_skeleton_subtrees(nodes, skeletons.write[skel_i], non_joints);
	}

	for (FBXSkeletonIndex skel_i = 0; skel_i < skeletons.size(); ++skel_i) {
		Ref<FBXSkeleton> skeleton = skeletons.write[skel_i];

		for (int i = 0; i < skeleton->joints.size(); ++i) {
			const FBXNodeIndex node_i = skeleton->joints[i];
			Ref<FBXNode> node = nodes[node_i];

			ERR_FAIL_COND_V(!node->joint, ERR_PARSE_ERROR);
			ERR_FAIL_COND_V(node->skeleton >= 0, ERR_PARSE_ERROR);
			node->skeleton = skel_i;
		}

		ERR_FAIL_COND_V(FBXSkinUtility::_determine_skeleton_roots(nodes, skeletons, skel_i), ERR_PARSE_ERROR);
	}

	return OK;
}
Error FBXSkinUtility::_reparent_non_joint_skeleton_subtrees(
		Vector<Ref<FBXNode>> &nodes,
		Ref<FBXSkeleton> p_skeleton,
		const Vector<FBXNodeIndex> &p_non_joints) {
	DisjointSet<FBXNodeIndex> subtree_set;

	// Populate the disjoint set with ONLY non joints that are in the skeleton hierarchy (non_joints vector)
	// This way we can find any joints that lie in between joints, as the current glTF specification
	// mentions nothing about non-joints being in between joints of the same skin. Hopefully one day we
	// can remove this code.

	// skinD depicted here explains this issue:
	// https://github.com/KhronosGroup/glTF-Asset-Generator/blob/master/Output/Positive/Animation_Skin

	for (int i = 0; i < p_non_joints.size(); ++i) {
		const FBXNodeIndex node_i = p_non_joints[i];

		subtree_set.insert(node_i);

		const FBXNodeIndex parent_i = nodes[node_i]->parent;
		if (parent_i >= 0 && p_non_joints.find(parent_i) >= 0 && !nodes[parent_i]->joint) {
			subtree_set.create_union(parent_i, node_i);
		}
	}

	// Find all the non joint subtrees and re-parent them to a new "fake" joint

	Vector<FBXNodeIndex> non_joint_subtree_roots;
	subtree_set.get_representatives(non_joint_subtree_roots);

	for (int root_i = 0; root_i < non_joint_subtree_roots.size(); ++root_i) {
		const FBXNodeIndex subtree_root = non_joint_subtree_roots[root_i];

		Vector<FBXNodeIndex> subtree_nodes;
		subtree_set.get_members(subtree_nodes, subtree_root);

		for (int subtree_i = 0; subtree_i < subtree_nodes.size(); ++subtree_i) {
			Ref<FBXNode> node = nodes[subtree_nodes[subtree_i]];
			node->joint = true;
			// Add the joint to the skeletons joints
			p_skeleton->joints.push_back(subtree_nodes[subtree_i]);
		}
	}

	return OK;
}
Error FBXSkinUtility::_determine_skeleton_roots(
		Vector<Ref<FBXNode>> &nodes,
		Vector<Ref<FBXSkeleton>> &skeletons,
		const FBXSkeletonIndex p_skel_i) {
	DisjointSet<FBXNodeIndex> disjoint_set;

	for (FBXNodeIndex i = 0; i < nodes.size(); ++i) {
		const Ref<FBXNode> node = nodes[i];

		if (node->skeleton != p_skel_i) {
			continue;
		}

		disjoint_set.insert(i);

		if (node->parent >= 0 && nodes[node->parent]->skeleton == p_skel_i) {
			disjoint_set.create_union(node->parent, i);
		}
	}

	Ref<FBXSkeleton> skeleton = skeletons.write[p_skel_i];

	Vector<FBXNodeIndex> representatives;
	disjoint_set.get_representatives(representatives);

	Vector<FBXNodeIndex> roots;

	for (int i = 0; i < representatives.size(); ++i) {
		Vector<FBXNodeIndex> set;
		disjoint_set.get_members(set, representatives[i]);
		const FBXNodeIndex root = _find_highest_node(nodes, set);
		ERR_FAIL_COND_V(root < 0, FAILED);
		roots.push_back(root);
	}

	roots.sort();

	skeleton->roots = roots;

	if (roots.size() == 0) {
		return FAILED;
	} else if (roots.size() == 1) {
		return OK;
	}

	// Check that the subtrees have the same parent root
	const FBXNodeIndex parent = nodes[roots[0]]->parent;
	for (int i = 1; i < roots.size(); ++i) {
		if (nodes[roots[i]]->parent != parent) {
			return FAILED;
		}
	}

	return OK;
}
Error FBXSkinUtility::_create_skeletons(
		HashSet<String> &unique_names,
		Vector<Ref<FBXSkin>> &skins,
		Vector<Ref<FBXNode>> &nodes,
		HashMap<ObjectID, FBXSkeletonIndex> &skeleton3d_to_fbx_skeleton,
		Vector<Ref<FBXSkeleton>> &skeletons,
		HashMap<FBXNodeIndex, Node *> &scene_nodes) {
	for (FBXSkeletonIndex skel_i = 0; skel_i < skeletons.size(); ++skel_i) {
		Ref<FBXSkeleton> fbx_skeleton = skeletons.write[skel_i];

		Skeleton3D *skeleton = memnew(Skeleton3D);
		fbx_skeleton->godot_skeleton = skeleton;
		skeleton3d_to_fbx_skeleton[skeleton->get_instance_id()] = skel_i;

		// Make a unique name, no gltf node represents this skeleton
		skeleton->set_name("Skeleton3D");

		List<FBXNodeIndex> bones;

		for (int i = 0; i < fbx_skeleton->roots.size(); ++i) {
			bones.push_back(fbx_skeleton->roots[i]);
		}

		// Make the skeleton creation deterministic by going through the roots in
		// a sorted order, and DEPTH FIRST
		bones.sort();

		while (!bones.is_empty()) {
			const FBXNodeIndex node_i = bones.front()->get();
			bones.pop_front();

			Ref<FBXNode> node = nodes[node_i];
			ERR_FAIL_COND_V(node->skeleton != skel_i, FAILED);

			{ // Add all child nodes to the stack (deterministically)
				Vector<FBXNodeIndex> child_nodes;
				for (int i = 0; i < node->children.size(); ++i) {
					const FBXNodeIndex child_i = node->children[i];
					if (nodes[child_i]->skeleton == skel_i) {
						child_nodes.push_back(child_i);
					}
				}

				// Depth first insertion
				child_nodes.sort();
				for (int i = child_nodes.size() - 1; i >= 0; --i) {
					bones.push_front(child_nodes[i]);
				}
			}

			const int bone_index = skeleton->get_bone_count();

			if (node->get_name().is_empty()) {
				node->set_name("bone");
			}

			node->set_name(_gen_unique_bone_name(unique_names, node->get_name()));

			skeleton->add_bone(node->get_name());
			skeleton->set_bone_rest(bone_index, node->xform);
			skeleton->set_bone_pose_position(bone_index, node->position);
			skeleton->set_bone_pose_rotation(bone_index, node->rotation.normalized());
			skeleton->set_bone_pose_scale(bone_index, node->scale);

			if (node->parent >= 0 && nodes[node->parent]->skeleton == skel_i) {
				const int bone_parent = skeleton->find_bone(nodes[node->parent]->get_name());
				ERR_FAIL_COND_V(bone_parent < 0, FAILED);
				skeleton->set_bone_parent(bone_index, skeleton->find_bone(nodes[node->parent]->get_name()));
			}

			scene_nodes.insert(node_i, skeleton);
		}
	}

	ERR_FAIL_COND_V(_map_skin_joints_indices_to_skeleton_bone_indices(skins, skeletons, nodes), ERR_PARSE_ERROR);

	return OK;
}
Error FBXSkinUtility::_map_skin_joints_indices_to_skeleton_bone_indices(
		Vector<Ref<FBXSkin>> &skins,
		Vector<Ref<FBXSkeleton>> &skeletons,
		Vector<Ref<FBXNode>> &nodes) {
	for (FBXSkinIndex skin_i = 0; skin_i < skins.size(); ++skin_i) {
		Ref<FBXSkin> skin = skins.write[skin_i];

		Ref<FBXSkeleton> skeleton = skeletons[skin->skeleton];

		for (int joint_index = 0; joint_index < skin->joints_original.size(); ++joint_index) {
			const FBXNodeIndex node_i = skin->joints_original[joint_index];
			const Ref<FBXNode> node = nodes[node_i];

			const int bone_index = skeleton->godot_skeleton->find_bone(node->get_name());
			ERR_FAIL_COND_V(bone_index < 0, FAILED);

			skin->joint_i_to_bone_i.insert(joint_index, bone_index);
		}
	}

	return OK;
}
Error FBXSkinUtility::_create_skins(Vector<Ref<FBXSkin>> &skins, Vector<Ref<FBXNode>> &nodes, bool use_named_skin_binds, HashSet<String> &unique_names) {
	for (FBXSkinIndex skin_i = 0; skin_i < skins.size(); ++skin_i) {
		Ref<FBXSkin> gltf_skin = skins.write[skin_i];

		Ref<Skin> skin;
		skin.instantiate();

		// Some skins don't have IBM's! What absolute monsters!
		const bool has_ibms = !gltf_skin->inverse_binds.is_empty();

		for (int joint_i = 0; joint_i < gltf_skin->joints_original.size(); ++joint_i) {
			FBXNodeIndex node = gltf_skin->joints_original[joint_i];
			String bone_name = nodes[node]->get_name();

			Transform3D xform;
			if (has_ibms) {
				xform = gltf_skin->inverse_binds[joint_i];
			}

			if (use_named_skin_binds) {
				skin->add_named_bind(bone_name, xform);
			} else {
				int32_t bone_i = gltf_skin->joint_i_to_bone_i[joint_i];
				skin->add_bind(bone_i, xform);
			}
		}

		gltf_skin->godot_skin = skin;
	}

	// Purge the duplicates!
	_remove_duplicate_skins(skins);

	// Create unique names now, after removing duplicates
	for (FBXSkinIndex skin_i = 0; skin_i < skins.size(); ++skin_i) {
		Ref<Skin> skin = skins.write[skin_i]->godot_skin;
		if (skin->get_name().is_empty()) {
			// Make a unique name, no node represents this skin
			skin->set_name(FBXDocument::_gen_unique_name(unique_names, "Skin"));
		}
	}

	return OK;
}
bool FBXSkinUtility::_skins_are_same(const Ref<Skin> p_skin_a, const Ref<Skin> p_skin_b) {
	if (p_skin_a->get_bind_count() != p_skin_b->get_bind_count()) {
		return false;
	}

	for (int i = 0; i < p_skin_a->get_bind_count(); ++i) {
		if (p_skin_a->get_bind_bone(i) != p_skin_b->get_bind_bone(i)) {
			return false;
		}
		if (p_skin_a->get_bind_name(i) != p_skin_b->get_bind_name(i)) {
			return false;
		}

		Transform3D a_xform = p_skin_a->get_bind_pose(i);
		Transform3D b_xform = p_skin_b->get_bind_pose(i);

		if (a_xform != b_xform) {
			return false;
		}
	}

	return true;
}
void FBXSkinUtility::_remove_duplicate_skins(Vector<Ref<FBXSkin>> &r_skins) {
	for (int i = 0; i < r_skins.size(); ++i) {
		for (int j = i + 1; j < r_skins.size(); ++j) {
			const Ref<Skin> skin_i = r_skins[i]->godot_skin;
			const Ref<Skin> skin_j = r_skins[j]->godot_skin;

			if (_skins_are_same(skin_i, skin_j)) {
				// replace it and delete the old
				r_skins.write[j]->godot_skin = skin_i;
			}
		}
	}
}
String FBXSkinUtility::_gen_unique_bone_name(HashSet<String> unique_names, const String &p_name) {
	String s_name = _sanitize_bone_name(p_name);
	if (s_name.is_empty()) {
		s_name = "bone";
	}
	String u_name;
	int index = 1;
	while (true) {
		u_name = s_name;

		if (index > 1) {
			u_name += "_" + itos(index);
		}
		if (!unique_names.has(u_name)) {
			break;
		}
		index++;
	}

	unique_names.insert(u_name);

	return u_name;
}
Error FBXSkinUtility::asset_parse_skins(
		const Vector<FBXNodeIndex> &input_skin_indices,
		const Vector<Ref<FBXSkin>> &input_skins,
		const Vector<Ref<FBXNode>> &input_nodes,
		Vector<FBXNodeIndex> &output_skin_indices,
		Vector<Ref<FBXSkin>> &output_skins,
		HashMap<FBXNodeIndex, bool> &joint_mapping) {
	output_skin_indices.clear();
	output_skins.clear();
	joint_mapping.clear();

	for (int i = 0; i < input_skin_indices.size(); ++i) {
		FBXNodeIndex skin_index = input_skin_indices[i];
		if (skin_index >= 0 && skin_index < input_skins.size()) {
			output_skin_indices.push_back(skin_index);
			output_skins.push_back(input_skins[skin_index]);
			Ref<FBXSkin> skin = input_skins[skin_index];
			Vector<FBXNodeIndex> skin_joints = skin->get_joints();
			for (int j = 0; j < skin_joints.size(); ++j) {
				FBXNodeIndex joint_index = skin_joints[j];
				joint_mapping[joint_index] = true;
			}
		}
	}

	return OK;
}
String FBXSkinUtility::_sanitize_bone_name(const String &p_name) {
	String bone_name = p_name;
	bone_name = bone_name.replace(":", "_");
	bone_name = bone_name.replace("/", "_");
	return bone_name;
}
