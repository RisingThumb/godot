#ifndef VRM_SECONDARY_GIZMO_H
#define VRM_SECONDARY_GIZMO_H

#include "modules/vrm/vrm_collidergroup.h"
#include "modules/vrm/vrm_secondary.h"
#include "scene/3d/mesh_instance_3d.h"
#include "scene/3d/skeleton_3d.h"
#include "scene/resources/immediate_mesh.h"
#include "scene/resources/material.h"
#include "vrm_spring_bone_logic.h"
#include "vrm_springbone.h"
#include "vrm_toplevel.h"

class SecondaryGizmo : public MeshInstance3D {
	GDCLASS(SecondaryGizmo, MeshInstance3D);

protected:
	static void _bind_methods();

public:
	VRMSecondary *secondary_node = nullptr;
	Ref<StandardMaterial3D> m;
	Ref<ImmediateMesh> mesh;

	~SecondaryGizmo();

	SecondaryGizmo(Node *parent = nullptr) {
		mesh.instantiate();
		secondary_node = cast_to<VRMSecondary>(parent);
		m->set_depth_draw_mode(BaseMaterial3D::DEPTH_DRAW_DISABLED);
		m->set_shading_mode(BaseMaterial3D::SHADING_MODE_UNSHADED);
		m->set_flag(StandardMaterial3D::FLAG_ALBEDO_FROM_VERTEX_COLOR, true);
		m->set_transparency(BaseMaterial3D::TRANSPARENCY_ALPHA);
	}

	void draw_in_editor() {
		mesh->clear_surfaces();
		if (secondary_node && Object::cast_to<VRMTopLevel>(secondary_node->get_parent())) {
			draw_spring_bones(Object::cast_to<VRMTopLevel>(secondary_node->get_parent())->get_gizmo_spring_bone_color());
			draw_collider_groups();
		}
	}

	void draw_in_game() {
		mesh->clear_surfaces();
		if (secondary_node && Object::cast_to<VRMTopLevel>(secondary_node->get_parent())) {
			draw_spring_bones(Object::cast_to<VRMTopLevel>(secondary_node->get_parent())->get_gizmo_spring_bone_color());
			draw_collider_groups();
		}
	}

	void draw_spring_bones(const Color &color);
	void draw_collider_groups();
	void draw_line(Vector3 begin_pos, Vector3 end_pos, Color color);
	void draw_sphere(Basis bas, Vector3 center, float radius, Color color);
};

#endif // VRM_SECONDARY_GIZMO_H