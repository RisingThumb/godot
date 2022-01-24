// Copyright 2019 The Draco Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#ifndef DRACO_SCENE_SCENE_H_
#define DRACO_SCENE_SCENE_H_

#include "draco/draco_features.h"

#ifdef DRACO_TRANSCODER_SUPPORTED
#include <memory>

#include "draco/animation/animation.h"
#include "draco/animation/skin.h"
#include "draco/mesh/mesh.h"
#include "draco/scene/mesh_group.h"
#include "draco/scene/scene_indices.h"
#include "draco/scene/scene_node.h"

namespace draco {

// Class used to hold all of the geometry to create a scene. A scene is
// comprised of one or more meshes, one or more scene nodes, one or more
// mesh groups, and a material library. The meshes are defined in their
// local space. A mesh group is a list of meshes. The scene nodes create
// a scene hierarchy to transform meshes in their local space into scene space.
// The material library contains all of the materials and textures used by the
// meshes in this scene.
class Scene {
 public:
  Scene() {}

  void Copy(const Scene &s);

  // Adds a Draco |mesh| to the scene. Returns the index to the stored mesh or
  // |kInvalidMeshIndex| if the mesh is a nullptr.
  MeshIndex AddMesh(std::unique_ptr<Mesh> mesh) {
    if (mesh == nullptr) {
      return kInvalidMeshIndex;
    }
    meshes_.push_back(std::move(mesh));
    return MeshIndex(meshes_.size() - 1);
  }

  // Removes base mesh and corresponding material at |index|, removes references
  // to removed base mesh and corresponding materials from mesh groups, and
  // updates references to remaining base meshes in mesh groups.
  Status RemoveMesh(MeshIndex index);

  // Returns the number of meshes in a scene before instancing is applied.
  int NumMeshes() const { return meshes_.size(); }

  // Returns a mesh in the scene before instancing is applied. The mesh
  // coordinates are local to the mesh.
  Mesh &GetMesh(MeshIndex index) { return *meshes_[index]; }
  const Mesh &GetMesh(MeshIndex index) const { return *meshes_[index]; }

  // Creates a mesh group and returns the index to the mesh group.
  MeshGroupIndex AddMeshGroup() {
    std::unique_ptr<MeshGroup> mesh(new MeshGroup());
    mesh_groups_.push_back(std::move(mesh));
    return MeshGroupIndex(mesh_groups_.size() - 1);
  }

  // Removes mesh group at |index|, invalidates references to removed mesh group
  // in scene nodes, and updates references to remaining mesh groups in scene
  // nodes.
  Status RemoveMeshGroup(MeshGroupIndex index);

  // Removes unused material at |index| and updates references to materials at
  // indices greater than |index|. Returns error status when |index| is out of
  // valid range and when material at |index| is used in the scene.
  Status RemoveMaterial(int index);

  // Returns the number of mesh groups in a scene.
  int NumMeshGroups() const { return mesh_groups_.size(); }

  // Returns a mesh group in the scene.
  MeshGroup *GetMeshGroup(MeshGroupIndex index) {
    return mesh_groups_[index].get();
  }
  const MeshGroup *GetMeshGroup(MeshGroupIndex index) const {
    return mesh_groups_[index].get();
  }

  // Creates a scene node and returns the index to the node.
  SceneNodeIndex AddNode() {
    std::unique_ptr<SceneNode> node(new SceneNode());
    nodes_.push_back(std::move(node));
    return SceneNodeIndex(nodes_.size() - 1);
  }

  // Returns the number of nodes in a scene.
  int NumNodes() const { return nodes_.size(); }

  // Returns a node in the scene.
  SceneNode *GetNode(SceneNodeIndex index) { return nodes_[index].get(); }
  const SceneNode *GetNode(SceneNodeIndex index) const {
    return nodes_[index].get();
  }

  // Returns the number of root node indices in a scene.
  int NumRootNodes() const { return root_node_indices_.size(); }
  SceneNodeIndex GetRootNodeIndex(int i) const { return root_node_indices_[i]; }
  void AddRootNodeIndex(SceneNodeIndex index) {
    root_node_indices_.push_back(index);
  }
  void SetRootNodeIndex(int i, SceneNodeIndex index) {
    root_node_indices_[i] = index;
  }

  const MaterialLibrary &GetMaterialLibrary() const {
    return material_library_;
  }
  MaterialLibrary &GetMaterialLibrary() { return material_library_; }

  // Creates an animation and returns the index to the animation.
  AnimationIndex AddAnimation() {
    std::unique_ptr<Animation> animation(new Animation());
    animations_.push_back(std::move(animation));
    return AnimationIndex(animations_.size() - 1);
  }

  // Returns the number of animations in a scene.
  int NumAnimations() const { return animations_.size(); }

  // Returns an animation in the scene.
  Animation *GetAnimation(AnimationIndex index) {
    return animations_[index].get();
  }
  const Animation *GetAnimation(AnimationIndex index) const {
    return animations_[index].get();
  }

  // Creates a skin and returns the index to the skin.
  SkinIndex AddSkin() {
    std::unique_ptr<Skin> skin(new Skin());
    skins_.push_back(std::move(skin));
    return SkinIndex(skins_.size() - 1);
  }

  // Returns the number of skins in a scene.
  int NumSkins() const { return skins_.size(); }

  // Returns a skin in the scene.
  Skin *GetSkin(SkinIndex index) { return skins_[index].get(); }
  const Skin *GetSkin(SkinIndex index) const { return skins_[index].get(); }

 private:
  IndexTypeVector<MeshIndex, std::unique_ptr<Mesh>> meshes_;
  IndexTypeVector<MeshGroupIndex, std::unique_ptr<MeshGroup>> mesh_groups_;
  IndexTypeVector<SceneNodeIndex, std::unique_ptr<SceneNode>> nodes_;
  std::vector<SceneNodeIndex> root_node_indices_;
  IndexTypeVector<AnimationIndex, std::unique_ptr<Animation>> animations_;
  IndexTypeVector<SkinIndex, std::unique_ptr<Skin>> skins_;

  // Materials used by this scene.
  MaterialLibrary material_library_;
};

}  // namespace draco

#endif  // DRACO_TRANSCODER_SUPPORTED
#endif  // DRACO_SCENE_SCENE_H_
