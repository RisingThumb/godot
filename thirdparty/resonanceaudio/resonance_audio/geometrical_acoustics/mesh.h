/*
Copyright 2018 Google Inc. All Rights Reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS-IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

// Collection of data structures useful to define a mesh.

#ifndef RESONANCE_AUDIO_GEOMETRICAL_ACOUSTICS_MESH_H_
#define RESONANCE_AUDIO_GEOMETRICAL_ACOUSTICS_MESH_H_

namespace vraudio {

// A simple vertex data structure.
struct Vertex {
  float x;
  float y;
  float z;
};

// A simple triangle data structure defined as the 3 indices of vertices.
struct Triangle {
  int v0;
  int v1;
  int v2;
};

}  // namespace vraudio

#endif  // RESONANCE_AUDIO_GEOMETRICAL_ACOUSTICS_MESH_H_
