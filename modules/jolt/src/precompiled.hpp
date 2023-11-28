#pragma once

// NOLINTBEGIN(readability-duplicate-include): We need to be able to include certain headers
// multiple times when they're conditionally included through multiple preprocessor definitions that
// might not be perfectly mutually exclusive.

#ifdef _MSC_VER
// HACK(mihe): CMake's Visual Studio generator doesn't support system include paths
#pragma warning(push, 0)

// Pushing level 0 doesn't seem to disable the ones we've explicitly enabled
// C4245: conversion from 'type1' to 'type2', signed/unsigned mismatch
// C4365: conversion from 'type1' to 'type2', signed/unsigned mismatch
#pragma warning(disable : 4245 4365)
#endif // _MSC_VER


#ifdef GDJ_CONFIG_EDITOR
#endif // GDJ_CONFIG_EDITOR

#ifdef JPH_DEBUG_RENDERER

#endif // JPH_DEBUG_RENDERER

#include <Jolt/Jolt.h>

#include <Jolt/Core/Factory.h>
#include <Jolt/Core/FixedSizeFreeList.h>
#include <Jolt/Core/IssueReporting.h>
#include <Jolt/Core/JobSystemWithBarrier.h>
#include <Jolt/Core/TempAllocator.h>
#include <Jolt/Geometry/ConvexSupport.h>
#include <Jolt/Geometry/GJKClosestPoint.h>
#include <Jolt/Physics/Body/BodyCreationSettings.h>
#include <Jolt/Physics/Body/BodyID.h>
#include <Jolt/Physics/Collision/BroadPhase/BroadPhaseLayer.h>
#include <Jolt/Physics/Collision/BroadPhase/BroadPhaseQuery.h>
#include <Jolt/Physics/Collision/CastResult.h>
#include <Jolt/Physics/Collision/CollidePointResult.h>
#include <Jolt/Physics/Collision/CollideShape.h>
#include <Jolt/Physics/Collision/CollisionDispatch.h>
#include <Jolt/Physics/Collision/CollisionGroup.h>
#include <Jolt/Physics/Collision/ContactListener.h>
#include <Jolt/Physics/Collision/EstimateCollisionResponse.h>
#include <Jolt/Physics/Collision/GroupFilter.h>
#include <Jolt/Physics/Collision/NarrowPhaseQuery.h>
#include <Jolt/Physics/Collision/ObjectLayer.h>
#include <Jolt/Physics/Collision/RayCast.h>
#include <Jolt/Physics/Collision/Shape/BoxShape.h>
#include <Jolt/Physics/Collision/Shape/CapsuleShape.h>
#include <Jolt/Physics/Collision/Shape/ConvexHullShape.h>
#include <Jolt/Physics/Collision/Shape/CylinderShape.h>
#include <Jolt/Physics/Collision/Shape/HeightFieldShape.h>
#include <Jolt/Physics/Collision/Shape/MeshShape.h>
#include <Jolt/Physics/Collision/Shape/MutableCompoundShape.h>
#include <Jolt/Physics/Collision/Shape/OffsetCenterOfMassShape.h>
#include <Jolt/Physics/Collision/Shape/RotatedTranslatedShape.h>
#include <Jolt/Physics/Collision/Shape/ScaledShape.h>
#include <Jolt/Physics/Collision/Shape/SphereShape.h>
#include <Jolt/Physics/Collision/Shape/StaticCompoundShape.h>
#include <Jolt/Physics/Constraints/FixedConstraint.h>
#include <Jolt/Physics/Constraints/HingeConstraint.h>
#include <Jolt/Physics/Constraints/PointConstraint.h>
#include <Jolt/Physics/Constraints/SixDOFConstraint.h>
#include <Jolt/Physics/Constraints/SliderConstraint.h>
#include <Jolt/Physics/Constraints/SwingTwistConstraint.h>
#include <Jolt/Physics/PhysicsScene.h>
#include <Jolt/Physics/PhysicsSystem.h>
#include <Jolt/RegisterTypes.h>

#ifdef JPH_DEBUG_RENDERER

#include <Jolt/Renderer/DebugRenderer.h>

#endif // JPH_DEBUG_RENDERER

#include <algorithm>
#include <atomic>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <thread>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

#ifdef _MSC_VER
#pragma warning(pop)
#endif // _MSC_VER

#include "modules/jolt/src/containers/free_list.hpp"
#include "modules/jolt/src/containers/hash_map.hpp"
#include "modules/jolt/src/containers/hash_set.hpp"
#include "modules/jolt/src/containers/inline_vector.hpp"
#include "modules/jolt/src/containers/local_vector.hpp"
#include "modules/jolt/src/containers/rid_owner.hpp"
#include "modules/jolt/src/misc/bind_macros.hpp"
#include "modules/jolt/src/misc/error_macros.hpp"
#include "modules/jolt/src/misc/gdclass_macros.hpp"
#include "modules/jolt/src/misc/jolt_stream_wrappers.hpp"
#include "modules/jolt/src/misc/math.hpp"
#include "modules/jolt/src/misc/scope_guard.hpp"
#include "modules/jolt/src/misc/type_conversions.hpp"
#include "modules/jolt/src/misc/utility_functions.hpp"
#include "modules/jolt/src/misc/utility_macros.hpp"

// NOLINTEND(readability-duplicate-include)
