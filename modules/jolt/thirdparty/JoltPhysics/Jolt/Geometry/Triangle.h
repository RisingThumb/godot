// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

JPH_NAMESPACE_BEGIN

/// A simple triangle and its material
class Triangle
{
public:
	JPH_OVERRIDE_NEW_DELETE

	/// Constructor
					Triangle() = default;
					Triangle(const Float3 &inV1, const Float3 &inV2, const Float3 &inV3) : mV { inV1, inV2, inV3 } { }
					Triangle(const Float3 &inV1, const Float3 &inV2, const Float3 &inV3, uint32 inMaterialIndex) : Triangle(inV1, inV2, inV3) { mMaterialIndex = inMaterialIndex; }
					Triangle(Vec3Arg inV1, Vec3Arg inV2, Vec3Arg inV3) { inV1.StoreFloat3(&mV[0]); inV2.StoreFloat3(&mV[1]); inV3.StoreFloat3(&mV[2]); }

	/// Get center of triangle
	Vec3			GetCentroid() const
	{
		return (Vec3::sLoadFloat3Unsafe(mV[0]) + Vec3::sLoadFloat3Unsafe(mV[1]) + Vec3::sLoadFloat3Unsafe(mV[2])) * (1.0f / 3.0f);
	}

	/// Vertices
	Float3			mV[3];
	uint32			mMaterialIndex = 0;			///< Follows mV[3] so that we can read mV as 4 vectors
};

using TriangleList = Array<Triangle>;

JPH_NAMESPACE_END
