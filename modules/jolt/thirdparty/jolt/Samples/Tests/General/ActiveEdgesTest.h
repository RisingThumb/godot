// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Tests/Test.h>

// This test checks if coplanar triangles are properly merged for the simulation in order to avoid collisions with invisible edges.
class ActiveEdgesTest : public Test
{
public:
	JPH_DECLARE_RTTI_VIRTUAL(JPH_NO_EXPORT, ActiveEdgesTest)

	// See: Test
	virtual void		Initialize() override;
};
