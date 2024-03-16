// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2023 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Jolt.h>

#include <Jolt/Physics/SoftBody/SoftBodyMotionProperties.h>
#include <Jolt/Physics/SoftBody/SoftBodyCreationSettings.h>
#include <Jolt/Physics/SoftBody/SoftBodyContactListener.h>
#include <Jolt/Physics/SoftBody/SoftBodyManifold.h>
#include <Jolt/Physics/PhysicsSystem.h>
#ifdef JPH_DEBUG_RENDERER
	#include <Jolt/Renderer/DebugRenderer.h>
#endif // JPH_DEBUG_RENDERER

JPH_NAMESPACE_BEGIN

void SoftBodyMotionProperties::CalculateMassAndInertia()
{
	MassProperties mp;

	for (const Vertex &v : mVertices)
		if (v.mInvMass > 0.0f)
		{
			Vec3 pos = v.mPosition;

			// Accumulate mass
			float mass = 1.0f / v.mInvMass;
			mp.mMass += mass;

			// Inertia tensor, diagonal
			// See equations https://en.wikipedia.org/wiki/Moment_of_inertia section 'Inertia Tensor'
			for (int i = 0; i < 3; ++i)
				mp.mInertia(i, i) += mass * (Square(pos[(i + 1) % 3]) + Square(pos[(i + 2) % 3]));

			// Inertia tensor off diagonal
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j)
					if (i != j)
						mp.mInertia(i, j) -= mass * pos[i] * pos[j];
		}
		else
		{
			// If one vertex is kinematic, the entire body will have infinite mass and inertia
			SetInverseMass(0.0f);
			SetInverseInertia(Vec3::sZero(), Quat::sIdentity());
			return;
		}

	SetMassProperties(EAllowedDOFs::All, mp);
}

void SoftBodyMotionProperties::Initialize(const SoftBodyCreationSettings &inSettings)
{
	// Store settings
	mSettings = inSettings.mSettings;
	mNumIterations = inSettings.mNumIterations;
	mPressure = inSettings.mPressure;
	mUpdatePosition = inSettings.mUpdatePosition;

	// Initialize vertices
	mVertices.resize(inSettings.mSettings->mVertices.size());
	Mat44 rotation = inSettings.mMakeRotationIdentity? Mat44::sRotation(inSettings.mRotation) : Mat44::sIdentity();
	for (Array<Vertex>::size_type v = 0, s = mVertices.size(); v < s; ++v)
	{
		const SoftBodySharedSettings::Vertex &in_vertex = inSettings.mSettings->mVertices[v];
		Vertex &out_vertex = mVertices[v];
		out_vertex.mPreviousPosition = out_vertex.mPosition = rotation * Vec3(in_vertex.mPosition);
		out_vertex.mVelocity = rotation.Multiply3x3(Vec3(in_vertex.mVelocity));
		out_vertex.mCollidingShapeIndex = -1;
		out_vertex.mHasContact = false;
		out_vertex.mLargestPenetration = -FLT_MAX;
		out_vertex.mInvMass = in_vertex.mInvMass;
		mLocalBounds.Encapsulate(out_vertex.mPosition);
	}

	// Allocate space for skinned vertices
	if (!inSettings.mSettings->mSkinnedConstraints.empty())
		mSkinState.resize(mVertices.size());

	// We don't know delta time yet, so we can't predict the bounds and use the local bounds as the predicted bounds
	mLocalPredictedBounds = mLocalBounds;

	CalculateMassAndInertia();
}

float SoftBodyMotionProperties::GetVolumeTimesSix() const
{
	float six_volume = 0.0f;
	for (const Face &f : mSettings->mFaces)
	{
		Vec3 x1 = mVertices[f.mVertex[0]].mPosition;
		Vec3 x2 = mVertices[f.mVertex[1]].mPosition;
		Vec3 x3 = mVertices[f.mVertex[2]].mPosition;
		six_volume += x1.Cross(x2).Dot(x3); // We pick zero as the origin as this is the center of the bounding box so should give good accuracy
	}
	return six_volume;
}

void SoftBodyMotionProperties::DetermineCollidingShapes(const SoftBodyUpdateContext &inContext, const PhysicsSystem &inSystem, const BodyLockInterface &inBodyLockInterface)
{
	JPH_PROFILE_FUNCTION();

	struct Collector : public CollideShapeBodyCollector
	{
									Collector(const SoftBodyUpdateContext &inContext, const PhysicsSystem &inSystem, const BodyLockInterface &inBodyLockInterface, Array<CollidingShape> &ioHits) :
										mContext(inContext),
										mInverseTransform(inContext.mCenterOfMassTransform.InversedRotationTranslation()),
										mBodyLockInterface(inBodyLockInterface),
										mCombineFriction(inSystem.GetCombineFriction()),
										mCombineRestitution(inSystem.GetCombineRestitution()),
										mHits(ioHits)
		{
		}

		virtual void				AddHit(const BodyID &inResult) override
		{
			BodyLockRead lock(mBodyLockInterface, inResult);
			if (lock.Succeeded())
			{
				const Body &soft_body = *mContext.mBody;
				const Body &body = lock.GetBody();
				if (body.IsRigidBody() // TODO: We should support soft body vs soft body
					&& soft_body.GetCollisionGroup().CanCollide(body.GetCollisionGroup()))
				{
					// Call the contact listener to see if we should accept this contact
					// If there is no contact listener then we can ignore the contact if the other body is a sensor
					SoftBodyContactSettings settings;
					settings.mIsSensor = body.IsSensor();
					if (mContext.mContactListener == nullptr? !settings.mIsSensor : mContext.mContactListener->OnSoftBodyContactValidate(soft_body, body, settings) == SoftBodyValidateResult::AcceptContact)
					{
						CollidingShape cs;
						cs.mCenterOfMassTransform = (mInverseTransform * body.GetCenterOfMassTransform()).ToMat44();
						cs.mShape = body.GetShape();
						cs.mBodyID = inResult;
						cs.mMotionType = body.GetMotionType();
						cs.mIsSensor = settings.mIsSensor;
						cs.mUpdateVelocities = false;
						cs.mFriction = mCombineFriction(soft_body, SubShapeID(), body, SubShapeID());
						cs.mRestitution = mCombineRestitution(soft_body, SubShapeID(), body, SubShapeID());
						if (cs.mMotionType == EMotionType::Dynamic)
						{
							const MotionProperties *mp = body.GetMotionProperties();
							cs.mInvMass = settings.mInvMassScale2 * mp->GetInverseMass();
							cs.mInvInertia = settings.mInvInertiaScale2 * mp->GetInverseInertiaForRotation(cs.mCenterOfMassTransform.GetRotation());
							cs.mSoftBodyInvMassScale = settings.mInvMassScale1;
							cs.mOriginalLinearVelocity = cs.mLinearVelocity = mInverseTransform.Multiply3x3(mp->GetLinearVelocity());
							cs.mOriginalAngularVelocity = cs.mAngularVelocity = mInverseTransform.Multiply3x3(mp->GetAngularVelocity());
						}
						mHits.push_back(cs);
					}
				}
			}
		}

	private:
		const SoftBodyUpdateContext &mContext;
		RMat44						mInverseTransform;
		const BodyLockInterface &	mBodyLockInterface;
		ContactConstraintManager::CombineFunction mCombineFriction;
		ContactConstraintManager::CombineFunction mCombineRestitution;
		Array<CollidingShape> &		mHits;
	};

	Collector collector(inContext, inSystem, inBodyLockInterface, mCollidingShapes);
	AABox bounds = mLocalBounds;
	bounds.Encapsulate(mLocalPredictedBounds);
	bounds = bounds.Transformed(inContext.mCenterOfMassTransform);
	bounds.ExpandBy(Vec3::sReplicate(mSettings->mVertexRadius));
	ObjectLayer layer = inContext.mBody->GetObjectLayer();
	DefaultBroadPhaseLayerFilter broadphase_layer_filter = inSystem.GetDefaultBroadPhaseLayerFilter(layer);
	DefaultObjectLayerFilter object_layer_filter = inSystem.GetDefaultLayerFilter(layer);
	inSystem.GetBroadPhaseQuery().CollideAABox(bounds, collector, broadphase_layer_filter, object_layer_filter);
}

void SoftBodyMotionProperties::DetermineCollisionPlanes(const SoftBodyUpdateContext &inContext, uint inVertexStart, uint inNumVertices)
{
	JPH_PROFILE_FUNCTION();

	// Generate collision planes
	for (const CollidingShape &cs : mCollidingShapes)
		cs.mShape->CollideSoftBodyVertices(cs.mCenterOfMassTransform, Vec3::sReplicate(1.0f), mVertices.data() + inVertexStart, inNumVertices, inContext.mDeltaTime, inContext.mDisplacementDueToGravity, int(&cs - mCollidingShapes.data()));
}

void SoftBodyMotionProperties::ApplyPressure(const SoftBodyUpdateContext &inContext)
{
	JPH_PROFILE_FUNCTION();

	float dt = inContext.mSubStepDeltaTime;
	float pressure_coefficient = mPressure;
	if (pressure_coefficient > 0.0f)
	{
		// Calculate total volume
		float six_volume = GetVolumeTimesSix();
		if (six_volume > 0.0f)
		{
			// Apply pressure
			// p = F / A = n R T / V (see https://en.wikipedia.org/wiki/Pressure)
			// Our pressure coefficient is n R T so the impulse is:
			// P = F dt = pressure_coefficient / V * A * dt
			float coefficient = pressure_coefficient * dt / six_volume; // Need to still multiply by 6 for the volume
			for (const Face &f : mSettings->mFaces)
			{
				Vec3 x1 = mVertices[f.mVertex[0]].mPosition;
				Vec3 x2 = mVertices[f.mVertex[1]].mPosition;
				Vec3 x3 = mVertices[f.mVertex[2]].mPosition;

				Vec3 impulse = coefficient * (x2 - x1).Cross(x3 - x1); // Area is half the cross product so need to still divide by 2
				for (uint32 i : f.mVertex)
				{
					Vertex &v = mVertices[i];
					v.mVelocity += v.mInvMass * impulse; // Want to divide by 3 because we spread over 3 vertices
				}
			}
		}
	}
}

void SoftBodyMotionProperties::IntegratePositions(const SoftBodyUpdateContext &inContext)
{
	JPH_PROFILE_FUNCTION();

	float dt = inContext.mSubStepDeltaTime;
	float linear_damping = max(0.0f, 1.0f - GetLinearDamping() * dt); // See: MotionProperties::ApplyForceTorqueAndDragInternal

	// Integrate
	Vec3 sub_step_gravity = inContext.mGravity * dt;
	for (Vertex &v : mVertices)
		if (v.mInvMass > 0.0f)
		{
			// Gravity
			v.mVelocity += sub_step_gravity;

			// Damping
			v.mVelocity *= linear_damping;

			// Integrate
			v.mPreviousPosition = v.mPosition;
			v.mPosition += v.mVelocity * dt;
		}
		else
		{
			// Integrate
			v.mPreviousPosition = v.mPosition;
			v.mPosition += v.mVelocity * dt;
		}
}

void SoftBodyMotionProperties::ApplyVolumeConstraints(const SoftBodyUpdateContext &inContext)
{
	JPH_PROFILE_FUNCTION();

	float inv_dt_sq = 1.0f / Square(inContext.mSubStepDeltaTime);

	// Satisfy volume constraints
	for (const Volume &v : mSettings->mVolumeConstraints)
	{
		Vertex &v1 = mVertices[v.mVertex[0]];
		Vertex &v2 = mVertices[v.mVertex[1]];
		Vertex &v3 = mVertices[v.mVertex[2]];
		Vertex &v4 = mVertices[v.mVertex[3]];

		Vec3 x1 = v1.mPosition;
		Vec3 x2 = v2.mPosition;
		Vec3 x3 = v3.mPosition;
		Vec3 x4 = v4.mPosition;

		// Calculate constraint equation
		Vec3 x1x2 = x2 - x1;
		Vec3 x1x3 = x3 - x1;
		Vec3 x1x4 = x4 - x1;
		float c = abs(x1x2.Cross(x1x3).Dot(x1x4)) - v.mSixRestVolume;

		// Calculate gradient of constraint equation
		Vec3 d1c = (x4 - x2).Cross(x3 - x2);
		Vec3 d2c = x1x3.Cross(x1x4);
		Vec3 d3c = x1x4.Cross(x1x2);
		Vec3 d4c = x1x2.Cross(x1x3);

		float w1 = v1.mInvMass;
		float w2 = v2.mInvMass;
		float w3 = v3.mInvMass;
		float w4 = v4.mInvMass;
		JPH_ASSERT(w1 > 0.0f || w2 > 0.0f || w3 > 0.0f || w4 > 0.0f);

		// Apply correction
		float lambda = -c / (w1 * d1c.LengthSq() + w2 * d2c.LengthSq() + w3 * d3c.LengthSq() + w4 * d4c.LengthSq() + v.mCompliance * inv_dt_sq);
		v1.mPosition += lambda * w1 * d1c;
		v2.mPosition += lambda * w2 * d2c;
		v3.mPosition += lambda * w3 * d3c;
		v4.mPosition += lambda * w4 * d4c;
	}
}

void SoftBodyMotionProperties::ApplySkinConstraints([[maybe_unused]] const SoftBodyUpdateContext &inContext)
{
	// Early out if nothing to do
	if (mSettings->mSkinnedConstraints.empty())
		return;

	JPH_ASSERT(mSkinStateTransform == inContext.mCenterOfMassTransform, "Skinning state is stale, artifacts will show!");

	// Apply the constraints
	Vertex *vertices = mVertices.data();
	const SkinState *skin_states = mSkinState.data();
	for (const Skinned &s : mSettings->mSkinnedConstraints)
	{
		Vertex &vertex = vertices[s.mVertex];
		const SkinState &skin_state = skin_states[s.mVertex];
		if (s.mMaxDistance > 0.0f)
		{
			// Move vertex if it violated the back stop
			if (s.mBackStopDistance < s.mMaxDistance)
			{
				// Center of the back stop sphere
				Vec3 center = skin_state.mPosition - skin_state.mNormal * (s.mBackStopDistance + s.mBackStopRadius);

				// Check if we're inside the back stop sphere
				Vec3 delta = vertex.mPosition - center;
				float delta_len_sq = delta.LengthSq();
				if (delta_len_sq < Square(s.mBackStopRadius))
				{
					// Push the vertex to the surface of the back stop sphere
					float delta_len = sqrt(delta_len_sq);
					vertex.mPosition = delta_len > 0.0f?
						center + delta * (s.mBackStopRadius / delta_len)
						: center + skin_state.mNormal * s.mBackStopRadius;
				}
			}

			// Clamp vertex distance to max distance from skinned position
			if (s.mMaxDistance < FLT_MAX)
			{
				Vec3 delta = vertex.mPosition - skin_state.mPosition;
				float delta_len_sq = delta.LengthSq();
				float max_distance_sq = Square(s.mMaxDistance);
				if (delta_len_sq > max_distance_sq)
					vertex.mPosition = skin_state.mPosition + delta * sqrt(max_distance_sq / delta_len_sq);
			}
		}
		else
		{
			// Kinematic: Just update the vertex position
			vertex.mPosition = skin_state.mPosition;
		}
	}
}

void SoftBodyMotionProperties::ApplyEdgeConstraints(const SoftBodyUpdateContext &inContext, uint inStartIndex, uint inEndIndex)
{
	JPH_PROFILE_FUNCTION();

	float inv_dt_sq = 1.0f / Square(inContext.mSubStepDeltaTime);

	// Satisfy edge constraints
	const Array<Edge> &edge_constraints = mSettings->mEdgeConstraints;
	for (uint i = inStartIndex; i < inEndIndex; ++i)
	{
		const Edge &e = edge_constraints[i];
		Vertex &v0 = mVertices[e.mVertex[0]];
		Vertex &v1 = mVertices[e.mVertex[1]];
		JPH_ASSERT(v0.mInvMass > 0.0f || v1.mInvMass > 0.0f);

		// Calculate current length
		Vec3 delta = v1.mPosition - v0.mPosition;
		float length = delta.Length();
		if (length > 0.0f)
		{
			// Apply correction
			Vec3 correction = delta * (length - e.mRestLength) / (length * (v0.mInvMass + v1.mInvMass + e.mCompliance * inv_dt_sq));
			v0.mPosition += v0.mInvMass * correction;
			v1.mPosition -= v1.mInvMass * correction;
		}
	}
}

void SoftBodyMotionProperties::ApplyLRAConstraints()
{
	JPH_PROFILE_FUNCTION();

	// Satisfy LRA constraints
	Vertex *vertices = mVertices.data();
	for (const LRA &lra : mSettings->mLRAConstraints)
	{
		JPH_ASSERT(lra.mVertex[0] < mVertices.size());
		JPH_ASSERT(lra.mVertex[1] < mVertices.size());
		const Vertex &vertex0 = vertices[lra.mVertex[0]];
		Vertex &vertex1 = vertices[lra.mVertex[1]];

		Vec3 delta = vertex1.mPosition - vertex0.mPosition;
		float delta_len_sq = delta.LengthSq();
		if (delta_len_sq > Square(lra.mMaxDistance))
			vertex1.mPosition = vertex0.mPosition + delta * lra.mMaxDistance / sqrt(delta_len_sq);
	}
}

void SoftBodyMotionProperties::ApplyCollisionConstraintsAndUpdateVelocities(const SoftBodyUpdateContext &inContext)
{
	JPH_PROFILE_FUNCTION();

	float dt = inContext.mSubStepDeltaTime;
	float restitution_treshold = -2.0f * inContext.mGravity.Length() * dt;
	float vertex_radius = mSettings->mVertexRadius;
	for (Vertex &v : mVertices)
		if (v.mInvMass > 0.0f)
		{
			// Remember previous velocity for restitution calculations
			Vec3 prev_v = v.mVelocity;

			// XPBD velocity update
			v.mVelocity = (v.mPosition - v.mPreviousPosition) / dt;

			// Satisfy collision constraint
			if (v.mCollidingShapeIndex >= 0)
			{
				// Check if there is a collision
				float projected_distance = -v.mCollisionPlane.SignedDistance(v.mPosition) + vertex_radius;
				if (projected_distance > 0.0f)
				{
					// Remember that there was a collision
					v.mHasContact = true;
					mHasContact = true;

					// Sensors should not have a collision response
					CollidingShape &cs = mCollidingShapes[v.mCollidingShapeIndex];
					if (!cs.mIsSensor)
					{
						// Note that we already calculated the velocity, so this does not affect the velocity (next iteration starts by setting previous position to current position)
						Vec3 contact_normal = v.mCollisionPlane.GetNormal();
						v.mPosition += contact_normal * projected_distance;

						// Apply friction as described in Detailed Rigid Body Simulation with Extended Position Based Dynamics - Matthias Muller et al.
						// See section 3.6:
						// Inverse mass: w1 = 1 / m1, w2 = 1 / m2 + (r2 x n)^T I^-1 (r2 x n) = 0 for a static object
						// r2 are the contact point relative to the center of mass of body 2
						// Lagrange multiplier for contact: lambda = -c / (w1 + w2)
						// Where c is the constraint equation (the distance to the plane, negative because penetrating)
						// Contact normal force: fn = lambda / dt^2
						// Delta velocity due to friction dv = -vt / |vt| * min(dt * friction * fn * (w1 + w2), |vt|) = -vt * min(-friction * c / (|vt| * dt), 1)
						// Note that I think there is an error in the paper, I added a mass term, see: https://github.com/matthias-research/pages/issues/29
						// Relative velocity: vr = v1 - v2 - omega2 x r2
						// Normal velocity: vn = vr . contact_normal
						// Tangential velocity: vt = vr - contact_normal * vn
						// Impulse: p = dv / (w1 + w2)
						// Changes in particle velocities:
						// v1 = v1 + p / m1
						// v2 = v2 - p / m2 (no change when colliding with a static body)
						// w2 = w2 - I^-1 (r2 x p) (no change when colliding with a static body)
						if (cs.mMotionType == EMotionType::Dynamic)
						{
							// Calculate normal and tangential velocity (equation 30)
							Vec3 r2 = v.mPosition - cs.mCenterOfMassTransform.GetTranslation();
							Vec3 v2 = cs.GetPointVelocity(r2);
							Vec3 relative_velocity = v.mVelocity - v2;
							Vec3 v_normal = contact_normal * contact_normal.Dot(relative_velocity);
							Vec3 v_tangential = relative_velocity - v_normal;
							float v_tangential_length = v_tangential.Length();

							// Calculate resulting inverse mass of vertex
							float vertex_inv_mass = cs.mSoftBodyInvMassScale * v.mInvMass;

							// Calculate inverse effective mass
							Vec3 r2_cross_n = r2.Cross(contact_normal);
							float w2 = cs.mInvMass + r2_cross_n.Dot(cs.mInvInertia * r2_cross_n);
							float w1_plus_w2 = vertex_inv_mass + w2;

							// Calculate delta relative velocity due to friction (modified equation 31)
							Vec3 dv;
							if (v_tangential_length > 0.0f)
								dv = v_tangential * min(cs.mFriction * projected_distance / (v_tangential_length * dt), 1.0f);
							else
								dv = Vec3::sZero();

							// Calculate delta relative velocity due to restitution (equation 35)
							dv += v_normal;
							float prev_v_normal = (prev_v - v2).Dot(contact_normal);
							if (prev_v_normal < restitution_treshold)
								dv += cs.mRestitution * prev_v_normal * contact_normal;

							// Calculate impulse
							Vec3 p = dv / w1_plus_w2;

							// Apply impulse to particle
							v.mVelocity -= p * vertex_inv_mass;

							// Apply impulse to rigid body
							cs.mLinearVelocity += p * cs.mInvMass;
							cs.mAngularVelocity += cs.mInvInertia * r2.Cross(p);

							// Mark that the velocities of the body we hit need to be updated
							cs.mUpdateVelocities = true;
						}
						else
						{
							// Body is not movable, equations are simpler

							// Calculate normal and tangential velocity (equation 30)
							Vec3 v_normal = contact_normal * contact_normal.Dot(v.mVelocity);
							Vec3 v_tangential = v.mVelocity - v_normal;
							float v_tangential_length = v_tangential.Length();

							// Apply friction (modified equation 31)
							if (v_tangential_length > 0.0f)
								v.mVelocity -= v_tangential * min(cs.mFriction * projected_distance / (v_tangential_length * dt), 1.0f);

							// Apply restitution (equation 35)
							v.mVelocity -= v_normal;
							float prev_v_normal = prev_v.Dot(contact_normal);
							if (prev_v_normal < restitution_treshold)
								v.mVelocity -= cs.mRestitution * prev_v_normal * contact_normal;
						}
					}
				}
			}
		}
}

void SoftBodyMotionProperties::UpdateSoftBodyState(SoftBodyUpdateContext &ioContext, const PhysicsSettings &inPhysicsSettings)
{
	JPH_PROFILE_FUNCTION();

	// Contact callback
	if (mHasContact && ioContext.mContactListener != nullptr)
		ioContext.mContactListener->OnSoftBodyContactAdded(*ioContext.mBody, SoftBodyManifold(this));

	// Loop through vertices once more to update the global state
	float dt = ioContext.mDeltaTime;
	float max_linear_velocity_sq = Square(GetMaxLinearVelocity());
	float max_v_sq = 0.0f;
	Vec3 linear_velocity = Vec3::sZero(), angular_velocity = Vec3::sZero();
	mLocalPredictedBounds = mLocalBounds = { };
	mHasContact = false;
	for (Vertex &v : mVertices)
	{
		// Calculate max square velocity
		float v_sq = v.mVelocity.LengthSq();
		max_v_sq = max(max_v_sq, v_sq);

		// Clamp if velocity is too high
		if (v_sq > max_linear_velocity_sq)
			v.mVelocity *= sqrt(max_linear_velocity_sq / v_sq);

		// Calculate local linear/angular velocity
		linear_velocity += v.mVelocity;
		angular_velocity += v.mPosition.Cross(v.mVelocity);

		// Update local bounding box
		mLocalBounds.Encapsulate(v.mPosition);

		// Create predicted position for the next frame in order to detect collisions before they happen
		mLocalPredictedBounds.Encapsulate(v.mPosition + v.mVelocity * dt + ioContext.mDisplacementDueToGravity);

		// Reset collision data for the next iteration
		v.mCollidingShapeIndex = -1;
		v.mHasContact = false;
		v.mLargestPenetration = -FLT_MAX;
	}

	// Calculate linear/angular velocity of the body by averaging all vertices and bringing the value to world space
	float num_vertices_divider = float(max(int(mVertices.size()), 1));
	SetLinearVelocity(ioContext.mCenterOfMassTransform.Multiply3x3(linear_velocity / num_vertices_divider));
	SetAngularVelocity(ioContext.mCenterOfMassTransform.Multiply3x3(angular_velocity / num_vertices_divider));

	if (mUpdatePosition)
	{
		// Shift the body so that the position is the center of the local bounds
		Vec3 delta = mLocalBounds.GetCenter();
		ioContext.mDeltaPosition = ioContext.mCenterOfMassTransform.Multiply3x3(delta);
		for (Vertex &v : mVertices)
			v.mPosition -= delta;

		// Offset bounds to match new position
		mLocalBounds.Translate(-delta);
		mLocalPredictedBounds.Translate(-delta);
	}
	else
		ioContext.mDeltaPosition = Vec3::sZero();

	// Test if we should go to sleep
	if (GetAllowSleeping())
	{
		if (max_v_sq > inPhysicsSettings.mPointVelocitySleepThreshold)
		{
			ResetSleepTestTimer();
			ioContext.mCanSleep = ECanSleep::CannotSleep;
		}
		else
			ioContext.mCanSleep = AccumulateSleepTime(dt, inPhysicsSettings.mTimeBeforeSleep);
	}
	else
		ioContext.mCanSleep = ECanSleep::CannotSleep;
}

void SoftBodyMotionProperties::UpdateRigidBodyVelocities(const SoftBodyUpdateContext &inContext, BodyInterface &inBodyInterface)
{
	JPH_PROFILE_FUNCTION();

	// Write back velocity deltas
	for (const CollidingShape &cs : mCollidingShapes)
		if (cs.mUpdateVelocities)
			inBodyInterface.AddLinearAndAngularVelocity(cs.mBodyID, inContext.mCenterOfMassTransform.Multiply3x3(cs.mLinearVelocity - cs.mOriginalLinearVelocity), inContext.mCenterOfMassTransform.Multiply3x3(cs.mAngularVelocity - cs.mOriginalAngularVelocity));

	// Clear colliding shapes to avoid hanging on to references to shapes
	mCollidingShapes.clear();
}

void SoftBodyMotionProperties::InitializeUpdateContext(float inDeltaTime, Body &inSoftBody, const PhysicsSystem &inSystem, SoftBodyUpdateContext &ioContext)
{
	JPH_PROFILE_FUNCTION();

	// Store body
	ioContext.mBody = &inSoftBody;
	ioContext.mMotionProperties = this;
	ioContext.mContactListener = inSystem.GetSoftBodyContactListener();

	// Convert gravity to local space
	ioContext.mCenterOfMassTransform = inSoftBody.GetCenterOfMassTransform();
	ioContext.mGravity = ioContext.mCenterOfMassTransform.Multiply3x3Transposed(GetGravityFactor() * inSystem.GetGravity());

	// Calculate delta time for sub step
	ioContext.mDeltaTime = inDeltaTime;
	ioContext.mSubStepDeltaTime = inDeltaTime / mNumIterations;

	// Calculate total displacement we'll have due to gravity over all sub steps
	// The total displacement as produced by our integrator can be written as: Sum(i * g * dt^2, i = 0..mNumIterations).
	// This is bigger than 0.5 * g * dt^2 because we first increment the velocity and then update the position
	// Using Sum(i, i = 0..n) = n * (n + 1) / 2 we can write this as:
	ioContext.mDisplacementDueToGravity = (0.5f * mNumIterations * (mNumIterations + 1) * Square(ioContext.mSubStepDeltaTime)) * ioContext.mGravity;
}

void SoftBodyMotionProperties::StartNextIteration(const SoftBodyUpdateContext &ioContext)
{
	ApplyPressure(ioContext);

	IntegratePositions(ioContext);

	ApplyVolumeConstraints(ioContext);
}

SoftBodyMotionProperties::EStatus SoftBodyMotionProperties::ParallelDetermineCollisionPlanes(SoftBodyUpdateContext &ioContext)
{
	// Do a relaxed read first to see if there is any work to do (this prevents us from doing expensive atomic operations and also prevents us from continuously incrementing the counter and overflowing it)
	uint num_vertices = (uint)mVertices.size();
	if (ioContext.mNextCollisionVertex.load(memory_order_relaxed) < num_vertices)
	{
		// Fetch next batch of vertices to process
		uint next_vertex = ioContext.mNextCollisionVertex.fetch_add(SoftBodyUpdateContext::cVertexCollisionBatch, memory_order_acquire);
		if (next_vertex < num_vertices)
		{
			// Process collision planes
			uint num_vertices_to_process = min(SoftBodyUpdateContext::cVertexCollisionBatch, num_vertices - next_vertex);
			DetermineCollisionPlanes(ioContext, next_vertex, num_vertices_to_process);
			uint vertices_processed = ioContext.mNumCollisionVerticesProcessed.fetch_add(SoftBodyUpdateContext::cVertexCollisionBatch, memory_order_release) + num_vertices_to_process;
			if (vertices_processed >= num_vertices)
			{
				// Start the first iteration
				JPH_IF_ENABLE_ASSERTS(uint iteration =) ioContext.mNextIteration.fetch_add(1, memory_order_relaxed);
				JPH_ASSERT(iteration == 0);
				StartNextIteration(ioContext);
				ioContext.mState.store(SoftBodyUpdateContext::EState::ApplyEdgeConstraints, memory_order_release);
			}
			return EStatus::DidWork;
		}
	}
	return EStatus::NoWork;
}

SoftBodyMotionProperties::EStatus SoftBodyMotionProperties::ParallelApplyEdgeConstraints(SoftBodyUpdateContext &ioContext, const PhysicsSettings &inPhysicsSettings)
{
	// Do a relaxed read first to see if there is any work to do (this prevents us from doing expensive atomic operations and also prevents us from continuously incrementing the counter and overflowing it)
	uint num_groups = (uint)mSettings->mEdgeGroupEndIndices.size();
	JPH_ASSERT(num_groups > 0, "SoftBodySharedSettings::Optimize should have been called!");
	uint32 edge_group, edge_start_idx;
	SoftBodyUpdateContext::sGetEdgeGroupAndStartIdx(ioContext.mNextEdgeConstraint.load(memory_order_relaxed), edge_group, edge_start_idx);
	if (edge_group < num_groups)
	{
		uint edge_group_size = mSettings->GetEdgeGroupSize(edge_group);
		if (edge_start_idx < edge_group_size || edge_group_size == 0)
		{
			// Fetch the next batch of edges to process
			uint64 next_edge_batch = ioContext.mNextEdgeConstraint.fetch_add(SoftBodyUpdateContext::cEdgeConstraintBatch, memory_order_acquire);
			SoftBodyUpdateContext::sGetEdgeGroupAndStartIdx(next_edge_batch, edge_group, edge_start_idx);
			if (edge_group < num_groups)
			{
				bool non_parallel_group = edge_group == num_groups - 1; // Last group is the non-parallel group and goes as a whole
				edge_group_size = mSettings->GetEdgeGroupSize(edge_group);
				if (non_parallel_group? edge_start_idx == 0 : edge_start_idx < edge_group_size)
				{
					// Process edges
					uint num_edges_to_process = non_parallel_group? edge_group_size : min(SoftBodyUpdateContext::cEdgeConstraintBatch, edge_group_size - edge_start_idx);
					if (edge_group > 0)
						edge_start_idx += mSettings->mEdgeGroupEndIndices[edge_group - 1];
					ApplyEdgeConstraints(ioContext, edge_start_idx, edge_start_idx + num_edges_to_process);

					// Test if we're at the end of this group
					uint edge_constraints_processed = ioContext.mNumEdgeConstraintsProcessed.fetch_add(num_edges_to_process, memory_order_relaxed) + num_edges_to_process;
					if (edge_constraints_processed >= edge_group_size)
					{
						// Non parallel group is the last group (which is also the only group that can be empty)
						if (non_parallel_group || mSettings->GetEdgeGroupSize(edge_group + 1) == 0)
						{
							// Finish the iteration
							ApplyLRAConstraints();

							ApplyCollisionConstraintsAndUpdateVelocities(ioContext);

							ApplySkinConstraints(ioContext);

							uint iteration = ioContext.mNextIteration.fetch_add(1, memory_order_relaxed);
							if (iteration < mNumIterations)
							{
								// Start a new iteration
								StartNextIteration(ioContext);

								// Reset next edge to process
								ioContext.mNumEdgeConstraintsProcessed.store(0, memory_order_relaxed);
								ioContext.mNextEdgeConstraint.store(0, memory_order_release);
							}
							else
							{
								// On final iteration we update the state
								UpdateSoftBodyState(ioContext, inPhysicsSettings);

								ioContext.mState.store(SoftBodyUpdateContext::EState::Done, memory_order_release);
								return EStatus::Done;
							}
						}
						else
						{
							// Next group
							ioContext.mNumEdgeConstraintsProcessed.store(0, memory_order_relaxed);
							ioContext.mNextEdgeConstraint.store(SoftBodyUpdateContext::sGetEdgeGroupStart(edge_group + 1), memory_order_release);
						}
					}
					return EStatus::DidWork;
				}
			}
		}
	}
	return EStatus::NoWork;
}

SoftBodyMotionProperties::EStatus SoftBodyMotionProperties::ParallelUpdate(SoftBodyUpdateContext &ioContext, const PhysicsSettings &inPhysicsSettings)
{
	switch (ioContext.mState.load(memory_order_relaxed))
	{
	case SoftBodyUpdateContext::EState::DetermineCollisionPlanes:
		return ParallelDetermineCollisionPlanes(ioContext);

	case SoftBodyUpdateContext::EState::ApplyEdgeConstraints:
		return ParallelApplyEdgeConstraints(ioContext, inPhysicsSettings);

	case SoftBodyUpdateContext::EState::Done:
		return EStatus::Done;

	default:
		JPH_ASSERT(false);
		return EStatus::NoWork;
	}
}

void SoftBodyMotionProperties::SkinVertices(RMat44Arg inCenterOfMassTransform, const Mat44 *inJointMatrices, [[maybe_unused]] uint inNumJoints, bool inHardSkinAll, TempAllocator &ioTempAllocator)
{
	// Calculate the skin matrices
	uint num_skin_matrices = uint(mSettings->mInvBindMatrices.size());
	uint skin_matrices_size = num_skin_matrices * sizeof(Mat44);
	Mat44 *skin_matrices = (Mat44 *)ioTempAllocator.Allocate(skin_matrices_size);
	const Mat44 *skin_matrices_end = skin_matrices + num_skin_matrices;
	const InvBind *inv_bind_matrix = mSettings->mInvBindMatrices.data();
	for (Mat44 *s = skin_matrices; s < skin_matrices_end; ++s, ++inv_bind_matrix)
		*s = inJointMatrices[inv_bind_matrix->mJointIndex] * inv_bind_matrix->mInvBind;

	// Skin the vertices
	mSkinStateTransform = inCenterOfMassTransform;
	JPH_IF_ENABLE_ASSERTS(uint num_vertices = uint(mSettings->mVertices.size());)
	JPH_ASSERT(mSkinState.size() == num_vertices);
	const SoftBodySharedSettings::Vertex *in_vertices = mSettings->mVertices.data();
	for (const Skinned &s : mSettings->mSkinnedConstraints)
	{
		// Get bind pose
		JPH_ASSERT(s.mVertex < num_vertices);
		Vec3 bind_pos = Vec3::sLoadFloat3Unsafe(in_vertices[s.mVertex].mPosition);

		// Skin vertex
		Vec3 pos = Vec3::sZero();
		for (const SkinWeight &w : s.mWeights)
		{
			// We assume that the first zero weight is the end of the list
			if (w.mWeight == 0.0f)
				break;

			JPH_ASSERT(w.mInvBindIndex < num_skin_matrices);
			pos += w.mWeight * (skin_matrices[w.mInvBindIndex] * bind_pos);
		}
		mSkinState[s.mVertex].mPosition = pos;
	}

	// Calculate the normals
	for (const Skinned &s : mSettings->mSkinnedConstraints)
	{
		Vec3 normal = Vec3::sZero();
		uint32 num_faces = s.mNormalInfo >> 24;
		if (num_faces > 0)
		{
			// Calculate normal
			const uint32 *f = &mSettings->mSkinnedConstraintNormals[s.mNormalInfo & 0xffffff];
			const uint32 *f_end = f + num_faces;
			while (f < f_end)
			{
				const Face &face = mSettings->mFaces[*f];
				Vec3 v0 = mSkinState[face.mVertex[0]].mPosition;
				Vec3 v1 = mSkinState[face.mVertex[1]].mPosition;
				Vec3 v2 = mSkinState[face.mVertex[2]].mPosition;
				normal += (v1 - v0).Cross(v2 - v0).NormalizedOr(Vec3::sZero());
				++f;
			}
			normal = normal.NormalizedOr(Vec3::sZero());
		}
		mSkinState[s.mVertex].mNormal = normal;
	}

	ioTempAllocator.Free(skin_matrices, skin_matrices_size);

	if (inHardSkinAll)
	{
		// Hard skin all vertices and reset their velocities
		for (const Skinned &s : mSettings->mSkinnedConstraints)
		{
			Vertex &vertex = mVertices[s.mVertex];
			vertex.mPosition = mSkinState[s.mVertex].mPosition;
			vertex.mVelocity = Vec3::sZero();
		}
	}
}

void SoftBodyMotionProperties::CustomUpdate(float inDeltaTime, Body &ioSoftBody, PhysicsSystem &inSystem)
{
	JPH_PROFILE_FUNCTION();

	// Create update context
	SoftBodyUpdateContext context;
	InitializeUpdateContext(inDeltaTime, ioSoftBody, inSystem, context);

	// Determine bodies we're colliding with
	DetermineCollidingShapes(context, inSystem, inSystem.GetBodyLockInterface());

	// Call the internal update until it finishes
	EStatus status;
	const PhysicsSettings &settings = inSystem.GetPhysicsSettings();
	while ((status = ParallelUpdate(context, settings)) == EStatus::DidWork)
		continue;
	JPH_ASSERT(status == EStatus::Done);

	// Update the state of the bodies we've collided with
	UpdateRigidBodyVelocities(context, inSystem.GetBodyInterface());

	// Update position of the soft body
	if (mUpdatePosition)
		inSystem.GetBodyInterface().SetPosition(ioSoftBody.GetID(), ioSoftBody.GetPosition() + context.mDeltaPosition, EActivation::DontActivate);
}

#ifdef JPH_DEBUG_RENDERER

void SoftBodyMotionProperties::DrawVertices(DebugRenderer *inRenderer, RMat44Arg inCenterOfMassTransform) const
{
	for (const Vertex &v : mVertices)
		inRenderer->DrawMarker(inCenterOfMassTransform * v.mPosition, v.mInvMass > 0.0f? Color::sGreen : Color::sRed, 0.05f);
}

void SoftBodyMotionProperties::DrawVertexVelocities(DebugRenderer *inRenderer, RMat44Arg inCenterOfMassTransform) const
{
	for (const Vertex &v : mVertices)
		inRenderer->DrawArrow(inCenterOfMassTransform * v.mPosition, inCenterOfMassTransform * (v.mPosition + v.mVelocity), Color::sYellow, 0.01f);
}

void SoftBodyMotionProperties::DrawEdgeConstraints(DebugRenderer *inRenderer, RMat44Arg inCenterOfMassTransform) const
{
	for (const Edge &e : mSettings->mEdgeConstraints)
		inRenderer->DrawLine(inCenterOfMassTransform * mVertices[e.mVertex[0]].mPosition, inCenterOfMassTransform * mVertices[e.mVertex[1]].mPosition, Color::sWhite);
}

void SoftBodyMotionProperties::DrawVolumeConstraints(DebugRenderer *inRenderer, RMat44Arg inCenterOfMassTransform) const
{
	for (const Volume &v : mSettings->mVolumeConstraints)
	{
		RVec3 x1 = inCenterOfMassTransform * mVertices[v.mVertex[0]].mPosition;
		RVec3 x2 = inCenterOfMassTransform * mVertices[v.mVertex[1]].mPosition;
		RVec3 x3 = inCenterOfMassTransform * mVertices[v.mVertex[2]].mPosition;
		RVec3 x4 = inCenterOfMassTransform * mVertices[v.mVertex[3]].mPosition;

		inRenderer->DrawTriangle(x1, x3, x2, Color::sYellow, DebugRenderer::ECastShadow::On);
		inRenderer->DrawTriangle(x2, x3, x4, Color::sYellow, DebugRenderer::ECastShadow::On);
		inRenderer->DrawTriangle(x1, x4, x3, Color::sYellow, DebugRenderer::ECastShadow::On);
		inRenderer->DrawTriangle(x1, x2, x4, Color::sYellow, DebugRenderer::ECastShadow::On);
	}
}

void SoftBodyMotionProperties::DrawSkinConstraints(DebugRenderer *inRenderer, RMat44Arg inCenterOfMassTransform) const
{
	for (const Skinned &s : mSettings->mSkinnedConstraints)
	{
		const SkinState &skin_state = mSkinState[s.mVertex];
		DebugRenderer::sInstance->DrawArrow(mSkinStateTransform * skin_state.mPosition, mSkinStateTransform * (skin_state.mPosition + 0.1f * skin_state.mNormal), Color::sOrange, 0.01f);
		DebugRenderer::sInstance->DrawLine(mSkinStateTransform * skin_state.mPosition, inCenterOfMassTransform * mVertices[s.mVertex].mPosition, Color::sBlue);
	}
}

void SoftBodyMotionProperties::DrawLRAConstraints(DebugRenderer *inRenderer, RMat44Arg inCenterOfMassTransform) const
{
	for (const LRA &l : mSettings->mLRAConstraints)
		inRenderer->DrawLine(inCenterOfMassTransform * mVertices[l.mVertex[0]].mPosition, inCenterOfMassTransform * mVertices[l.mVertex[1]].mPosition, Color::sGrey);
}

void SoftBodyMotionProperties::DrawPredictedBounds(DebugRenderer *inRenderer, RMat44Arg inCenterOfMassTransform) const
{
	inRenderer->DrawWireBox(inCenterOfMassTransform, mLocalPredictedBounds, Color::sRed);
}

#endif // JPH_DEBUG_RENDERER

void SoftBodyMotionProperties::SaveState(StateRecorder &inStream) const
{
	MotionProperties::SaveState(inStream);

	for (const Vertex &v : mVertices)
	{
		inStream.Write(v.mPreviousPosition);
		inStream.Write(v.mPosition);
		inStream.Write(v.mVelocity);
	}

	inStream.Write(mLocalBounds.mMin);
	inStream.Write(mLocalBounds.mMax);
	inStream.Write(mLocalPredictedBounds.mMin);
	inStream.Write(mLocalPredictedBounds.mMax);
}

void SoftBodyMotionProperties::RestoreState(StateRecorder &inStream)
{
	MotionProperties::RestoreState(inStream);

	for (Vertex &v : mVertices)
	{
		inStream.Read(v.mPreviousPosition);
		inStream.Read(v.mPosition);
		inStream.Read(v.mVelocity);
	}

	inStream.Read(mLocalBounds.mMin);
	inStream.Read(mLocalBounds.mMax);
	inStream.Read(mLocalPredictedBounds.mMin);
	inStream.Read(mLocalPredictedBounds.mMax);
}

JPH_NAMESPACE_END
