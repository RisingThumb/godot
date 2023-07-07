/**************************************************************************/
/*  crit_spring_damper.h                                                  */
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

#ifndef CRIT_SPRING_DAMPER_H
#define CRIT_SPRING_DAMPER_H

#include "core/config/project_settings.h"
#include "core/io/resource.h"
#include "core/object/class_db.h"
#include "core/object/method_bind.h"
#include "core/object/ref_counted.h"
#include "core/variant/dictionary.h"
#include "core/variant/variant.h"

#include <cmath>

struct CritDampSpring : public RefCounted {
	GDCLASS(CritDampSpring, RefCounted)

	static constexpr float Ln2 = 0.69314718056;

	static float inline square(float x) {
		return x * x;
	}

	static Vector3 damp_adjustment_exact(Vector3 g, float halflife, float dt, float eps = 1e-8) {
		float factor = 1.0 - fast_negexp((CritDampSpring::Ln2 * dt) / (halflife + eps));
		return g * factor;
	}

	static Quaternion damp_adjustment_exact_quat(Quaternion g, float halflife, float dt, float eps = 1e-8) {
		float factor = 1.0 - fast_negexp((CritDampSpring::Ln2 * dt) / (halflife + eps));
		return Quaternion().slerp(g, factor);
	}
	static Variant damper_exponential(Variant variable, Variant goal, float damping, float dt);

	static inline float fast_negexp(float x);

	static inline Variant damper_exact(Variant variable, Variant goal, float halflife, float dt, float eps = 1e-5f);

	static inline float halflife_to_damping(float halflife, float eps = 1e-5f);

	static inline float damping_to_halflife(float damping, float eps = 1e-5f);

	static inline float frequency_to_stiffness(float frequency);

	static inline float stiffness_to_frequency(float stiffness);

	static inline float critical_halflife(float frequency);

	static inline float critical_frequency(float halflife);

	static inline float damping_ratio_to_stiffness(float ratio, float damping);

	static inline float damping_ratio_to_damping(float ratio, float stiffness);

	static inline float maximum_spring_velocity_to_halflife(float x, float x_goal, float v_max);

	static void _spring_damper_exact(
			float &x,
			float &v,
			float x_goal,
			float v_goal,
			float damping_ratio,
			float halflife,
			float dt,
			float eps = 1e-5f);

	static void _critical_spring_damper_exact(
			float &x,
			float &v,
			float x_goal,
			float v_goal,
			float halflife,
			float dt);

	static inline PackedFloat32Array critical_spring_damper_exact(float x, float v, float x_goal, float v_goal, float halflife, float dt);

	static void _simple_spring_damper_exact(float &x, float &v, float x_goal, float halflife, float dt);
	static inline PackedFloat32Array simple_spring_damper_exact(float x, float v, float x_goal, float halflife, float dt);

	static inline void _decay_spring_damper_exact(float &x, float &v, float halflife, float dt);
	static inline PackedFloat32Array decay_spring_damper_exact(float x, float v, float halflife, float dt);

	//	Reach the x_goal at timed t_goal in the future
	//	Apprehension parameter controls how far into the future we try to track the linear interpolation
	static void _timed_spring_damper_exact(
			float &x, float &v,
			const float xi,
			const float &x_goal, const float &t_goal,
			const float &halflife, const float &dt,
			float apprehension = 2.0f);

	static inline PackedFloat32Array timed_spring_damper_exact(float x, float v,
			const float xi,
			const float x_goal, const float t_goal,
			const float halflife, const float dt,
			const float apprehension = 2.0f);

	static Dictionary character_update(
			Vector3 pos,
			Vector3 vel,
			Vector3 acc,
			Quaternion quaternion,
			Vector3 angular_velocity,
			Vector3 v_goal,
			Quaternion q_goal,
			float halflife_vel,
			float halflife_rot,
			float dt);

	static Dictionary character_predict(
			Vector3 x, Vector3 v, Vector3 a,
			Quaternion q, Vector3 angular_v,
			Vector3 v_goal, Quaternion q_goal,
			float halflife_v, float halflife_q,
			const PackedFloat32Array dts);

protected:
	static void _bind_methods();
};

#endif // CRIT_SPRING_DAMPER_H