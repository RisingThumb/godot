/*************************************************************************/
/*  xr_positional_tracker.h                                              */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                      https://godotengine.org                          */
/*************************************************************************/
/* Copyright (c) 2007-2022 Juan Linietsky, Ariel Manzur.                 */
/* Copyright (c) 2014-2022 Godot Engine contributors (cf. AUTHORS.md).   */
/*                                                                       */
/* Permission is hereby granted, free of charge, to any person obtaining */
/* a copy of this software and associated documentation files (the       */
/* "Software"), to deal in the Software without restriction, including   */
/* without limitation the rights to use, copy, modify, merge, publish,   */
/* distribute, sublicense, and/or sell copies of the Software, and to    */
/* permit persons to whom the Software is furnished to do so, subject to */
/* the following conditions:                                             */
/*                                                                       */
/* The above copyright notice and this permission notice shall be        */
/* included in all copies or substantial portions of the Software.       */
/*                                                                       */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,       */
/* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF    */
/* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.*/
/* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY  */
/* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,  */
/* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE     */
/* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                */
/*************************************************************************/

#ifndef XR_POSITIONAL_TRACKER_H
#define XR_POSITIONAL_TRACKER_H

#include "core/os/thread_safe.h"
#include "scene/resources/mesh.h"
#include "servers/xr/xr_pose.h"
#include "servers/xr_server.h"

/**
	The positional tracker object as an object that represents the position and orientation of a tracked object like a controller or headset.
	An AR/VR Interface will registered the trackers it manages with our AR/VR server and update its position and orientation.
	This is where potentially additional AR/VR interfaces may be active as there are AR/VR SDKs that solely deal with positional tracking.
*/

class XRPositionalTracker : public RefCounted {
	GDCLASS(XRPositionalTracker, RefCounted);
	_THREAD_SAFE_CLASS_

public:
	enum TrackerHand {
		TRACKER_HAND_UNKNOWN, /* unknown or not applicable */
		TRACKER_HAND_LEFT, /* controller is the left hand controller */
		TRACKER_HAND_RIGHT /* controller is the right hand controller */
	};

private:
	XRServer::TrackerType type; // type of tracker
	StringName name; // (unique) name of the tracker
	String description; // description of the tracker, this is interface dependent, for OpenXR this will be the interaction profile bound for to the tracker
	TrackerHand hand; // if known, the hand this tracker is held in

	Map<StringName, Ref<XRPose>> poses;
	Map<StringName, Variant> inputs;

protected:
	static void _bind_methods();

public:
	void set_tracker_type(XRServer::TrackerType p_type);
	XRServer::TrackerType get_tracker_type() const;
	void set_tracker_name(const StringName &p_name);
	StringName get_tracker_name() const;
	void set_tracker_desc(const String &p_desc);
	String get_tracker_desc() const;
	XRPositionalTracker::TrackerHand get_tracker_hand() const;
	void set_tracker_hand(const XRPositionalTracker::TrackerHand p_hand);

	bool has_pose(const StringName &p_action_name) const;
	Ref<XRPose> get_pose(const StringName &p_action_name) const;
	void invalidate_pose(const StringName &p_action_name);
	void set_pose(const StringName &p_action_name, const Transform3D &p_transform, const Vector3 &p_linear_velocity, const Vector3 &p_angular_velocity, const XRPose::TrackingConfidence p_tracking_confidence = XRPose::XR_TRACKING_CONFIDENCE_HIGH);

	Variant get_input(const StringName &p_action_name) const;
	void set_input(const StringName &p_action_name, const Variant &p_value);

	XRPositionalTracker();
	~XRPositionalTracker() {}
};

VARIANT_ENUM_CAST(XRPositionalTracker::TrackerHand);

#endif // XR_POSITIONAL_TRACKER_H
