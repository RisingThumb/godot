/*************************************************************************/
/*  resonance_audio_wrapper.h                                            */
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

#ifndef RESONANCE_AUDIO_WRAPPER_H
#define RESONANCE_AUDIO_WRAPPER_H
#include "core/math/audio_frame.h"
#include "core/object/class_db.h"
#include "core/object/object.h"

#ifdef RESONANCEAUDIO_ENABLED
#include "thirdparty/resonanceaudio/resonance_audio/api/resonance_audio_api.h"

struct AudioSourceId {
	vraudio::ResonanceAudioApi::SourceId id = -1;
};

class ResonanceAudioWrapper : public Object {
	GDCLASS(ResonanceAudioWrapper, Object);
	static ResonanceAudioWrapper *singleton;
	static vraudio::ResonanceAudioApi *resonance_api;

public:
	ResonanceAudioWrapper();
	~ResonanceAudioWrapper() {
		if (resonance_api) {
			delete resonance_api;
			resonance_api = nullptr;
		}
	}
	static ResonanceAudioWrapper *get_singleton();

	AudioSourceId register_audio_source();
	AudioSourceId register_stero_audio_source();
	void unregister_audio_source(AudioSourceId audio_source);
	void set_source_transform(AudioSourceId source, Transform3D source_transform);
	void set_linear_source_volume(AudioSourceId audio_source, real_t volume);
	void set_head_transform(Transform3D head_transform);

	void push_source_buffer(AudioSourceId source, int num_frames, AudioFrame *frames);
	bool pull_listener_buffer(int num_frames, AudioFrame *frames);

	void set_source_attenuation(AudioSourceId source, float attenuation_linear);
};
#endif
#endif
