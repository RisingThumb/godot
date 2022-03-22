/*************************************************************************/
/*  resonance_audio_wrapper.cpp                                          */
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

#ifdef RESONANCEAUDIO_ENABLED
#include "resonance_audio_wrapper.h"
#include "core/error/error_macros.h"
#include "servers/audio_server.h"

ResonanceAudioWrapper *ResonanceAudioWrapper::singleton = nullptr;
vraudio::ResonanceAudioApi *ResonanceAudioWrapper::resonance_api = nullptr;

ResonanceAudioWrapper::ResonanceAudioWrapper() {
	singleton = this;
	if (!AudioServer::get_singleton()) {
		return;
	}
	if (resonance_api) {
		return;
	}
	resonance_api = vraudio::CreateResonanceAudioApi(
			/* num_channels= */ 2, AudioServer::get_singleton()->thread_get_mix_buffer_size(), AudioServer::get_singleton()->get_mix_rate());
}

ResonanceAudioWrapper *ResonanceAudioWrapper::get_singleton() {
	return singleton;
}

AudioSourceId ResonanceAudioWrapper::register_audio_source() {
	AudioSourceId new_source;
	ERR_FAIL_NULL_V(resonance_api, new_source);
	new_source.id = resonance_api->CreateSoundObjectSource(vraudio::RenderingMode::kBinauralHighQuality);
	resonance_api->SetSourceDistanceModel(
			new_source.id,
			vraudio::DistanceRolloffModel::kNone,
			/* min_distance= */ 0,
			/* max_distance= */ 0);

	return new_source;
}

void ResonanceAudioWrapper::unregister_audio_source(AudioSourceId audio_source) {
	ERR_FAIL_NULL(resonance_api);
	resonance_api->DestroySource(audio_source.id);
}

void ResonanceAudioWrapper::set_source_transform(AudioSourceId source, Transform3D source_transform) {
	ERR_FAIL_NULL(resonance_api);
	Quaternion source_rotation = Quaternion(source_transform.basis);
	resonance_api->SetSourcePosition(source.id, source_transform.origin.x, source_transform.origin.y, source_transform.origin.z);
	resonance_api->SetSourceRotation(source.id, source_rotation.x, source_rotation.y, source_rotation.z, source_rotation.w);
}

void ResonanceAudioWrapper::set_head_transform(Transform3D head_transform) {
	ERR_FAIL_NULL(resonance_api);
	Quaternion head_rotation = Quaternion(head_transform.basis);
	resonance_api->SetHeadPosition(head_transform.origin.x, head_transform.origin.y, head_transform.origin.z);
	resonance_api->SetHeadRotation(head_rotation.x, head_rotation.y, head_rotation.z, head_rotation.w);
}

void ResonanceAudioWrapper::push_source_buffer(AudioSourceId source, int num_frames, AudioFrame *frames) {
	ERR_FAIL_NULL(resonance_api);
	// Frames are just interleaved floats.
	resonance_api->SetInterleavedBuffer(source.id, (const float *)frames, /* num_channels= */ 2, num_frames);
}

bool ResonanceAudioWrapper::pull_listener_buffer(int num_frames, AudioFrame *frames) {
	ERR_FAIL_NULL_V(resonance_api, false);
	// Frames are just interleaved floats.
	bool success = resonance_api->FillInterleavedOutputBuffer(/* num_channels= */ 2, num_frames, (float *)frames);
	if (!success) {
		// Zero out the array because Resonance Audio fills buffers with garbage on error under some circumstances.
		memset(frames, 0, num_frames * sizeof(AudioFrame));
	}
	return success;
}

void ResonanceAudioWrapper::set_source_attenuation(AudioSourceId source, float attenuation_linear) {
	ERR_FAIL_NULL(resonance_api);
	resonance_api->SetSourceDistanceAttenuation(source.id, attenuation_linear);
}

AudioSourceId ResonanceAudioWrapper::ResonanceAudioWrapper::register_stero_audio_source() {
	AudioSourceId new_source;
	ERR_FAIL_NULL_V(resonance_api, new_source);
	new_source.id = resonance_api->CreateStereoSource(2);
	return new_source;
}

void ResonanceAudioWrapper::set_linear_source_volume(AudioSourceId audio_source, real_t volume) {
	ERR_FAIL_NULL(resonance_api);
	resonance_api->SetSourceVolume(audio_source.id, volume);
}
#endif
