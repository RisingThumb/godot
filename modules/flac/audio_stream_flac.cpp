/**************************************************************************/
/*  audio_stream_flac.cpp                                                 */
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

#include <cstddef>
#define DR_FLAC_IMPLEMENTATION
#define DR_FLAC_NO_STDIO
#define DR_FLAC_NO_OGG
#include "audio_stream_flac.h"

#include "thirdparty/dr_flac/dr_flac.h"

#include "core/io/file_access.h"

int AudioStreamPlaybackFLAC::_mix_internal(AudioFrame *p_buffer, int p_frames) {
	if (!active) {
		return 0;
	}

	int todo = p_frames;
	int start_buffer = 0;

	int frames_mixed_this_step = p_frames;

	int beat_length_frames = -1;
	bool beat_loop = flac_stream->has_loop() && flac_stream->get_bpm() > 0 && flac_stream->get_beat_count() > 0;
	if (beat_loop) {
		beat_length_frames = flac_stream->get_beat_count() * flac_stream->sample_rate * 60 / flac_stream->get_bpm();
	}

	while (todo && active) {
		float *buffer = (float *)p_buffer;
		if (start_buffer > 0) {
			buffer = (buffer + start_buffer * 2);
		}
		int mixed = drflac_read_pcm_frames_f32(pFlac, todo, buffer);

		for (int i = 0; i < mixed; i++) {
			if (loop_fade_remaining < FADE_SIZE) {
				p_buffer[p_frames - todo] += loop_fade[loop_fade_remaining] * (float(FADE_SIZE - loop_fade_remaining) / float(FADE_SIZE));
				loop_fade_remaining++;
			}

			--todo;
			++frames_mixed;

			if (beat_loop && (int)frames_mixed >= beat_length_frames) {
				for (int i = 0; i < FADE_SIZE; i++) {
					p_buffer[i] = AudioFrame(buffer[i * flac_stream->channels], buffer[i * flac_stream->channels + flac_stream->channels - 1]);
				}
				loop_fade_remaining = 0;
				seek(flac_stream->loop_offset);
				loops++;
			}
		}

		if (todo) {
			//end of file!
			if (flac_stream->loop) {
				//loop
				seek(flac_stream->loop_offset);
				loops++;
				// we still have buffer to fill, start from this element in the next iteration.
				start_buffer = p_frames - todo;
			} else {
				frames_mixed_this_step = p_frames - todo;

				for (int i = p_frames - todo; i < p_frames; i++) {
					p_buffer[i] = AudioFrame(0, 0);
				}
				active = false;
				todo = 0;
			}
		}
	}

	return frames_mixed_this_step;
}

float AudioStreamPlaybackFLAC::get_stream_sampling_rate() {
	return flac_stream->sample_rate;
}

void AudioStreamPlaybackFLAC::start(double p_from_pos) {
	active = true;
	seek(p_from_pos);
	loops = 0;
	begin_resample();
}

void AudioStreamPlaybackFLAC::stop() {
	active = false;
}

bool AudioStreamPlaybackFLAC::is_playing() const {
	return active;
}

int AudioStreamPlaybackFLAC::get_loop_count() const {
	return loops;
}

double AudioStreamPlaybackFLAC::get_playback_position() const {
	return double(frames_mixed) / double(flac_stream->sample_rate);
}

void AudioStreamPlaybackFLAC::seek(double p_time) {
	if (!active) {
		return;
	}

	if (p_time >= flac_stream->get_length()) {
		p_time = 0;
	}
	frames_mixed = flac_stream->sample_rate * p_time;
	drflac_seek_to_pcm_frame(pFlac, frames_mixed);
}

void AudioStreamPlaybackFLAC::tag_used_streams() {
	flac_stream->tag_used(get_playback_position());
}

AudioStreamPlaybackFLAC::~AudioStreamPlaybackFLAC() {
	if (pFlac) {
		drflac_close(pFlac);
		pFlac = nullptr;
	}
}

Ref<AudioStreamPlayback> AudioStreamFLAC::instantiate_playback() {
	Ref<AudioStreamPlaybackFLAC> flacs;

	ERR_FAIL_COND_V_MSG(data.is_empty(), flacs,
			"This AudioStreamFLAC does not have an audio file assigned "
			"to it. AudioStreamFLAC should not be created from the "
			"inspector or with `.new()`. Instead, load an audio file.");

	flacs.instantiate();
	flacs->flac_stream = Ref<AudioStreamFLAC>(this);

	flacs->pFlac = drflac_open_memory(data.ptr(), data_len, nullptr);

	flacs->frames_mixed = 0;
	flacs->active = false;
	flacs->loops = 0;

	if (!flacs->pFlac) {
		ERR_FAIL_COND_V(!flacs->pFlac, Ref<AudioStreamPlaybackFLAC>());
	}

	return flacs;
}

String AudioStreamFLAC::get_stream_name() const {
	return "";
}

void AudioStreamFLAC::clear_data() {
	data.clear();
}

void AudioStreamFLAC::set_data(const Vector<uint8_t> &p_data) {
	int src_data_len = p_data.size();

	const uint8_t *src_datar = p_data.ptr();

	drflac *pflac = drflac_open_memory(src_datar, src_data_len, nullptr);
	ERR_FAIL_COND(pflac == nullptr);

	if (pflac->channels != 2) {
		ERR_FAIL_MSG("Number of channels must be exactly 2.");
	}

	channels = pflac->channels;
	sample_rate = pflac->sampleRate;
	length = float(pflac->totalPCMFrameCount) / float(sample_rate);

	clear_data();

	data.resize(src_data_len);
	memcpy(data.ptrw(), src_datar, src_data_len);
	data_len = src_data_len;
}

Vector<uint8_t> AudioStreamFLAC::get_data() const {
	return data;
}

void AudioStreamFLAC::set_loop(bool p_enable) {
	loop = p_enable;
}

bool AudioStreamFLAC::has_loop() const {
	return loop;
}

void AudioStreamFLAC::set_loop_offset(double p_seconds) {
	loop_offset = p_seconds;
}

double AudioStreamFLAC::get_loop_offset() const {
	return loop_offset;
}

double AudioStreamFLAC::get_length() const {
	return length;
}

bool AudioStreamFLAC::is_monophonic() const {
	return false;
}

void AudioStreamFLAC::set_bpm(double p_bpm) {
	ERR_FAIL_COND(p_bpm < 0);
	bpm = p_bpm;
	emit_changed();
}

double AudioStreamFLAC::get_bpm() const {
	return bpm;
}

void AudioStreamFLAC::set_beat_count(int p_beat_count) {
	ERR_FAIL_COND(p_beat_count < 0);
	beat_count = p_beat_count;
	emit_changed();
}

int AudioStreamFLAC::get_beat_count() const {
	return beat_count;
}

void AudioStreamFLAC::set_bar_beats(int p_bar_beats) {
	ERR_FAIL_COND(p_bar_beats < 0);
	bar_beats = p_bar_beats;
	emit_changed();
}

int AudioStreamFLAC::get_bar_beats() const {
	return bar_beats;
}

void AudioStreamFLAC::_bind_methods() {
	ClassDB::bind_method(D_METHOD("set_data", "data"), &AudioStreamFLAC::set_data);
	ClassDB::bind_method(D_METHOD("get_data"), &AudioStreamFLAC::get_data);

	ClassDB::bind_method(D_METHOD("set_loop", "enable"), &AudioStreamFLAC::set_loop);
	ClassDB::bind_method(D_METHOD("has_loop"), &AudioStreamFLAC::has_loop);

	ClassDB::bind_method(D_METHOD("set_loop_offset", "seconds"), &AudioStreamFLAC::set_loop_offset);
	ClassDB::bind_method(D_METHOD("get_loop_offset"), &AudioStreamFLAC::get_loop_offset);

	ClassDB::bind_method(D_METHOD("set_bpm", "bpm"), &AudioStreamFLAC::set_bpm);
	ClassDB::bind_method(D_METHOD("get_bpm"), &AudioStreamFLAC::get_bpm);

	ClassDB::bind_method(D_METHOD("set_beat_count", "count"), &AudioStreamFLAC::set_beat_count);
	ClassDB::bind_method(D_METHOD("get_beat_count"), &AudioStreamFLAC::get_beat_count);

	ClassDB::bind_method(D_METHOD("set_bar_beats", "count"), &AudioStreamFLAC::set_bar_beats);
	ClassDB::bind_method(D_METHOD("get_bar_beats"), &AudioStreamFLAC::get_bar_beats);

	ADD_PROPERTY(PropertyInfo(Variant::PACKED_BYTE_ARRAY, "data", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_NO_EDITOR), "set_data", "get_data");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "bpm", PROPERTY_HINT_RANGE, "0,400,0.01,or_greater"), "set_bpm", "get_bpm");
	ADD_PROPERTY(PropertyInfo(Variant::INT, "beat_count", PROPERTY_HINT_RANGE, "0,512,1,or_greater"), "set_beat_count", "get_beat_count");
	ADD_PROPERTY(PropertyInfo(Variant::INT, "bar_beats", PROPERTY_HINT_RANGE, "2,32,1,or_greater"), "set_bar_beats", "get_bar_beats");
	ADD_PROPERTY(PropertyInfo(Variant::BOOL, "loop", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_NO_EDITOR), "set_loop", "has_loop");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "loop_offset", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_NO_EDITOR), "set_loop_offset", "get_loop_offset");
}

AudioStreamFLAC::AudioStreamFLAC() {
}

AudioStreamFLAC::~AudioStreamFLAC() {
	clear_data();
}
