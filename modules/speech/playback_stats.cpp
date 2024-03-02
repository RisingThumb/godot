#include "playback_stats.h"
#include "speech_processor.h"

Dictionary PlaybackStats::get_playback_stats() {
	double playback_pushed_frames = playback_pushed_calls * (buffer_frame_count * 1.0);
	double playback_discarded_frames = playback_discarded_calls * (buffer_frame_count * 1.0);
	Dictionary dict;
	dict["playback_ring_limit_s"] = playback_ring_buffer_length / double(SpeechProcessor::SPEECH_SETTING_VOICE_PACKET_SAMPLE_RATE);
	dict["playback_ring_current_size_s"] = playback_ring_current_size / double(SpeechProcessor::SPEECH_SETTING_VOICE_PACKET_SAMPLE_RATE);
	dict["playback_ring_max_size_s"] = playback_ring_max_size / double(SpeechProcessor::SPEECH_SETTING_VOICE_PACKET_SAMPLE_RATE);
	dict["playback_ring_mean_size_s"] = 0;
	if (playback_push_buffer_calls > 0) {
		dict["playback_ring_mean_size_s"] = playback_ring_size_sum / playback_push_buffer_calls / double(SpeechProcessor::SPEECH_SETTING_VOICE_PACKET_SAMPLE_RATE);
	} else {
		dict["playback_ring_mean_size_s"] = 0;
	}
	dict["jitter_buffer_current_size_s"] = float(jitter_buffer_current_size) * SpeechProcessor::SPEECH_SETTING_PACKET_DELTA_TIME;
	dict["jitter_buffer_max_size_s"] = float(jitter_buffer_max_size) * SpeechProcessor::SPEECH_SETTING_PACKET_DELTA_TIME;
	dict["jitter_buffer_mean_size_s"] = 0;
	if (jitter_buffer_calls > 0) {
		dict["jitter_buffer_mean_size_s"] = float(jitter_buffer_size_sum) / jitter_buffer_calls * SpeechProcessor::SPEECH_SETTING_PACKET_DELTA_TIME;
	}
	dict["jitter_buffer_calls"] = jitter_buffer_calls;
	dict["playback_position_s"] = playback_position;
	dict["playback_get_percent"] = 0;
	dict["playback_discard_percent"] = 0;
	if (playback_pushed_frames > 0) {
		dict["playback_get_percent"] = 100.0 * playback_get_frames / playback_pushed_frames;
		dict["playback_discard_percent"] = 100.0 * playback_discarded_frames / playback_pushed_frames;
	}
	dict["playback_get_s"] = playback_get_frames / double(SpeechProcessor::SPEECH_SETTING_VOICE_PACKET_SAMPLE_RATE);
	dict["playback_pushed_s"] = playback_pushed_frames / double(SpeechProcessor::SPEECH_SETTING_VOICE_PACKET_SAMPLE_RATE);
	dict["playback_discarded_s"] = playback_discarded_frames / double(SpeechProcessor::SPEECH_SETTING_VOICE_PACKET_SAMPLE_RATE);
	dict["playback_push_buffer_calls"] = floor(playback_push_buffer_calls);
	dict["playback_blank_s"] = playback_blank_push_calls * SpeechProcessor::SPEECH_SETTING_PACKET_DELTA_TIME;
	dict["playback_blank_percent"] = 0;
	if (playback_push_buffer_calls > 0) {
		dict["playback_blank_percent"] = 100.0 * playback_blank_push_calls / playback_push_buffer_calls;
	}
	dict["playback_skips"] = floor(playback_skips);
	return dict;
}

void PlaybackStats::_bind_methods() {
	ClassDB::bind_method(D_METHOD("get_playback_stats"),
			&PlaybackStats::get_playback_stats);
}
