/*************************************************************************/
/*  os_osx.h                                                             */
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

#ifndef OS_OSX_H
#define OS_OSX_H

#include "core/input/input.h"
#include "crash_handler_osx.h"
#include "drivers/coreaudio/audio_driver_coreaudio.h"
#include "drivers/coremidi/midi_driver_coremidi.h"
#include "drivers/unix/os_unix.h"
#include "joypad_osx.h"
#include "servers/audio_server.h"

class OS_OSX : public OS_Unix {
	bool force_quit = false;

	JoypadOSX *joypad_osx = nullptr;

#ifdef COREAUDIO_ENABLED
	AudioDriverCoreAudio audio_driver;
#endif
#ifdef COREMIDI_ENABLED
	MIDIDriverCoreMidi midi_driver;
#endif

	CrashHandler crash_handler;

	CFRunLoopObserverRef pre_wait_observer;

	MainLoop *main_loop = nullptr;

	List<String> launch_service_args;

	static _FORCE_INLINE_ String get_framework_executable(const String &p_path);
	static void pre_wait_observer_cb(CFRunLoopObserverRef p_observer, CFRunLoopActivity p_activiy, void *p_context);

protected:
	virtual void initialize_core() override;
	virtual void initialize() override;
	virtual void finalize() override;

	virtual void initialize_joypads() override;

	virtual void set_main_loop(MainLoop *p_main_loop) override;
	virtual void delete_main_loop() override;

public:
	virtual void set_cmdline_platform_args(const List<String> &p_args);
	virtual List<String> get_cmdline_platform_args() const override;

	virtual String get_name() const override;

	virtual void alert(const String &p_alert, const String &p_title = "ALERT!") override;

	virtual Error open_dynamic_library(const String p_path, void *&p_library_handle, bool p_also_set_library_path = false, String *r_resolved_path = nullptr) override;

	virtual MainLoop *get_main_loop() const override;

	virtual String get_config_path() const override;
	virtual String get_data_path() const override;
	virtual String get_cache_path() const override;
	virtual String get_bundle_resource_dir() const override;
	virtual String get_bundle_icon_path() const override;
	virtual String get_godot_dir_name() const override;

	virtual String get_system_dir(SystemDir p_dir, bool p_shared_storage = true) const override;

	virtual Error shell_open(String p_uri) override;
	virtual Error register_protocol(String p_protocol) override;
	virtual Error unregister_protocol(String p_protocol) override;

	virtual String get_locale() const override;

	virtual String get_executable_path() const override;
	virtual Error create_process(const String &p_path, const List<String> &p_arguments, ProcessID *r_child_id = nullptr, bool p_open_console = false) override;
	virtual Error create_instance(const List<String> &p_arguments, ProcessID *r_child_id = nullptr) override;

	virtual String get_unique_id() const override;
	virtual String get_processor_name() const override;

	virtual bool _check_internal_feature_support(const String &p_feature) override;

	virtual void disable_crash_handler() override;
	virtual bool is_disable_crash_handler() const override;

	virtual Error move_to_trash(const String &p_path) override;

	void run();

	OS_OSX();
	~OS_OSX();
};

#endif
