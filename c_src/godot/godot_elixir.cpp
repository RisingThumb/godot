#include "godot_elixir.h"

#include "core/string/string_builder.h"
#include "main/main.h"
#include "modules/gdscript/gdscript.h"
#include "scene/main/node.h"

#include <locale.h>
#include <stdlib.h>
#include <unistd.h>

static OS_LinuxBSD os;

UNIFEX_TERM init(UnifexEnv *env, MyState *state, char **in_strings, unsigned int list_length) {
	int err = OK;
	if (state) {
		return init_result_fail(env, state, "Godot is already initialized.");
	}
	state = unifex_alloc_state(env);
	err = Main::setup(in_strings[0], list_length - 1, &in_strings[1]);
	if (err != OK) {		
		return init_result_fail(env, state, "Godot can't be setup.");
	}
	if (!Main::start()) {
		return init_result_fail(env, state, "Godot can't start.");
	}
	return init_result_ok(env, state, err);
}

UNIFEX_TERM call(UnifexEnv *env, MyState *state, char *method) {
	if (!state) {
		return init_result_fail(env, state, "Godot is not initialized.");
	}
	if (!os.get_main_loop()->get_script_instance()) {
		return init_result_fail(env, state, "Godot does not have a script instance.");
	}
	Callable::CallError call_error;
	Variant res = os.get_main_loop()->call(method, nullptr, 0, call_error);
	switch (res.get_type()) {
		case Variant::NIL: {
			return init_result_fail(env, state, "Call is invalid.");
		} break;
		case Variant::BOOL: {
			int res_int = res;
			return call_result_ok_bool(env, state, res_int);
		} break;
		case Variant::INT: {
			int res_int = res;
			return call_result_ok_int(env, state, res_int);
		} break;
		case Variant::FLOAT: {
			const CharString res_char_string = String(res).utf8();
			const char *res_string = res_char_string.get_data();
			return call_result_ok_string(env, state, res_string);
			// return call_result_ok_float(env, state, res.get_type(), res_float);
		} break;
		// Variant::DOUBLE
		case Variant::STRING: {
			const CharString res_char_string = String(res).utf8();
			const char *res_string = res_char_string.get_data();
			return call_result_ok_string(env, state, res_string);
		} break;
		case Variant::ARRAY: {
			const Array array = res;
			StringBuilder builder;
			for (int32_t array_i = 0; array_i < array.size(); array_i++) {
				builder.append(array[array_i]);
			}
			const CharString res_string = String(builder.as_string()).utf8();
			return call_result_ok_string(env, state, res_string);
		} break;
		default: {
			const CharString res_char_string = (String("Unsupported result: ") + String(res)).utf8();
			const char *res_string = res_char_string.get_data();
			return init_result_fail(env, state, res_string);
		} break;
	}
	return init_result_fail(env, state, "Call is invalid.");
}

void handle_destroy_state(UnifexEnv *env, MyState *state) {
	UNIFEX_UNUSED(env);
	UNIFEX_UNUSED(state);
	os.get_main_loop()->finalize();
	Main::cleanup();
}
