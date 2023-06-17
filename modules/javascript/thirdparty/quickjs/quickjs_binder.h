#ifndef QUICKJS_BINDING_HELPER_H
#define QUICKJS_BINDING_HELPER_H

#include "core/io/resource.h"
#include "core/os/memory.h"
#include "core/os/thread.h"
#include "core/variant/callable.h"

#include "../../javascript_binder.h"
#include "./quickjs/quickjs.h"
#include "quickjs_builtin_binder.h"

#ifdef QUICKJS_WITH_DEBUGGER
#include "quickjs_debugger.h"
#endif

#define JS_HIDDEN_SYMBOL(x) ("\xFF" x)
#define BINDING_DATA_FROM_JS(ctx, p_val) (JavaScriptGCHandler *)JS_GetOpaque((p_val), QuickJSBinder::get_origin_class_id((ctx)))
#define GET_JSVALUE(p_gc_handler) JS_MKPTR(JS_TAG_OBJECT, (p_gc_handler).javascript_object)
#define NO_MODULE_EXPORT_SUPPORT 0
#define MODULE_HAS_REFCOUNT 0 // The module doesn't seem to follow the reference count rule in quickjs.
#define MAX_ARGUMENT_COUNT 50
#define PROP_NAME_CONSOLE_LOG_OBJECT_TO_JSON "LOG_OBJECT_TO_JSON"
#define ENDL "\r\n"

class QuickJSBinder : public JavaScriptBinder {
	friend class QuickJSBuiltinBinder;
	friend class QuickJSWorker;
	QuickJSBuiltinBinder builtin_binder;

protected:
	static SafeNumeric<uint32_t> global_context_id;
	static SafeNumeric<uint64_t> global_transfer_id;
	JSRuntime *runtime;
	JSContext *ctx;
	JSMallocFunctions godot_allocator;
	uint32_t context_id;
	Thread::ID thread_id;
#ifdef QUICKJS_WITH_DEBUGGER
	Ref<QuickJSDebugger> debugger;
#endif

public:
	struct PtrHasher {
		static _FORCE_INLINE_ uint32_t hash(const void *p_ptr) {
			union {
				const void *p;
				uint64_t i;
			} u;
			u.p = p_ptr;
			return HashMapHasherDefault::hash(u.i);
		}
	};

	enum JavaScriptModuleFlags {
		MODULE_FLAG_SCRIPT = 1,
		MODULE_FLAG_NATIVE = 2,
		MODULE_FLAG_EVALUATED = 1 << 3,
		MODULE_FLAG_RESOURCE = 1 << 4,
	};

	struct ModuleCache {
		int flags = 0;
		JSModuleDef *module = NULL;
		uint32_t hash = 0;
		JSValue res_value;
		Ref<Resource> res;
	};

	struct CommonJSModule {
		int flags;
		JSValue exports;
		String md5;
		Ref<Resource> res;
	};

	enum {
		PROP_DEF_DEFAULT = JS_PROP_ENUMERABLE | JS_PROP_CONFIGURABLE,
	};

	enum {
		__JS_ATOM_NULL = JS_ATOM_NULL,
#if !defined(EMSCRIPTEN)
#define CONFIG_ATOMICS
#endif
#define DEF(name, str) JS_ATOM_##name,
#include "quickjs/quickjs-atom.h"
#undef DEF
#ifdef CONFIG_ATOMICS
#undef CONFIG_ATOMICS
#endif
		JS_ATOM_END,
	};

	enum ConsoleLogType {
		CONSOLE_LOG,
		CONSOLE_WARN,
		CONSOLE_ERROR,
		CONSOLE_TRACE,
	};

protected:
	JSAtom js_key_godot_classid;
	JSAtom js_key_godot_classname;
	JSAtom js_key_godot_tooled;
	JSAtom js_key_godot_icon_path;
	JSAtom js_key_godot_exports;
	JSAtom js_key_godot_signals;

	JSValue global_object;
	JSValue godot_object;
	JSValue console_object;
	JSValue empty_function;
	JSValue js_operators;
	JSValue js_operators_create;
	Vector<JSValue> godot_singletons;

	_FORCE_INLINE_ static void *js_binder_malloc(JSMallocState *s, size_t size) { return memalloc(size); }
	_FORCE_INLINE_ static void js_binder_free(JSMallocState *s, void *ptr) {
		if (ptr)
			memfree(ptr);
	}
	_FORCE_INLINE_ static void *js_binder_realloc(JSMallocState *s, void *ptr, size_t size) { return memrealloc(ptr, size); }

	static String resolve_module_file(const String &file);
	static JSModuleDef *js_module_loader(JSContext *ctx, const char *module_name, void *opaque);
	static JSModuleDef *js_make_module(JSContext *ctx, const String &p_id, const JSValueConst &p_value);
	ModuleCache *js_compile_and_cache_module(JSContext *ctx, const String &p_code, const String &p_filename, JavaScriptError *r_error);
	ModuleCache *js_compile_and_cache_module(JSContext *ctx, const Vector<uint8_t> &p_bytecode, const String &p_filename, JavaScriptError *r_error);
	ModuleCache js_compile_module(JSContext *ctx, const String &p_code, const String &p_filename, JavaScriptError *r_error);
	static Error js_evalute_module(JSContext *ctx, ModuleCache *p_module, JavaScriptError *r_error);
	static int resource_module_initializer(JSContext *ctx, JSModuleDef *m);

	struct ClassBindData {
		JSClassID class_id;
		CharString class_name;
		JSValue prototype;
		JSValue constructor;
		JSClassDef jsclass;
		const ClassDB::ClassInfo *gdclass;
		const ClassBindData *base_class;
	};
	ClassBindData godot_origin_class;
	const ClassBindData *godot_object_class;
	const ClassBindData *godot_reference_class;
	HashMap<JSClassID, ClassBindData> class_bindings;
	HashMap<StringName, const ClassBindData *> classname_bindings;
	HashMap<String, ModuleCache> module_cache;
	HashMap<String, CommonJSModule> commonjs_module_cache;
	ClassBindData worker_class_data;
	List<JavaScriptGCHandler *> workers;
	Vector<MethodBind *> godot_methods;
	int internal_godot_method_id;
	Vector<const ClassDB::PropertySetGet *> godot_object_indexed_properties;
	int internal_godot_indexed_property_id;
#if NO_MODULE_EXPORT_SUPPORT
	String parsing_script_file;
#endif

	JSClassID register_class(const ClassDB::ClassInfo *p_cls);
	void add_godot_origin();
	void add_godot_classes();
	void add_godot_globals();
	void add_global_console();
	void add_global_properties();
	void add_global_worker();

	static JSValue object_constructor(JSContext *ctx, JSValueConst this_val, int argc, JSValueConst *argv, int class_id);
	static void initialize_properties(JSContext *ctx, const JavaScriptClassInfo *p_class, JSValue p_object);

	static void object_finalizer(JavaScriptGCHandler *p_bind);
	static void origin_finalizer(JSRuntime *rt, JSValue val);

	static JSValue object_free(JSContext *ctx, JSValueConst this_val, int argc, JSValueConst *argv);
	static JSValue object_method(JSContext *ctx, JSValueConst this_val, int argc, JSValueConst *argv, int method_id);
	static JSValue object_indexed_property(JSContext *ctx, JSValueConst this_val, int argc, JSValueConst *argv, int property_id);
	static JSValue godot_to_string(JSContext *ctx, JSValueConst this_val, int argc, JSValueConst *argv);
	static JSValue godot_get_type(JSContext *ctx, JSValueConst this_val, int argc, JSValueConst *argv);
	static JSValue godot_load(JSContext *ctx, JSValueConst this_val, int argc, JSValueConst *argv);
	static JSValue godot_instance_from_id(JSContext *ctx, JSValueConst this_val, int argc, JSValueConst *argv);

	static JSValue godot_object_method_connect(JSContext *ctx, JSValueConst this_val, int argc, JSValueConst *argv);

	static void add_debug_binding_info(JSContext *ctx, JSValueConst p_obj, const JavaScriptGCHandler *p_bind);

	const JavaScriptClassInfo *register_javascript_class(const JSValueConst &p_constructor, const String &p_path);
	void free_javascript_class(const JavaScriptClassInfo &p_class);
	static JSValue godot_register_class(JSContext *ctx, JSValueConst this_val, int argc, JSValueConst *argv);
	static JSValue godot_register_signal(JSContext *ctx, JSValueConst this_val, int argc, JSValueConst *argv);
	static JSValue godot_register_property(JSContext *ctx, JSValueConst this_val, int argc, JSValueConst *argv);

	enum {
		SCRIPT_META_TOOLED,
		SCRIPT_META_ICON,
	};
	static JSValue godot_set_script_meta(JSContext *ctx, JSValueConst this_val, int argc, JSValueConst *argv, int magic);

	static JSValue console_functions(JSContext *ctx, JSValue this_val, int argc, JSValue *argv, int magic);
	static JSValue global_request_animation_frame(JSContext *ctx, JSValueConst this_val, int argc, JSValueConst *argv);
	static JSValue global_cancel_animation_frame(JSContext *ctx, JSValueConst this_val, int argc, JSValueConst *argv);

	static JSValue worker_constructor(JSContext *ctx, JSValueConst this_val, int argc, JSValueConst *argv);
	static void worker_finializer(JSRuntime *rt, JSValue val);
	static JSValue worker_post_message(JSContext *ctx, JSValueConst this_val, int argc, JSValueConst *argv);
	static JSValue worker_terminate(JSContext *ctx, JSValueConst this_val, int argc, JSValueConst *argv);
	static JSValue godot_abandon_value(JSContext *ctx, JSValueConst this_val, int argc, JSValueConst *argv);
	static JSValue godot_adopt_value(JSContext *ctx, JSValueConst this_val, int argc, JSValueConst *argv);

	_FORCE_INLINE_ static JSValue js_empty_func(JSContext *ctx, JSValueConst this_val, int argc, JSValueConst *argv) { return JS_UNDEFINED; }
	_FORCE_INLINE_ static JSValue js_empty_consturctor(JSContext *ctx, JSValueConst this_val, int argc, JSValueConst *argv) { return JS_NewObject(ctx); }
	static JSValue godot_builtin_function(JSContext *ctx, JSValueConst this_val, int argc, JSValueConst *argv, int magic);

	static int get_js_array_length(JSContext *ctx, JSValue p_val);
	static void get_own_property_names(JSContext *ctx, JSValue p_object, HashSet<String> *r_list);

	static JSAtom get_atom(JSContext *ctx, const StringName &p_key);
	static HashMap<uint64_t, Variant> transfer_deopot;
	static HashMap<String, const char *> class_remap;
#ifdef TOOLS_ENABLED
	Dictionary modified_api;
#endif
public:
	static Error define_operators(JSContext *ctx, JSValue p_prototype, JSValue *p_operators, int p_size);

public:
	static JSValue variant_to_var(JSContext *ctx, const Variant p_var);
	static Variant var_to_variant(JSContext *ctx, JSValue p_val);
	static bool validate_type(JSContext *ctx, Variant::Type p_type, JSValueConst &p_val);
	static void dump_exception(JSContext *ctx, const JSValueConst &p_exception, JavaScriptError *r_error);
	virtual String error_to_string(const JavaScriptError &p_error) override;
	virtual Error get_stacks(List<JavaScriptStackInfo> &r_stacks) override;
	virtual String get_backtrace_message(const List<JavaScriptStackInfo> &stacks) override;

	static Dictionary js_to_dictionary(JSContext *ctx, const JSValueConst &p_val, List<void *> &stack);

	_FORCE_INLINE_ static real_t js_to_number(JSContext *ctx, const JSValueConst &p_val) {
		double_t v = 0;
		JS_ToFloat64(ctx, &v, p_val);
		return real_t(v);
	}
	_FORCE_INLINE_ static String js_to_string(JSContext *ctx, const JSValueConst &p_val) {
		String ret;
		size_t len = 0;
		const char *utf8 = JS_ToCStringLen(ctx, &len, p_val);
		ret.parse_utf8(utf8, len);
		JS_FreeCString(ctx, utf8);
		return ret;
	}
	_FORCE_INLINE_ static bool js_to_bool(JSContext *ctx, const JSValueConst &p_val) {
		return JS_ToBool(ctx, p_val);
	}
	_FORCE_INLINE_ static int32_t js_to_int(JSContext *ctx, const JSValueConst &p_val) {
		int32_t i;
		JS_ToInt32(ctx, &i, p_val);
		return i;
	}
	_FORCE_INLINE_ static uint32_t js_to_uint(JSContext *ctx, const JSValueConst &p_val) {
		uint32_t u;
		JS_ToUint32(ctx, &u, p_val);
		return u;
	}
	_FORCE_INLINE_ static int64_t js_to_int64(JSContext *ctx, const JSValueConst &p_val) {
		int64_t i;
		JS_ToInt64(ctx, &i, p_val);
		return i;
	}
	_FORCE_INLINE_ static uint64_t js_to_uint64(JSContext *ctx, const JSValueConst &p_val) {
		uint64_t i;
		JS_ToIndex(ctx, &i, p_val);
		return i;
	}
	_FORCE_INLINE_ static JSValue to_js_number(JSContext *ctx, real_t p_val) {
		return JS_NewFloat64(ctx, double(p_val));
	}
	_FORCE_INLINE_ static JSValue to_js_string(JSContext *ctx, const String &text) {
		CharString utf8 = text.utf8();
		return JS_NewStringLen(ctx, utf8.get_data(), utf8.length());
	}
	_FORCE_INLINE_ static JSValue to_js_bool(JSContext *ctx, bool p_val) {
		return JS_NewBool(ctx, p_val);
	}

public:

	QuickJSBinder() {
		context_id = QuickJSBinder::global_context_id.increment();
		internal_godot_method_id = 0;
		internal_godot_indexed_property_id = 0;
		godot_allocator.js_malloc = QuickJSBinder::js_binder_malloc;
		godot_allocator.js_free = QuickJSBinder::js_binder_free;
		godot_allocator.js_realloc = QuickJSBinder::js_binder_realloc;
		godot_allocator.js_malloc_usable_size = nullptr;
		godot_object_class = nullptr;
		godot_reference_class = nullptr;
	}

	virtual ~QuickJSBinder() {
	}

	_FORCE_INLINE_ static QuickJSBinder *get_context_binder(JSContext *ctx) {
		return static_cast<QuickJSBinder *>(JS_GetContextOpaque(ctx));
	}

	virtual JavaScriptBinder *get_context_binder(void *p_context) override {
		return QuickJSBinder::get_context_binder((JSContext *)p_context);
	}

	virtual Thread::ID get_thread_id() const override { return thread_id; }

	_FORCE_INLINE_ static QuickJSBinder *get_runtime_binder(JSRuntime *rt) {
		return static_cast<QuickJSBinder *>(JS_GetMollocState(rt)->opaque);
	}

	_FORCE_INLINE_ static JavaScriptGCHandler *new_gc_handler(JSContext *ctx) {
		JavaScriptGCHandler *h = memnew(JavaScriptGCHandler);
		h->context = ctx;
		return h;
	}

	_FORCE_INLINE_ QuickJSBuiltinBinder &get_builtin_binder() { return builtin_binder; }

	_FORCE_INLINE_ JSClassID get_origin_class_id() { return godot_origin_class.class_id; }
	_FORCE_INLINE_ const ClassBindData get_origin_class() const { return godot_origin_class; }
	_FORCE_INLINE_ static JSClassID get_origin_class_id(JSContext *ctx) { return get_context_binder(ctx)->godot_origin_class.class_id; }

	virtual void initialize() override;
	virtual void uninitialize() override;
	virtual void language_finalize() override;
	virtual void frame() override;

	virtual JavaScriptGCHandler *alloc_object_binding_data(Object *p_object) override;
	virtual void free_object_binding_data(JavaScriptGCHandler *p_gc_handle) override;
	virtual void godot_refcount_incremented(JavaScriptGCHandler *p_gc_handle) override;
	virtual bool godot_refcount_decremented(JavaScriptGCHandler *p_gc_handle) override;

	static Error bind_gc_object(JSContext *ctx, JavaScriptGCHandler *data, Object *p_object);

	virtual Error eval_string(const String &p_source, EvalType type, const String &p_path, JavaScriptGCHandler &r_ret) override;
	virtual Error safe_eval_text(const String &p_source, EvalType type, const String &p_path, String &r_error, JavaScriptGCHandler &r_ret) override;

	virtual Error compile_to_bytecode(const String &p_code, const String &p_file, Vector<uint8_t> &r_bytecode) override;
	virtual Error load_bytecode(const Vector<uint8_t> &p_bytecode, const String &p_file, JavaScriptGCHandler *r_module) override;

	virtual const JavaScriptClassInfo *parse_javascript_class(const String &p_code, const String &p_path, bool ignore_cacehe, JavaScriptError *r_error) override;
	virtual const JavaScriptClassInfo *parse_javascript_class(const Vector<uint8_t> &p_bytecode, const String &p_path, bool ignore_cacehe, JavaScriptError *r_error) override;
	const JavaScriptClassInfo *parse_javascript_class_from_module(ModuleCache *p_module, const String &p_path, JavaScriptError *r_error);

	virtual JavaScriptGCHandler create_js_instance_for_godot_object(const JavaScriptClassInfo *p_class, Object *p_object) override;

	virtual bool get_instance_property(const JavaScriptGCHandler &p_object, const StringName &p_name, Variant &r_ret) override;
	virtual bool set_instance_property(const JavaScriptGCHandler &p_object, const StringName &p_name, const Variant &p_value) override;
	virtual bool has_method(const JavaScriptGCHandler &p_object, const StringName &p_name) override;
	virtual bool has_signal(const JavaScriptClassInfo *p_class, const StringName &p_signal) override;
	virtual bool validate(const String &p_code, const String &p_path, JavaScriptError *r_error) override;

	virtual Variant call_method(const JavaScriptGCHandler &p_object, const StringName &p_method, const Variant **p_args, int p_argcount, Callable::CallError &r_error) override;
	virtual Variant call(const JavaScriptGCHandler &p_fuction, const JavaScriptGCHandler &p_target, const Variant **p_args, int p_argcount, Callable::CallError &r_error) override;

#ifdef TOOLS_ENABLED
	virtual const Dictionary &get_modified_api() const override { return modified_api; }
#endif
};

#endif // QUICKJS_BINDING_HELPER_H
