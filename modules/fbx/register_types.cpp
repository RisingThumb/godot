/**************************************************************************/
/*  register_types.cpp                                                    */
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

#include "register_types.h"

#include "fbx_document.h"

#ifdef TOOLS_ENABLED
#include "editor/editor_scene_importer_fbx.h"

#include "core/config/project_settings.h"
#include "editor/editor_node.h"
#include "editor/editor_settings.h"

static void _editor_init() {
	Ref<EditorSceneFormatImporterUFBX> import_fbx;
	import_fbx.instantiate();
	ResourceImporterScene::add_importer(import_fbx);

	Ref<EditorSceneFormatImporterUFBX> importer;
	importer.instantiate();
	ResourceImporterScene::get_scene_singleton()->add_importer(importer);
}
#endif // TOOLS_ENABLED

#define FBX_REGISTER_DOCUMENT_EXTENSION(m_doc_ext_class) \
	Ref<m_doc_ext_class> extension_##m_doc_ext_class;    \
	extension_##m_doc_ext_class.instantiate();           \
	FBXDocument::register_fbx_document_extension(extension_##m_doc_ext_class);

void initialize_fbx_module(ModuleInitializationLevel p_level) {
	if (p_level == MODULE_INITIALIZATION_LEVEL_SCENE) {
		// glTF API available at runtime.
		GDREGISTER_CLASS(FBXAccessor);
		GDREGISTER_CLASS(FBXAnimation);
		GDREGISTER_CLASS(FBXBufferView);
		GDREGISTER_CLASS(FBXCamera);
		GDREGISTER_CLASS(FBXDocument);
		GDREGISTER_CLASS(FBXDocumentExtension);
		GDREGISTER_CLASS(FBXMesh);
		GDREGISTER_CLASS(FBXNode);
		GDREGISTER_CLASS(FBXSkeleton);
		GDREGISTER_CLASS(FBXSkin);
		GDREGISTER_CLASS(FBXState);
		GDREGISTER_CLASS(FBXTexture);
		GDREGISTER_CLASS(FBXTextureSampler);
	}

#ifdef TOOLS_ENABLED
	if (p_level == MODULE_INITIALIZATION_LEVEL_EDITOR) {
		// Editor-specific API.
		ClassDB::APIType prev_api = ClassDB::get_current_api();
		ClassDB::set_current_api(ClassDB::API_EDITOR);

		GDREGISTER_CLASS(EditorSceneFormatImporterUFBX);

		GDREGISTER_CLASS(EditorSceneFormatImporterUFBX);
		// Can't (a priori) run external app on these platforms.
		GLOBAL_DEF_RST("filesystem/import/blender/enabled.android", false);
		GLOBAL_DEF_RST("filesystem/import/blender/enabled.web", false);
		GLOBAL_DEF_RST("filesystem/import/fbx/enabled.android", false);
		GLOBAL_DEF_RST("filesystem/import/fbx/enabled.web", false);

		ClassDB::set_current_api(prev_api);
		EditorNode::add_init_callback(_editor_init);
	}

#endif // TOOLS_ENABLED
}

void uninitialize_fbx_module(ModuleInitializationLevel p_level) {
	if (p_level != MODULE_INITIALIZATION_LEVEL_SCENE) {
		return;
	}
	FBXDocument::unregister_all_fbx_document_extensions();
}
