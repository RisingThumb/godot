/*************************************************************************/
/*  navigation_polygon_editor_plugin.h                                   */
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

#ifndef NAVIGATION_POLYGON_EDITOR_PLUGIN_H
#define NAVIGATION_POLYGON_EDITOR_PLUGIN_H

#include "editor/plugins/abstract_polygon_2d_editor.h"
#include "scene/2d/navigation_region_2d.h"

class NavigationPolygonEditor : public AbstractPolygon2DEditor {
	GDCLASS(NavigationPolygonEditor, AbstractPolygon2DEditor);

	NavigationRegion2D *node = nullptr;

	Ref<NavigationPolygon> _ensure_navpoly() const;

protected:
	virtual Node2D *_get_node() const override;
	virtual void _set_node(Node *p_polygon) override;

	virtual int _get_polygon_count() const override;
	virtual Variant _get_polygon(int p_idx) const override;
	virtual void _set_polygon(int p_idx, const Variant &p_polygon) const override;

	virtual void _action_add_polygon(const Variant &p_polygon) override;
	virtual void _action_remove_polygon(int p_idx) override;
	virtual void _action_set_polygon(int p_idx, const Variant &p_previous, const Variant &p_polygon) override;

	virtual bool _has_resource() const override;
	virtual bool _resource_is_foreign() const override;
	virtual void _create_resource() override;

public:
	NavigationPolygonEditor();
};

class NavigationPolygonEditorPlugin : public AbstractPolygon2DEditorPlugin {
	GDCLASS(NavigationPolygonEditorPlugin, AbstractPolygon2DEditorPlugin);

public:
	NavigationPolygonEditorPlugin();
};

#endif // NAVIGATION_POLYGON_EDITOR_PLUGIN_H
