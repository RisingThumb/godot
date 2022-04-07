/*************************************************************************/
/*  animation_library.h                                                  */
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

#ifndef ANIMATION_LIBRARY_H
#define ANIMATION_LIBRARY_H

#include "core/variant/typed_array.h"
#include "scene/resources/animation.h"

class AnimationLibrary : public Resource {
	GDCLASS(AnimationLibrary, Resource)

	void _set_data(const Dictionary &p_data);
	Dictionary _get_data() const;

	TypedArray<StringName> _get_animation_list() const;

	friend class AnimationPlayer; //for faster access
	Map<StringName, Ref<Animation>> animations;

protected:
	static void _bind_methods();

public:
	Error add_animation(const StringName &p_name, const Ref<Animation> &p_animation);
	void remove_animation(const StringName &p_name);
	void rename_animation(const StringName &p_name, const StringName &p_new_name);
	bool has_animation(const StringName &p_name) const;
	Ref<Animation> get_animation(const StringName &p_name) const;
	void get_animation_list(List<StringName> *p_animations) const;

	AnimationLibrary();
};

#endif // ANIMATIONLIBRARY_H
