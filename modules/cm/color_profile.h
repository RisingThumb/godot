/*************************************************************************/
/*  color_profile.h                                                      */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                      https://godotengine.org                          */
/*************************************************************************/
/* Copyright (c) 2007-2021 Juan Linietsky, Ariel Manzur.                 */
/* Copyright (c) 2014-2021 Godot Engine contributors (cf. AUTHORS.md).   */
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

#ifndef COLOR_PROFILE_H
#define COLOR_PROFILE_H

#include "core/object/ref_counted.h"
#include "thirdparty/lcms2/lcms2.h"

class ColorProfile : public RefCounted {
	GDCLASS(ColorProfile, RefCounted)

	bool profile_valid;
	cmsHPROFILE profile;

	void _set_profile(cmsHPROFILE p_profile);

public:
	enum Predef {
		PREDEF_SRGB,
	};

	class Handle {
		friend class ColorProfile;
		friend class ColorTransform;

		cmsHPROFILE profile;

		Handle() { profile = nullptr; }
		Handle(cmsHPROFILE p_profile) { profile = p_profile; }
	};

	bool is_valid();
	Handle get_profile_handle();

	void load_predef(Predef p_predef);
	void load_from_file(const String &p_file);

	ColorProfile();
	ColorProfile(Predef p_predef);
	ColorProfile(const String &p_file);
	~ColorProfile();
};

#endif // COLOR_PROFILE_H
