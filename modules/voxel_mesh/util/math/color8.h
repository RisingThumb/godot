/**************************************************************************/
/*  color8.h                                                              */
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

#ifndef COLOR8_H
#define COLOR8_H

#include <core/math/color.h>

struct Color8 {
	union {
		struct {
			uint8_t r;
			uint8_t g;
			uint8_t b;
			uint8_t a;
		};
		uint8_t components[4];
	};

	Color8() :
			r(0), g(0), b(0), a(0) {}

	Color8(uint8_t p_r, uint8_t p_g, uint8_t p_b, uint8_t p_a) :
			r(p_r), g(p_g), b(p_b), a(p_a) {}

	Color8(Color c) {
		r = c.r * 255;
		g = c.g * 255;
		b = c.b * 255;
		a = c.a * 255;
	}

	static inline Color8 from_u8(uint8_t v) {
		// rrggbbaa
		return Color8(v >> 6, ((v >> 4) & 3), ((v >> 2) & 3), v & 3);
	}

	static inline Color8 from_u16(uint16_t v) {
		// rrrrgggg bbbbaaaa 🐐
		return Color8(v >> 12, ((v >> 8) & 0xf), ((v >> 4) & 0xf), v & 0xf);
	}

	static inline Color8 from_u32(uint32_t c) {
		// rrrrrrrr gggggggg bbbbbbbb aaaaaaaa
		return Color8(c >> 24, (c >> 16) & 0xff, (c >> 8) & 0xff, c & 0xff);
	}

	inline uint8_t to_u8() const {
		// Lossy
		return ((r >> 6) << 6) | ((g >> 6) << 4) | ((b >> 6) << 2) | (a >> 6);
	}

	inline uint16_t to_u16() const {
		// Lossy
		return ((r >> 4) << 12) | ((g >> 4) << 8) | ((b >> 4) << 4) | (a >> 4);
	}

	inline uint32_t to_u32() const {
		return (r << 24) | (g << 16) | (b << 8) | a;
	}

	operator Color() const {
		return Color(r / 255.f, g / 255.f, b / 255.f, a / 255.f);
	}
};

#endif // COLOR8_H
