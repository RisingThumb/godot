/// MIT License
///
/// Copyright (c) 2023 jolt_physics_sync
///
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to deal
/// in the Software without restriction, including without limitation the rights
/// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
/// copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
///
/// The above copyright notice and this permission notice shall be included in all
/// copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
/// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
/// SOFTWARE.

#include "jolt.h"

Jolt *Jolt::the_singleton = nullptr;

Jolt *Jolt::singleton() {
	return the_singleton;
}

void Jolt::_bind_methods() {
}

Jolt::Jolt() :
		Object() {
	if (the_singleton == nullptr) {
		the_singleton = this;
	} else {
		CRASH_NOW_MSG("You can't have two Jolt singletons at the same time.");
	}
}

Jolt::~Jolt() {
	if (the_singleton == this) {
		the_singleton = nullptr;
	}
}
