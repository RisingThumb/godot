// Copyright 2015 The Crashpad Authors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef CRASHPAD_COMPAT_WIN_STRINGS_H_
#define CRASHPAD_COMPAT_WIN_STRINGS_H_

#ifdef __cplusplus
extern "C" {
#endif

int strcasecmp(const char* s1, const char* s2);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // CRASHPAD_COMPAT_WIN_STRINGS_H_
