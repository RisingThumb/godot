/* -*- Mode: C; tab-width: 4 -*-
 *
 * Copyright (c) 2003-2004 Apple Computer, Inc. All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */



#ifndef _Firewall_h

#define _Firewall_h





#include "CommonServices.h"

#include "DebugServices.h"





#if defined(__cplusplus)

extern "C"

{

#endif





OSStatus

mDNSAddToFirewall

		(

		LPWSTR	executable,

		LPWSTR	name

		);


BOOL
mDNSIsFileAndPrintSharingEnabled( BOOL * retry );





#if defined(__cplusplus)

}

#endif





#endif

