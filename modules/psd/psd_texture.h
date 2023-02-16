/**************************************************************************/
/*  psd_texture.h                                                         */
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

#ifndef PSD_TEXTURE_H
#define PSD_TEXTURE_H

#include "core/io/resource_loader.h"

#include "Psd.h"

#include "PsdMallocAllocator.h"

#include "PsdChannel.h"
#include "PsdChannelType.h"
#include "PsdColorMode.h"
#include "PsdDocument.h"
#include "PsdExport.h"
#include "PsdExportDocument.h"
#include "PsdImageDataSection.h"
#include "PsdImageResourcesSection.h"
#include "PsdInterleave.h"
#include "PsdLayer.h"
#include "PsdLayerCanvasCopy.h"
#include "PsdLayerMask.h"
#include "PsdLayerMaskSection.h"
#include "PsdLayerType.h"
#include "PsdParseDocument.h"
#include "PsdParseImageDataSection.h"
#include "PsdParseImageResourcesSection.h"
#include "PsdParseLayerMaskSection.h"
#include "PsdPlanarImage.h"
#include "PsdVectorMask.h"

#include "core/io/image.h"
#include "scene/resources/texture.h"

#include <sstream>
#include <string>

#include "core/io/image_loader.h"

class PSDTexture;

class PSDTexture : public ImageFormatLoader {
	GDCLASS(PSDTexture, ImageFormatLoader);
	bool cropToCanvas = true;

	HashMap<String, Ref<Image>> layers;

	typedef enum {
		KRA,
		PSD
	} IMPORT_TYPE;

	typedef enum {
		MONOCHROME,
		RGB,
		RGBA
	} COLOR_SPACE_NAME;

protected:
	static void _bind_methods();

public:
	virtual Error load_image(Ref<Image> p_image, Ref<FileAccess> f, BitField<ImageFormatLoader::LoaderFlags> p_flags, float p_scale) override;
	virtual void get_recognized_extensions(List<String> *p_extensions) const override;

	PSDTexture();
	virtual ~PSDTexture();
};

#endif // PSD_TEXTURE_H
