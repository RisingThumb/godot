/*************************************************************************/
/*  fbx_material.cpp                                                     */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                      https://godotengine.org                          */
/*************************************************************************/
/* Copyright (c) 2007-2020 Juan Linietsky, Ariel Manzur.                 */
/* Copyright (c) 2014-2020 Godot Engine contributors (cf. AUTHORS.md).   */
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

#include "fbx_material.h"

#include "scene/resources/material.h"
#include "scene/resources/texture.h"

String FBXMaterial::get_material_name() const {
	return material_name;
}

void FBXMaterial::set_imported_material(const Assimp::FBX::Material *p_material) {
	material = p_material;
}

void FBXMaterial::add_search_string(String p_filename, String p_current_directory, String search_directory, Vector<String> &texture_search_paths) {
	if (search_directory.empty()) {
		texture_search_paths.push_back(p_current_directory.get_base_dir().plus_file(p_filename));
	} else {
		texture_search_paths.push_back(p_current_directory.get_base_dir().plus_file(search_directory + "/" + p_filename));
		texture_search_paths.push_back(p_current_directory.get_base_dir().plus_file("../" + search_directory + "/" + p_filename));
	}
}

String find_file(const String &p_base, const String &p_file_to_find) {
	_Directory dir;
	dir.open(p_base);

	dir.list_dir_begin();
	String n = dir.get_next();
	while (n != String()) {
		if (n == "." || n == "..") {
			n = dir.get_next();
			continue;
		}
		if (dir.current_is_dir()) {
			// Don't use `path_to` or the returned path will be wrong.
			const String f = find_file(p_base + "/" + n, p_file_to_find);
			if (f != "") {
				return f;
			}
		} else if (n == p_file_to_find) {
			return p_base + "/" + n;
		}
		n = dir.get_next();
	}
	dir.list_dir_end();

	return String();
}

// fbx will not give us good path information and let's not regex them to fix them
// no relative paths are in fbx generally they have a rel field but it's populated incorrectly by the SDK.
String FBXMaterial::find_texture_path_by_filename(const String p_filename, const String p_current_directory) {
	_Directory dir;
	Vector<String> paths;
	add_search_string(p_filename, p_current_directory, "", paths);
	add_search_string(p_filename, p_current_directory, "texture", paths);
	add_search_string(p_filename, p_current_directory, "textures", paths);
	add_search_string(p_filename, p_current_directory, "Textures", paths);
	add_search_string(p_filename, p_current_directory, "materials", paths);
	add_search_string(p_filename, p_current_directory, "mats", paths);
	add_search_string(p_filename, p_current_directory, "pictures", paths);
	add_search_string(p_filename, p_current_directory, "images", paths);

	for (int i = 0; i < paths.size(); i++) {
		if (dir.file_exists(paths[i])) {
			return paths[i];
		}
	}

	// We were not able to find the texture in the common locations,
	// try to find it into the project globally.
	// The common textures can be stored into one of those folders:
	// res://asset
	// res://texture
	// res://material
	// res://mat
	// res://image
	// res://picture
	//
	// Note the folders can also be called with custom names, like:
	// res://my_assets
	// since the keyword `asset` is into the directory name the textures will be
	// searched there too.

	dir.open("res://");
	dir.list_dir_begin();
	String n = dir.get_next();
	while (n != String()) {
		if (n == "." || n == "..") {
			n = dir.get_next();
			continue;
		}
		if (dir.current_is_dir()) {
			const String lower_n = n.to_lower();
			if (
					// Don't need to use plural.
					lower_n.find("asset") >= 0 ||
					lower_n.find("texture") >= 0 ||
					lower_n.find("material") >= 0 ||
					lower_n.find("mat") >= 0 ||
					lower_n.find("image") >= 0 ||
					lower_n.find("picture") >= 0) {
				// Don't use `path_to` or the returned path will be wrong.
				const String f = find_file(String("res://") + n, p_filename);
				if (f != "") {
					return f;
				}
			}
		}
		n = dir.get_next();
	}
	dir.list_dir_end();

	return "";
}

FBXMaterial::MaterialInfo FBXMaterial::extract_material_info(const Assimp::FBX::Material *material) const {
	MaterialInfo mat_info;

	// TODO Layered textures are a collection on textures stored into an array.
	// Extract layered textures is not yet supported. Per each texture in the
	// layered texture array you want to use the below method to extract those.

	for (std::pair<std::string, const Assimp::FBX::Texture *> texture : material->Textures()) {
		const std::string &fbx_mapping_name = texture.first;

		if (fbx_feature_mapping_desc.count(fbx_mapping_name) > 0) {
			// This is a feature not a normal texture.
			mat_info.features.push_back(fbx_feature_mapping_desc.at(fbx_mapping_name));
			continue;
		}

		ERR_CONTINUE_MSG(fbx_texture_mapping_desc.count(fbx_mapping_name) <= 0, "This FBX has a material with mapping name: " + String(fbx_mapping_name.c_str()) + " which is not yet supported by this importer. Consider open an issue so we can support it.");

		const String absoulte_fbx_file_path = texture.second->FileName().c_str();
		const String file_extension = absoulte_fbx_file_path.get_extension().to_upper();

		const String file_extension_uppercase = file_extension.to_upper();

		// TODO: we don't support EMBED for DDS and TGA.
		ERR_CONTINUE_MSG(
				file_extension_uppercase != "PNG" &&
						file_extension_uppercase != "JPEG" &&
						file_extension_uppercase != "JPG" &&
						file_extension_uppercase != "TGA" &&
						file_extension_uppercase != "WEBP" &&
						file_extension_uppercase != "DDS",
				"The FBX file contains a texture with an unrecognized extension: " + file_extension_uppercase);

		const String texture_name = absoulte_fbx_file_path.get_file();
		const SpatialMaterial::TextureParam mapping_mode = fbx_texture_mapping_desc.at(fbx_mapping_name);

		TextureFileMapping file_mapping;
		file_mapping.map_mode = mapping_mode;
		file_mapping.name = texture_name;
		file_mapping.texture = texture.second;
		mat_info.textures.push_back(file_mapping);

		// Make sure to active the various features.
		switch (mapping_mode) {
			case SpatialMaterial::TextureParam::TEXTURE_ALBEDO:
			case SpatialMaterial::TextureParam::TEXTURE_METALLIC:
			case SpatialMaterial::TextureParam::TEXTURE_ROUGHNESS:
			case SpatialMaterial::TextureParam::TEXTURE_FLOWMAP:
			case SpatialMaterial::TextureParam::TEXTURE_REFRACTION:
			case SpatialMaterial::TextureParam::TEXTURE_MAX:
				// No features required.
				break;
			case SpatialMaterial::TextureParam::TEXTURE_EMISSION:
				mat_info.features.push_back(SpatialMaterial::Feature::FEATURE_EMISSION);
				break;
			case SpatialMaterial::TextureParam::TEXTURE_NORMAL:
				mat_info.features.push_back(SpatialMaterial::Feature::FEATURE_NORMAL_MAPPING);
				break;
			case SpatialMaterial::TextureParam::TEXTURE_RIM:
				mat_info.features.push_back(SpatialMaterial::Feature::FEATURE_RIM);
				break;
			case SpatialMaterial::TextureParam::TEXTURE_CLEARCOAT:
				mat_info.features.push_back(SpatialMaterial::Feature::FEATURE_CLEARCOAT);
				break;
			case SpatialMaterial::TextureParam::TEXTURE_AMBIENT_OCCLUSION:
				mat_info.features.push_back(SpatialMaterial::Feature::FEATURE_AMBIENT_OCCLUSION);
				break;
			case SpatialMaterial::TextureParam::TEXTURE_DEPTH:
				mat_info.features.push_back(SpatialMaterial::Feature::FEATURE_DEPTH_MAPPING);
				break;
			case SpatialMaterial::TextureParam::TEXTURE_SUBSURFACE_SCATTERING:
				mat_info.features.push_back(SpatialMaterial::Feature::FEATURE_SUBSURACE_SCATTERING);
				break;
			case SpatialMaterial::TextureParam::TEXTURE_TRANSMISSION:
				mat_info.features.push_back(SpatialMaterial::Feature::FEATURE_TRANSMISSION);
				break;
			case SpatialMaterial::TextureParam::TEXTURE_DETAIL_ALBEDO:
			case SpatialMaterial::TextureParam::TEXTURE_DETAIL_MASK:
			case SpatialMaterial::TextureParam::TEXTURE_DETAIL_NORMAL:
				mat_info.features.push_back(SpatialMaterial::Feature::FEATURE_DETAIL);
				break;
		}
	}

	return mat_info;
}

template <class T>
T extract_from_prop(Assimp::FBX::PropertyPtr prop, const T &p_default, const std::string &p_name, const String &p_type) {
	ERR_FAIL_COND_V_MSG(prop == nullptr, p_default, "invalid property passed to extractor");
	const Assimp::FBX::TypedProperty<T> *val = dynamic_cast<const Assimp::FBX::TypedProperty<T> *>(prop);

	ERR_FAIL_COND_V_MSG(val == nullptr, p_default, "The FBX is corrupted, the property `" + String(p_name.c_str()) + "` is a `" + String(typeid(*prop).name()) + "` but should be a " + p_type);
	// Make sure to not lost any eventual opacity.
	return val->Value();
}

Ref<SpatialMaterial> FBXMaterial::import_material(ImportState &state) {

	ERR_FAIL_COND_V(material == nullptr, nullptr);

	const String p_fbx_current_directory = state.path;

	Ref<SpatialMaterial> spatial_material;

	// read the material file
	// is material two sided
	// read material name
	print_verbose("[material] material name: " + ImportUtils::FBXNodeToName(material->Name()));
	material_name = ImportUtils::FBXNodeToName(material->Name());

	// Extract info.
	MaterialInfo material_info = extract_material_info(material);

	// Extract other parameters info.
	for (Assimp::FBX::LazyPropertyMap::value_type iter : material->Props()->GetLazyProperties()) {
		const std::string name = iter.first;
		//const Assimp::FBX::ElementPtr element = iter.second;

		if (name.empty()) {
			continue;
		}

		PropertyDesc desc = PROPERTY_DESC_NOT_FOUND;
		if (fbx_properties_desc.count(name) > 0) {
			desc = fbx_properties_desc.at(name);
		}

		if (desc == PROPERTY_DESC_IGNORE) {
			print_verbose("The FBX material parameter: `" + String(name.c_str()) + "` is ignored.");
			continue;
		} else {
			print_verbose("FBX Material parameter: " + String(name.c_str()));
		}

		if (desc == PROPERTY_DESC_NOT_FOUND) {
			continue;
		}

		ERR_CONTINUE_MSG(desc == PROPERTY_DESC_NOT_FOUND, "The FBX material parameter: `" + String(name.c_str()) + "` was not recognized. Please open an issue so we can add the support to it.");

		Assimp::FBX::PropertyPtr prop = material->Props()->Get(name);

		//Assimp::FBX::PropertyPtr prop = prop.second.

		if (prop == nullptr) {
			continue;
		}
		ERR_CONTINUE_MSG(prop == nullptr, "This file may be corrupted because is not possible to extract the material parameter: " + String(name.c_str()));

		if (spatial_material.is_null()) {
			// Done here so if no data no material is created.
			spatial_material.instance();
		}

		switch (desc) {
			case PROPERTY_DESC_ALBEDO_COLOR: {
				const Vector3 color = extract_from_prop(prop, Vector3(0, 0, 0), name, "Vector3");
				// Make sure to not lost any eventual opacity.
				Color c = spatial_material->get_albedo();
				c[0] = color[0];
				c[1] = color[1];
				c[2] = color[2];
				spatial_material->set_albedo(c);
			} break;
			case PROPERTY_DESC_TRANSPARENT: {
				const real_t opacity = extract_from_prop(prop, 1.0f, name, "float");
				if (opacity < (1.0 - CMP_EPSILON)) {
					Color c = spatial_material->get_albedo();
					c[3] = opacity;
					spatial_material->set_albedo(c);
					material_info.features.push_back(SpatialMaterial::Feature::FEATURE_TRANSPARENT);
					spatial_material->set_depth_draw_mode(SpatialMaterial::DEPTH_DRAW_ALPHA_OPAQUE_PREPASS);
				}
			} break;
			case PROPERTY_DESC_METALLIC: {
				spatial_material->set_metallic(std::min(1.0f, extract_from_prop(prop, 1.0f, name, "float")));
			} break;
			case PROPERTY_DESC_ROUGHNESS: {
				spatial_material->set_roughness(std::min(1.0f, extract_from_prop(prop, 1.0f, name, "float")));
			} break;
			case PROPERTY_DESC_COAT: {
				spatial_material->set_clearcoat(extract_from_prop(prop, 1.0f, name, "float"));
				material_info.features.push_back(SpatialMaterial::Feature::FEATURE_CLEARCOAT);
			} break;
			case PROPERTY_DESC_COAT_ROUGHNESS: {
				spatial_material->set_clearcoat_gloss(1.0 - extract_from_prop(prop, 0.5f, name, "float"));
				material_info.features.push_back(SpatialMaterial::Feature::FEATURE_CLEARCOAT);
			} break;
			case PROPERTY_DESC_EMISSIVE: {
				const real_t emissive = extract_from_prop(prop, 0.0f, name, "float");
				if (emissive > CMP_EPSILON) {
					spatial_material->set_emission_energy(emissive);
					material_info.features.push_back(SpatialMaterial::Feature::FEATURE_EMISSION);
				}
			} break;
			case PROPERTY_DESC_EMISSIVE_COLOR: {
				const Vector3 color = extract_from_prop(prop, Vector3(0, 0, 0), name, "Vector3");
				Color c;
				c[0] = color[0];
				c[1] = color[1];
				c[2] = color[2];
				spatial_material->set_emission(c);
			} break;
			case PROPERTY_DESC_NOT_FOUND:
			case PROPERTY_DESC_IGNORE:
				// Already checked, can't happen.
				CRASH_NOW();
				break;
		}
	}

	// Set the material features.
	for (int x = 0; x < material_info.features.size(); x++) {
		if (spatial_material.is_null()) {
			// Done here so if no textures no material is created.
			spatial_material.instance();
		}
		spatial_material->set_feature(material_info.features[x], true);
	}

	// Set the textures.
	for (int x = 0; x < material_info.textures.size(); x++) {
		TextureFileMapping mapping = material_info.textures[x];
		Ref<Texture> texture;
		print_verbose("texture mapping name: " + mapping.name);

		if (state.cached_image_searches.has(mapping.name)) {
			texture = state.cached_image_searches[mapping.name];
		} else {
			String path = find_texture_path_by_filename(mapping.name, p_fbx_current_directory);
			if (!path.empty()) {
				Error err;
				Ref<Texture> image_texture = ResourceLoader::load(path, "Texture", false, &err);

				ERR_CONTINUE_MSG(err != OK, "unable to import image file not loaded yet: " + path);
				ERR_CONTINUE(image_texture == NULL || image_texture.is_null());

				texture = image_texture;
				state.cached_image_searches.insert(mapping.name, texture);
				print_verbose("Created texture from loaded image file.");

			} else if (mapping.texture != nullptr && mapping.texture->Media() != nullptr) {
				// This is an embedded texture. Extract it.
				Ref<Image> image;
				image.instance();

				const String extension = mapping.name.get_extension().to_upper();
				if (extension == "PNG") {

					// The stored file is a PNG.
					image = Image::_png_mem_loader_func(mapping.texture->Media()->Content(), mapping.texture->Media()->ContentLength());
					ERR_CONTINUE_MSG(image.is_valid() == false, "FBX Embedded PNG image load fail.");

				} else if (
						extension == "JPEG" ||
						extension == "JPG") {

					// The stored file is a JPEG.
					image = Image::_jpg_mem_loader_func(mapping.texture->Media()->Content(), mapping.texture->Media()->ContentLength());
					ERR_CONTINUE_MSG(image.is_valid() == false, "FBX Embedded JPEG image load fail.");

				} else if (extension == "TGA") {

					// The stored file is a TGA.
					//image = Image::_tga_mem_loader_func(mapping.texture->Media()->Content(), mapping.texture->Media()->ContentLength());
					//ERR_CONTINUE_MSG(image.is_valid() == false, "FBX Embedded TGA image load fail.");

				} else if (extension == "WEBP") {

					// The stored file is a WEBP.
					image = Image::_webp_mem_loader_func(mapping.texture->Media()->Content(), mapping.texture->Media()->ContentLength());
					ERR_CONTINUE_MSG(image.is_valid() == false, "FBX Embedded WEBP image load fail.");

					// } else if (extension == "DDS") {
					// 	// In this moment is not possible to extract a DDS from a buffer, TODO consider add it to godot. See `textureloader_dds.cpp::load().
					// 	// The stored file is a DDS.
				} else {
					ERR_CONTINUE_MSG(true, "The embedded image with extension: " + extension + " is not yet supported. Open an issue please.");
				}

				Ref<ImageTexture> image_texture;
				image_texture.instance();
				image_texture->create_from_image(image);

				const int32_t flags = Texture::FLAGS_DEFAULT;
				image_texture->set_flags(flags);

				texture = image_texture;
				state.cached_image_searches[mapping.name] = texture;
				print_verbose("Created texture from embedded image.");
			} else {
				ERR_CONTINUE_MSG(true, "The FBX texture, with name: `" + mapping.name + "`, is not found into the project nor is stored as embedded file. Make sure to insert the texture as embedded file or into the project, then reimport.");
			}
		}
		if (spatial_material.is_null()) {
			// Done here so if no textures no material is created.
			spatial_material.instance();
		}

		switch (mapping.map_mode) {
			case SpatialMaterial::TextureParam::TEXTURE_METALLIC:
				if (mapping.name.to_lower().find("ser") >= 0) {
					// SER shader.
					spatial_material->set_metallic_texture_channel(SpatialMaterial::TextureChannel::TEXTURE_CHANNEL_RED);
				} else {
					// Use grayscale as default.
					spatial_material->set_metallic_texture_channel(SpatialMaterial::TextureChannel::TEXTURE_CHANNEL_GRAYSCALE);
				}
				break;
			case SpatialMaterial::TextureParam::TEXTURE_ROUGHNESS:
				if (mapping.name.to_lower().find("ser") >= 0) {
					// SER shader.
					spatial_material->set_roughness_texture_channel(SpatialMaterial::TextureChannel::TEXTURE_CHANNEL_BLUE);
				} else {
					// Use grayscale as default.
					spatial_material->set_roughness_texture_channel(SpatialMaterial::TextureChannel::TEXTURE_CHANNEL_GRAYSCALE);
				}
				break;
			case SpatialMaterial::TextureParam::TEXTURE_AMBIENT_OCCLUSION:
				// Use grayscale as default.
				spatial_material->set_ao_texture_channel(SpatialMaterial::TextureChannel::TEXTURE_CHANNEL_GRAYSCALE);
				break;
			case SpatialMaterial::TextureParam::TEXTURE_REFRACTION:
				// Use grayscale as default.
				spatial_material->set_refraction_texture_channel(SpatialMaterial::TextureChannel::TEXTURE_CHANNEL_GRAYSCALE);
				break;
			default:
				// Nothing to do.
				break;
		}

		spatial_material->set_texture(mapping.map_mode, texture);
	}

	if (spatial_material.is_valid()) {
		spatial_material->set_name(material_name);
	}

	return spatial_material;
}
