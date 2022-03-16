/*************************************************************************/
/*  openxr_select_interaction_profile_dialog.h                           */
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

#ifndef OPENXR_SELECT_INTERACTION_PROFILE_DIALOG_H
#define OPENXR_SELECT_INTERACTION_PROFILE_DIALOG_H

#include "../action_map/openxr_defs.h"
#include "scene/gui/box_container.h"
#include "scene/gui/button.h"
#include "scene/gui/dialogs.h"
#include "scene/gui/label.h"
#include "scene/gui/line_edit.h"
#include "scene/gui/scroll_container.h"
#include "scene/gui/separator.h"
#include "scene/gui/text_edit.h"

class OpenXRSelectInteractionProfileDialog : public ConfirmationDialog {
	GDCLASS(OpenXRSelectInteractionProfileDialog, ConfirmationDialog);

private:
	String selected_interaction_profile;
	Dictionary ip_buttons;

	VBoxContainer *main_vb = nullptr;
	ScrollContainer *scroll = nullptr;

protected:
	static void _bind_methods();
	void _notification(int p_what);

public:
	void _on_select_interaction_profile(const String p_interaction_profile);
	void open(PackedStringArray p_do_not_include);
	virtual void ok_pressed() override;

	OpenXRSelectInteractionProfileDialog();
};

#endif // !OPENXR_SELECT_INTERACTION_PROFILE_DIALOG_H
