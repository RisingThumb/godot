/*************************************************************************/
/*  godot_window_delegate.mm                                             */
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

#include "godot_window_delegate.h"

#include "display_server_macos.h"
#include "godot_button_view.h"
#include "godot_window.h"

@implementation GodotWindowDelegate

- (void)setWindowID:(DisplayServer::WindowID)wid {
	window_id = wid;
}

- (BOOL)windowShouldClose:(id)sender {
	DisplayServerMacOS *ds = (DisplayServerMacOS *)DisplayServer::get_singleton();
	if (!ds || !ds->has_window(window_id)) {
		return YES;
	}

	ds->send_window_event(ds->get_window(window_id), DisplayServerMacOS::WINDOW_EVENT_CLOSE_REQUEST);
	return NO;
}

- (void)windowWillClose:(NSNotification *)notification {
	DisplayServerMacOS *ds = (DisplayServerMacOS *)DisplayServer::get_singleton();
	if (!ds || !ds->has_window(window_id)) {
		return;
	}

	ds->popup_close(window_id);

	DisplayServerMacOS::WindowData &wd = ds->get_window(window_id);
	while (wd.transient_children.size()) {
		ds->window_set_transient(*wd.transient_children.begin(), DisplayServerMacOS::INVALID_WINDOW_ID);
	}

	if (wd.transient_parent != DisplayServerMacOS::INVALID_WINDOW_ID) {
		ds->window_set_transient(window_id, DisplayServerMacOS::INVALID_WINDOW_ID);
	}

	ds->window_destroy(window_id);
}

- (NSArray<NSWindow *> *)customWindowsToEnterFullScreenForWindow:(NSWindow *)window {
	DisplayServerMacOS *ds = (DisplayServerMacOS *)DisplayServer::get_singleton();
	if (!ds || !ds->has_window(window_id)) {
		return nullptr;
	}

	old_frame = [window frame];
	old_style_mask = [window styleMask];

	NSMutableArray<NSWindow *> *windows = [[NSMutableArray alloc] init];
	[windows addObject:window];

	return windows;
}

- (void)window:(NSWindow *)window startCustomAnimationToEnterFullScreenWithDuration:(NSTimeInterval)duration {
	[(GodotWindow *)window setAnimDuration:duration];
	[window setFrame:[[window screen] frame] display:YES animate:YES];
}

- (void)windowDidEnterFullScreen:(NSNotification *)notification {
	DisplayServerMacOS *ds = (DisplayServerMacOS *)DisplayServer::get_singleton();
	if (!ds || !ds->has_window(window_id)) {
		return;
	}

	DisplayServerMacOS::WindowData &wd = ds->get_window(window_id);
	wd.fullscreen = true;

	// Reset window size limits.
	[wd.window_object setContentMinSize:NSMakeSize(0, 0)];
	[wd.window_object setContentMaxSize:NSMakeSize(FLT_MAX, FLT_MAX)];
	[(GodotWindow *)wd.window_object setAnimDuration:-1.0f];

	// Reset custom window buttons.
	if ([wd.window_object styleMask] & NSWindowStyleMaskFullSizeContentView) {
		ds->window_set_custom_window_buttons(wd, false);
	}

	// Force window resize event.
	[self windowDidResize:notification];
	ds->send_window_event(wd, DisplayServerMacOS::WINDOW_EVENT_TITLEBAR_CHANGE);
}

- (NSArray<NSWindow *> *)customWindowsToExitFullScreenForWindow:(NSWindow *)window {
	DisplayServerMacOS *ds = (DisplayServerMacOS *)DisplayServer::get_singleton();
	if (!ds || !ds->has_window(window_id)) {
		return nullptr;
	}

	DisplayServerMacOS::WindowData &wd = ds->get_window(window_id);

	// Restore custom window buttons.
	if ([wd.window_object styleMask] & NSWindowStyleMaskFullSizeContentView) {
		ds->window_set_custom_window_buttons(wd, true);
	}

	ds->send_window_event(wd, DisplayServerMacOS::WINDOW_EVENT_TITLEBAR_CHANGE);

	NSMutableArray<NSWindow *> *windows = [[NSMutableArray alloc] init];
	[windows addObject:wd.window_object];
	return windows;
}

- (void)window:(NSWindow *)window startCustomAnimationToExitFullScreenWithDuration:(NSTimeInterval)duration {
	[(GodotWindow *)window setAnimDuration:duration];
	[window setStyleMask:old_style_mask];
	[window setFrame:old_frame display:YES animate:YES];
}

- (void)windowDidExitFullScreen:(NSNotification *)notification {
	DisplayServerMacOS *ds = (DisplayServerMacOS *)DisplayServer::get_singleton();
	if (!ds || !ds->has_window(window_id)) {
		return;
	}

	DisplayServerMacOS::WindowData &wd = ds->get_window(window_id);
	wd.fullscreen = false;

	[(GodotWindow *)wd.window_object setAnimDuration:-1.0f];

	// Set window size limits.
	const float scale = ds->screen_get_max_scale();
	if (wd.min_size != Size2i()) {
		Size2i size = wd.min_size / scale;
		[wd.window_object setContentMinSize:NSMakeSize(size.x, size.y)];
	}
	if (wd.max_size != Size2i()) {
		Size2i size = wd.max_size / scale;
		[wd.window_object setContentMaxSize:NSMakeSize(size.x, size.y)];
	}

	// Restore resizability state.
	if (wd.resize_disabled) {
		[wd.window_object setStyleMask:[wd.window_object styleMask] & ~NSWindowStyleMaskResizable];
	}

	// Restore on-top state.
	if (wd.on_top) {
		[wd.window_object setLevel:NSFloatingWindowLevel];
	}

	// Force window resize event.
	[self windowDidResize:notification];
}

- (void)windowDidChangeBackingProperties:(NSNotification *)notification {
	DisplayServerMacOS *ds = (DisplayServerMacOS *)DisplayServer::get_singleton();
	if (!ds || !ds->has_window(window_id)) {
		return;
	}

	DisplayServerMacOS::WindowData &wd = ds->get_window(window_id);

	CGFloat new_scale_factor = [wd.window_object backingScaleFactor];
	CGFloat old_scale_factor = [[[notification userInfo] objectForKey:@"NSBackingPropertyOldScaleFactorKey"] doubleValue];

	if (new_scale_factor != old_scale_factor) {
		// Set new display scale and window size.
		const float scale = ds->screen_get_max_scale();
		const NSRect content_rect = [wd.window_view frame];

		wd.size.width = content_rect.size.width * scale;
		wd.size.height = content_rect.size.height * scale;

		ds->send_window_event(wd, DisplayServerMacOS::WINDOW_EVENT_DPI_CHANGE);

		CALayer *layer = [wd.window_view layer];
		if (layer) {
			layer.contentsScale = scale;
		}

		//Force window resize event
		[self windowDidResize:notification];
	}
}

- (void)windowWillStartLiveResize:(NSNotification *)notification {
	DisplayServerMacOS *ds = (DisplayServerMacOS *)DisplayServer::get_singleton();
	if (ds && ds->has_window(window_id)) {
		DisplayServerMacOS::WindowData &wd = ds->get_window(window_id);
		wd.last_frame_rect = [wd.window_object frame];
		ds->set_is_resizing(true);
	}
}

- (void)windowDidEndLiveResize:(NSNotification *)notification {
	DisplayServerMacOS *ds = (DisplayServerMacOS *)DisplayServer::get_singleton();
	if (ds) {
		ds->set_is_resizing(false);
	}
}

- (void)windowDidResize:(NSNotification *)notification {
	DisplayServerMacOS *ds = (DisplayServerMacOS *)DisplayServer::get_singleton();
	if (!ds || !ds->has_window(window_id)) {
		return;
	}

	DisplayServerMacOS::WindowData &wd = ds->get_window(window_id);
	const NSRect content_rect = [wd.window_view frame];
	const float scale = ds->screen_get_max_scale();
	wd.size.width = content_rect.size.width * scale;
	wd.size.height = content_rect.size.height * scale;

	CALayer *layer = [wd.window_view layer];
	if (layer) {
		layer.contentsScale = scale;
	}

	ds->window_resize(window_id, wd.size.width, wd.size.height);

	if (!wd.rect_changed_callback.is_null()) {
		Variant size = Rect2i(ds->window_get_position(window_id), ds->window_get_size(window_id));
		Variant *sizep = &size;
		Variant ret;
		Callable::CallError ce;
		wd.rect_changed_callback.callp((const Variant **)&sizep, 1, ret, ce);
	}
}

- (void)windowDidMove:(NSNotification *)notification {
	DisplayServerMacOS *ds = (DisplayServerMacOS *)DisplayServer::get_singleton();
	if (!ds || !ds->has_window(window_id)) {
		return;
	}

	DisplayServerMacOS::WindowData &wd = ds->get_window(window_id);
	ds->release_pressed_events();

	if (!wd.rect_changed_callback.is_null()) {
		Variant size = Rect2i(ds->window_get_position(window_id), ds->window_get_size(window_id));
		Variant *sizep = &size;
		Variant ret;
		Callable::CallError ce;
		wd.rect_changed_callback.callp((const Variant **)&sizep, 1, ret, ce);
	}
}

- (void)windowDidBecomeKey:(NSNotification *)notification {
	DisplayServerMacOS *ds = (DisplayServerMacOS *)DisplayServer::get_singleton();
	if (!ds || !ds->has_window(window_id)) {
		return;
	}

	DisplayServerMacOS::WindowData &wd = ds->get_window(window_id);

	if (wd.window_button_view) {
		[(GodotButtonView *)wd.window_button_view displayButtons];
	}

	if (ds->mouse_get_mode() == DisplayServer::MOUSE_MODE_CAPTURED) {
		const NSRect content_rect = [wd.window_view frame];
		NSRect point_in_window_rect = NSMakeRect(content_rect.size.width / 2, content_rect.size.height / 2, 0, 0);
		NSPoint point_on_screen = [[wd.window_view window] convertRectToScreen:point_in_window_rect].origin;
		CGPoint mouse_warp_pos = { point_on_screen.x, CGDisplayBounds(CGMainDisplayID()).size.height - point_on_screen.y };
		CGWarpMouseCursorPosition(mouse_warp_pos);
	} else {
		ds->update_mouse_pos(wd, [wd.window_object mouseLocationOutsideOfEventStream]);
	}

	ds->set_last_focused_window(window_id);
	ds->send_window_event(wd, DisplayServerMacOS::WINDOW_EVENT_FOCUS_IN);
}

- (void)windowDidResignKey:(NSNotification *)notification {
	DisplayServerMacOS *ds = (DisplayServerMacOS *)DisplayServer::get_singleton();
	if (!ds || !ds->has_window(window_id)) {
		return;
	}

	DisplayServerMacOS::WindowData &wd = ds->get_window(window_id);

	if (wd.window_button_view) {
		[(GodotButtonView *)wd.window_button_view displayButtons];
	}

	ds->release_pressed_events();
	ds->send_window_event(wd, DisplayServerMacOS::WINDOW_EVENT_FOCUS_OUT);
}

- (void)windowDidMiniaturize:(NSNotification *)notification {
	DisplayServerMacOS *ds = (DisplayServerMacOS *)DisplayServer::get_singleton();
	if (!ds || !ds->has_window(window_id)) {
		return;
	}

	DisplayServerMacOS::WindowData &wd = ds->get_window(window_id);

	ds->release_pressed_events();
	ds->send_window_event(wd, DisplayServerMacOS::WINDOW_EVENT_FOCUS_OUT);
}

- (void)windowDidDeminiaturize:(NSNotification *)notification {
	DisplayServerMacOS *ds = (DisplayServerMacOS *)DisplayServer::get_singleton();
	if (!ds || !ds->has_window(window_id)) {
		return;
	}

	DisplayServerMacOS::WindowData &wd = ds->get_window(window_id);

	ds->set_last_focused_window(window_id);
	ds->send_window_event(wd, DisplayServerMacOS::WINDOW_EVENT_FOCUS_IN);
}

@end
