#!/bin/sh

scons platform=linuxbsd target=template_debug arch=x86_64
scons platform=linuxbsd target=template_release arch=x86_64

mv ./bin/godot.linuxbsd.template_debug.x86_64 $HOME/.local/share/godot/export_templates/4.2.1.stable/linux_debug.x86_64
mv ./bin/godot.linuxbsd.template_release.x86_64 $HOME/.local/share/godot/export_templates/4.2.1.stable/linux_release.x86_64
