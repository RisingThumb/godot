#!/bin/sh

scons platform=linuxbsd target=template_debug arch=x86_64
scons platform=linuxbsd target=template_release arch=x86_64

mv ./bin/godot.linuxbsd.template_debug.double.x86_64 $HOME/.local/share/godot/export_templates/4.3.1.rc/linux_debug.double.x86_64
mv ./bin/godot.linuxbsd.template_release.double.x86_64 $HOME/.local/share/godot/export_templates/4.3.1.rc/linux_release.double.x86_64
