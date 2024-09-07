scons precision=double platform=windows target=template_debug arch=x86_64
scons precision=double platform=windows target=template_release arch=x86_64

mv ./bin/godot.windows.template_debug.x86_64.exe $HOME/.local/share/godot/export_templates/4.3.beta/windows_64_debug.exe
#mv ./bin/godot.windows.template_debug.x86_64.console.exe $HOME/.local/share/godot/export_templates/4.3.beta/windows_64_debug.exe
mv ./bin/godot.windows.template_release.x86_64.exe $HOME/.local/share/godot/export_templates/4.3.beta/windows_64_release.exe
#mv ./bin/godot.windows.template_release.x86_64.console.exe $HOME/.local/share/godot/export_templates/4.3.beta/windows_64_release.exe
