import subprocess
import os


def can_build(env, platform):
    if platform in ("ios", "web", "android"):
        return False

    try:
        subprocess.check_output(["rustup", "--version"], stderr=subprocess.STDOUT)
        subprocess.check_output(["cargo", "--version"], stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print("Cargo or Rustup not found. mvsqlite build skipped.")
        if platform == "windows":
            print("Use `scoop install rustup-gnu` and `rustup target add x86_64-pc-windows-gnu` to install rust.")
        return False

    if platform == "windows":
        use_mingw = env["use_mingw"]
        # Check if MSVC target is installed
        try:
            targets_output = subprocess.check_output(["rustup", "target", "list", "--installed"], stderr=subprocess.STDOUT)
            installed_targets = targets_output.decode().splitlines()
            if "x86_64-pc-windows-msvc" in installed_targets and not use_mingw:
                print("MSVC target is installed. mvsqlite build skipped.")
                return False
        except subprocess.CalledProcessError:
            print("Failed to list installed rustup targets.")
            return False

    return True


def configure(env):
    pass


def get_doc_classes():
    return [
        "MVSQLite",
        "MVSQLiteQuery",
    ]


def get_doc_path():
    return "doc_classes"
