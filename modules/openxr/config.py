def can_build(env, platform):
    print("platform " + platform)
    if platform == "android" or platform == "linuxbsd" or platform == "windows":
        return env["openxr"]
    else:
        # not supported on these platforms
        return False


def configure(env):
    pass


def get_doc_classes():
    return [
        "OpenXRInterface",
        "OpenXRAction",
        "OpenXRActionSet",
        "OpenXRActionMap",
        "OpenXRInteractionProfile",
        "OpenXRIPBinding",
    ]


def get_doc_path():
    return "doc_classes"
