def can_build(env, platform):
    return True


def configure(env):
    pass


def get_doc_path():
    return "doc_classes"


def get_doc_classes():
    return [
        "VGColor",
        "VGGradient",
        "VGLinearGradient",
        "VGMeshRenderer",
        "VGPaint",
        "VGPath",
        "VGRadialGradient",
        "VGRenderer",
        "EditorSceneImporterSVG",
    ]
