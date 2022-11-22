def can_build(env, platform):
    return not env["disable_3d"]


def configure(env):
    pass


def get_doc_classes():
    return [
        "EditorSceneFormatImporterBlend",
        "EditorSceneFormatImporterFBX",
        "EditorSceneFormatImporterGLTF",
        "GLTFAccessor",
        "GLTFAnimation",
        "GLTFBufferView",
        "GLTFCamera",
        "GLTFCollider",
        "GLTFDocument",
        "GLTFDocumentExtension",
        "GLTFDocumentExtensionConvertImporterMesh",
        "GLTFLight",
        "GLTFMesh",
        "GLTFNode",
        "GLTFPhysicsBody",
        "GLTFSkeleton",
        "GLTFSkin",
        "GLTFSpecGloss",
        "GLTFState",
        "GLTFTexture",
        "GLTFTextureSampler",
    ]


def get_doc_path():
    return "doc_classes"
