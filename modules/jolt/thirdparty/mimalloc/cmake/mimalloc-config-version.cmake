set(mi_version_major 2)
set(mi_version_minor 1)
set(mi_version_patch 2)
set(mi_version ${mi_version_major}.${mi_version_minor})

set(PACKAGE_VERSION ${mi_version})
if(PACKAGE_FIND_VERSION_MAJOR)
    if("${PACKAGE_FIND_VERSION_MAJOR}" EQUAL "${mi_version_major}")
        if ("${PACKAGE_FIND_VERSION_MINOR}" EQUAL "${mi_version_minor}")
            set(PACKAGE_VERSION_EXACT TRUE)
        elseif("${PACKAGE_FIND_VERSION_MINOR}" LESS "${mi_version_minor}")
            set(PACKAGE_VERSION_COMPATIBLE TRUE)
        else()
            set(PACKAGE_VERSION_UNSUITABLE TRUE)
        endif()
    else()
        set(PACKAGE_VERSION_UNSUITABLE TRUE)
    endif()
endif()
