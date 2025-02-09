# Avoid multiple calls to find_package to append duplicated properties to the targets
include_guard()########### VARIABLES #######################################################################
#############################################################################################
set(vorbis_FRAMEWORKS_FOUND_DEBUG "") # Will be filled later
conan_find_apple_frameworks(vorbis_FRAMEWORKS_FOUND_DEBUG "${vorbis_FRAMEWORKS_DEBUG}" "${vorbis_FRAMEWORK_DIRS_DEBUG}")

set(vorbis_LIBRARIES_TARGETS "") # Will be filled later


######## Create an interface target to contain all the dependencies (frameworks, system and conan deps)
if(NOT TARGET vorbis_DEPS_TARGET)
    add_library(vorbis_DEPS_TARGET INTERFACE IMPORTED)
endif()

set_property(TARGET vorbis_DEPS_TARGET
             APPEND PROPERTY INTERFACE_LINK_LIBRARIES
             $<$<CONFIG:Debug>:${vorbis_FRAMEWORKS_FOUND_DEBUG}>
             $<$<CONFIG:Debug>:${vorbis_SYSTEM_LIBS_DEBUG}>
             $<$<CONFIG:Debug>:Ogg::ogg;Vorbis::vorbis;Vorbis::vorbisenc;Vorbis::vorbisfile>)

####### Find the libraries declared in cpp_info.libs, create an IMPORTED target for each one and link the
####### vorbis_DEPS_TARGET to all of them
conan_package_library_targets("${vorbis_LIBS_DEBUG}"    # libraries
                              "${vorbis_LIB_DIRS_DEBUG}" # package_libdir
                              "${vorbis_BIN_DIRS_DEBUG}" # package_bindir
                              "${vorbis_LIBRARY_TYPE_DEBUG}"
                              "${vorbis_IS_HOST_WINDOWS_DEBUG}"
                              vorbis_DEPS_TARGET
                              vorbis_LIBRARIES_TARGETS  # out_libraries_targets
                              "_DEBUG"
                              "vorbis"    # package_name
                              "${vorbis_NO_SONAME_MODE_DEBUG}")  # soname

# FIXME: What is the result of this for multi-config? All configs adding themselves to path?
set(CMAKE_MODULE_PATH ${vorbis_BUILD_DIRS_DEBUG} ${CMAKE_MODULE_PATH})

########## COMPONENTS TARGET PROPERTIES Debug ########################################

    ########## COMPONENT vorbis::vorbisfile-alias #############

        set(vorbis_vorbis_vorbisfile-alias_FRAMEWORKS_FOUND_DEBUG "")
        conan_find_apple_frameworks(vorbis_vorbis_vorbisfile-alias_FRAMEWORKS_FOUND_DEBUG "${vorbis_vorbis_vorbisfile-alias_FRAMEWORKS_DEBUG}" "${vorbis_vorbis_vorbisfile-alias_FRAMEWORK_DIRS_DEBUG}")

        set(vorbis_vorbis_vorbisfile-alias_LIBRARIES_TARGETS "")

        ######## Create an interface target to contain all the dependencies (frameworks, system and conan deps)
        if(NOT TARGET vorbis_vorbis_vorbisfile-alias_DEPS_TARGET)
            add_library(vorbis_vorbis_vorbisfile-alias_DEPS_TARGET INTERFACE IMPORTED)
        endif()

        set_property(TARGET vorbis_vorbis_vorbisfile-alias_DEPS_TARGET
                     APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                     $<$<CONFIG:Debug>:${vorbis_vorbis_vorbisfile-alias_FRAMEWORKS_FOUND_DEBUG}>
                     $<$<CONFIG:Debug>:${vorbis_vorbis_vorbisfile-alias_SYSTEM_LIBS_DEBUG}>
                     $<$<CONFIG:Debug>:${vorbis_vorbis_vorbisfile-alias_DEPENDENCIES_DEBUG}>
                     )

        ####### Find the libraries declared in cpp_info.component["xxx"].libs,
        ####### create an IMPORTED target for each one and link the 'vorbis_vorbis_vorbisfile-alias_DEPS_TARGET' to all of them
        conan_package_library_targets("${vorbis_vorbis_vorbisfile-alias_LIBS_DEBUG}"
                              "${vorbis_vorbis_vorbisfile-alias_LIB_DIRS_DEBUG}"
                              "${vorbis_vorbis_vorbisfile-alias_BIN_DIRS_DEBUG}" # package_bindir
                              "${vorbis_vorbis_vorbisfile-alias_LIBRARY_TYPE_DEBUG}"
                              "${vorbis_vorbis_vorbisfile-alias_IS_HOST_WINDOWS_DEBUG}"
                              vorbis_vorbis_vorbisfile-alias_DEPS_TARGET
                              vorbis_vorbis_vorbisfile-alias_LIBRARIES_TARGETS
                              "_DEBUG"
                              "vorbis_vorbis_vorbisfile-alias"
                              "${vorbis_vorbis_vorbisfile-alias_NO_SONAME_MODE_DEBUG}")


        ########## TARGET PROPERTIES #####################################
        set_property(TARGET vorbis::vorbisfile-alias
                     APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                     $<$<CONFIG:Debug>:${vorbis_vorbis_vorbisfile-alias_OBJECTS_DEBUG}>
                     $<$<CONFIG:Debug>:${vorbis_vorbis_vorbisfile-alias_LIBRARIES_TARGETS}>
                     )

        if("${vorbis_vorbis_vorbisfile-alias_LIBS_DEBUG}" STREQUAL "")
            # If the component is not declaring any "cpp_info.components['foo'].libs" the system, frameworks etc are not
            # linked to the imported targets and we need to do it to the global target
            set_property(TARGET vorbis::vorbisfile-alias
                         APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                         vorbis_vorbis_vorbisfile-alias_DEPS_TARGET)
        endif()

        set_property(TARGET vorbis::vorbisfile-alias APPEND PROPERTY INTERFACE_LINK_OPTIONS
                     $<$<CONFIG:Debug>:${vorbis_vorbis_vorbisfile-alias_LINKER_FLAGS_DEBUG}>)
        set_property(TARGET vorbis::vorbisfile-alias APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                     $<$<CONFIG:Debug>:${vorbis_vorbis_vorbisfile-alias_INCLUDE_DIRS_DEBUG}>)
        set_property(TARGET vorbis::vorbisfile-alias APPEND PROPERTY INTERFACE_LINK_DIRECTORIES
                     $<$<CONFIG:Debug>:${vorbis_vorbis_vorbisfile-alias_LIB_DIRS_DEBUG}>)
        set_property(TARGET vorbis::vorbisfile-alias APPEND PROPERTY INTERFACE_COMPILE_DEFINITIONS
                     $<$<CONFIG:Debug>:${vorbis_vorbis_vorbisfile-alias_COMPILE_DEFINITIONS_DEBUG}>)
        set_property(TARGET vorbis::vorbisfile-alias APPEND PROPERTY INTERFACE_COMPILE_OPTIONS
                     $<$<CONFIG:Debug>:${vorbis_vorbis_vorbisfile-alias_COMPILE_OPTIONS_DEBUG}>)

    ########## COMPONENT vorbis::vorbisenc-alias #############

        set(vorbis_vorbis_vorbisenc-alias_FRAMEWORKS_FOUND_DEBUG "")
        conan_find_apple_frameworks(vorbis_vorbis_vorbisenc-alias_FRAMEWORKS_FOUND_DEBUG "${vorbis_vorbis_vorbisenc-alias_FRAMEWORKS_DEBUG}" "${vorbis_vorbis_vorbisenc-alias_FRAMEWORK_DIRS_DEBUG}")

        set(vorbis_vorbis_vorbisenc-alias_LIBRARIES_TARGETS "")

        ######## Create an interface target to contain all the dependencies (frameworks, system and conan deps)
        if(NOT TARGET vorbis_vorbis_vorbisenc-alias_DEPS_TARGET)
            add_library(vorbis_vorbis_vorbisenc-alias_DEPS_TARGET INTERFACE IMPORTED)
        endif()

        set_property(TARGET vorbis_vorbis_vorbisenc-alias_DEPS_TARGET
                     APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                     $<$<CONFIG:Debug>:${vorbis_vorbis_vorbisenc-alias_FRAMEWORKS_FOUND_DEBUG}>
                     $<$<CONFIG:Debug>:${vorbis_vorbis_vorbisenc-alias_SYSTEM_LIBS_DEBUG}>
                     $<$<CONFIG:Debug>:${vorbis_vorbis_vorbisenc-alias_DEPENDENCIES_DEBUG}>
                     )

        ####### Find the libraries declared in cpp_info.component["xxx"].libs,
        ####### create an IMPORTED target for each one and link the 'vorbis_vorbis_vorbisenc-alias_DEPS_TARGET' to all of them
        conan_package_library_targets("${vorbis_vorbis_vorbisenc-alias_LIBS_DEBUG}"
                              "${vorbis_vorbis_vorbisenc-alias_LIB_DIRS_DEBUG}"
                              "${vorbis_vorbis_vorbisenc-alias_BIN_DIRS_DEBUG}" # package_bindir
                              "${vorbis_vorbis_vorbisenc-alias_LIBRARY_TYPE_DEBUG}"
                              "${vorbis_vorbis_vorbisenc-alias_IS_HOST_WINDOWS_DEBUG}"
                              vorbis_vorbis_vorbisenc-alias_DEPS_TARGET
                              vorbis_vorbis_vorbisenc-alias_LIBRARIES_TARGETS
                              "_DEBUG"
                              "vorbis_vorbis_vorbisenc-alias"
                              "${vorbis_vorbis_vorbisenc-alias_NO_SONAME_MODE_DEBUG}")


        ########## TARGET PROPERTIES #####################################
        set_property(TARGET vorbis::vorbisenc-alias
                     APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                     $<$<CONFIG:Debug>:${vorbis_vorbis_vorbisenc-alias_OBJECTS_DEBUG}>
                     $<$<CONFIG:Debug>:${vorbis_vorbis_vorbisenc-alias_LIBRARIES_TARGETS}>
                     )

        if("${vorbis_vorbis_vorbisenc-alias_LIBS_DEBUG}" STREQUAL "")
            # If the component is not declaring any "cpp_info.components['foo'].libs" the system, frameworks etc are not
            # linked to the imported targets and we need to do it to the global target
            set_property(TARGET vorbis::vorbisenc-alias
                         APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                         vorbis_vorbis_vorbisenc-alias_DEPS_TARGET)
        endif()

        set_property(TARGET vorbis::vorbisenc-alias APPEND PROPERTY INTERFACE_LINK_OPTIONS
                     $<$<CONFIG:Debug>:${vorbis_vorbis_vorbisenc-alias_LINKER_FLAGS_DEBUG}>)
        set_property(TARGET vorbis::vorbisenc-alias APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                     $<$<CONFIG:Debug>:${vorbis_vorbis_vorbisenc-alias_INCLUDE_DIRS_DEBUG}>)
        set_property(TARGET vorbis::vorbisenc-alias APPEND PROPERTY INTERFACE_LINK_DIRECTORIES
                     $<$<CONFIG:Debug>:${vorbis_vorbis_vorbisenc-alias_LIB_DIRS_DEBUG}>)
        set_property(TARGET vorbis::vorbisenc-alias APPEND PROPERTY INTERFACE_COMPILE_DEFINITIONS
                     $<$<CONFIG:Debug>:${vorbis_vorbis_vorbisenc-alias_COMPILE_DEFINITIONS_DEBUG}>)
        set_property(TARGET vorbis::vorbisenc-alias APPEND PROPERTY INTERFACE_COMPILE_OPTIONS
                     $<$<CONFIG:Debug>:${vorbis_vorbis_vorbisenc-alias_COMPILE_OPTIONS_DEBUG}>)

    ########## COMPONENT Vorbis::vorbisfile #############

        set(vorbis_Vorbis_vorbisfile_FRAMEWORKS_FOUND_DEBUG "")
        conan_find_apple_frameworks(vorbis_Vorbis_vorbisfile_FRAMEWORKS_FOUND_DEBUG "${vorbis_Vorbis_vorbisfile_FRAMEWORKS_DEBUG}" "${vorbis_Vorbis_vorbisfile_FRAMEWORK_DIRS_DEBUG}")

        set(vorbis_Vorbis_vorbisfile_LIBRARIES_TARGETS "")

        ######## Create an interface target to contain all the dependencies (frameworks, system and conan deps)
        if(NOT TARGET vorbis_Vorbis_vorbisfile_DEPS_TARGET)
            add_library(vorbis_Vorbis_vorbisfile_DEPS_TARGET INTERFACE IMPORTED)
        endif()

        set_property(TARGET vorbis_Vorbis_vorbisfile_DEPS_TARGET
                     APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbisfile_FRAMEWORKS_FOUND_DEBUG}>
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbisfile_SYSTEM_LIBS_DEBUG}>
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbisfile_DEPENDENCIES_DEBUG}>
                     )

        ####### Find the libraries declared in cpp_info.component["xxx"].libs,
        ####### create an IMPORTED target for each one and link the 'vorbis_Vorbis_vorbisfile_DEPS_TARGET' to all of them
        conan_package_library_targets("${vorbis_Vorbis_vorbisfile_LIBS_DEBUG}"
                              "${vorbis_Vorbis_vorbisfile_LIB_DIRS_DEBUG}"
                              "${vorbis_Vorbis_vorbisfile_BIN_DIRS_DEBUG}" # package_bindir
                              "${vorbis_Vorbis_vorbisfile_LIBRARY_TYPE_DEBUG}"
                              "${vorbis_Vorbis_vorbisfile_IS_HOST_WINDOWS_DEBUG}"
                              vorbis_Vorbis_vorbisfile_DEPS_TARGET
                              vorbis_Vorbis_vorbisfile_LIBRARIES_TARGETS
                              "_DEBUG"
                              "vorbis_Vorbis_vorbisfile"
                              "${vorbis_Vorbis_vorbisfile_NO_SONAME_MODE_DEBUG}")


        ########## TARGET PROPERTIES #####################################
        set_property(TARGET Vorbis::vorbisfile
                     APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbisfile_OBJECTS_DEBUG}>
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbisfile_LIBRARIES_TARGETS}>
                     )

        if("${vorbis_Vorbis_vorbisfile_LIBS_DEBUG}" STREQUAL "")
            # If the component is not declaring any "cpp_info.components['foo'].libs" the system, frameworks etc are not
            # linked to the imported targets and we need to do it to the global target
            set_property(TARGET Vorbis::vorbisfile
                         APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                         vorbis_Vorbis_vorbisfile_DEPS_TARGET)
        endif()

        set_property(TARGET Vorbis::vorbisfile APPEND PROPERTY INTERFACE_LINK_OPTIONS
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbisfile_LINKER_FLAGS_DEBUG}>)
        set_property(TARGET Vorbis::vorbisfile APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbisfile_INCLUDE_DIRS_DEBUG}>)
        set_property(TARGET Vorbis::vorbisfile APPEND PROPERTY INTERFACE_LINK_DIRECTORIES
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbisfile_LIB_DIRS_DEBUG}>)
        set_property(TARGET Vorbis::vorbisfile APPEND PROPERTY INTERFACE_COMPILE_DEFINITIONS
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbisfile_COMPILE_DEFINITIONS_DEBUG}>)
        set_property(TARGET Vorbis::vorbisfile APPEND PROPERTY INTERFACE_COMPILE_OPTIONS
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbisfile_COMPILE_OPTIONS_DEBUG}>)

    ########## COMPONENT Vorbis::vorbisenc #############

        set(vorbis_Vorbis_vorbisenc_FRAMEWORKS_FOUND_DEBUG "")
        conan_find_apple_frameworks(vorbis_Vorbis_vorbisenc_FRAMEWORKS_FOUND_DEBUG "${vorbis_Vorbis_vorbisenc_FRAMEWORKS_DEBUG}" "${vorbis_Vorbis_vorbisenc_FRAMEWORK_DIRS_DEBUG}")

        set(vorbis_Vorbis_vorbisenc_LIBRARIES_TARGETS "")

        ######## Create an interface target to contain all the dependencies (frameworks, system and conan deps)
        if(NOT TARGET vorbis_Vorbis_vorbisenc_DEPS_TARGET)
            add_library(vorbis_Vorbis_vorbisenc_DEPS_TARGET INTERFACE IMPORTED)
        endif()

        set_property(TARGET vorbis_Vorbis_vorbisenc_DEPS_TARGET
                     APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbisenc_FRAMEWORKS_FOUND_DEBUG}>
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbisenc_SYSTEM_LIBS_DEBUG}>
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbisenc_DEPENDENCIES_DEBUG}>
                     )

        ####### Find the libraries declared in cpp_info.component["xxx"].libs,
        ####### create an IMPORTED target for each one and link the 'vorbis_Vorbis_vorbisenc_DEPS_TARGET' to all of them
        conan_package_library_targets("${vorbis_Vorbis_vorbisenc_LIBS_DEBUG}"
                              "${vorbis_Vorbis_vorbisenc_LIB_DIRS_DEBUG}"
                              "${vorbis_Vorbis_vorbisenc_BIN_DIRS_DEBUG}" # package_bindir
                              "${vorbis_Vorbis_vorbisenc_LIBRARY_TYPE_DEBUG}"
                              "${vorbis_Vorbis_vorbisenc_IS_HOST_WINDOWS_DEBUG}"
                              vorbis_Vorbis_vorbisenc_DEPS_TARGET
                              vorbis_Vorbis_vorbisenc_LIBRARIES_TARGETS
                              "_DEBUG"
                              "vorbis_Vorbis_vorbisenc"
                              "${vorbis_Vorbis_vorbisenc_NO_SONAME_MODE_DEBUG}")


        ########## TARGET PROPERTIES #####################################
        set_property(TARGET Vorbis::vorbisenc
                     APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbisenc_OBJECTS_DEBUG}>
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbisenc_LIBRARIES_TARGETS}>
                     )

        if("${vorbis_Vorbis_vorbisenc_LIBS_DEBUG}" STREQUAL "")
            # If the component is not declaring any "cpp_info.components['foo'].libs" the system, frameworks etc are not
            # linked to the imported targets and we need to do it to the global target
            set_property(TARGET Vorbis::vorbisenc
                         APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                         vorbis_Vorbis_vorbisenc_DEPS_TARGET)
        endif()

        set_property(TARGET Vorbis::vorbisenc APPEND PROPERTY INTERFACE_LINK_OPTIONS
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbisenc_LINKER_FLAGS_DEBUG}>)
        set_property(TARGET Vorbis::vorbisenc APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbisenc_INCLUDE_DIRS_DEBUG}>)
        set_property(TARGET Vorbis::vorbisenc APPEND PROPERTY INTERFACE_LINK_DIRECTORIES
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbisenc_LIB_DIRS_DEBUG}>)
        set_property(TARGET Vorbis::vorbisenc APPEND PROPERTY INTERFACE_COMPILE_DEFINITIONS
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbisenc_COMPILE_DEFINITIONS_DEBUG}>)
        set_property(TARGET Vorbis::vorbisenc APPEND PROPERTY INTERFACE_COMPILE_OPTIONS
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbisenc_COMPILE_OPTIONS_DEBUG}>)

    ########## COMPONENT Vorbis::vorbis #############

        set(vorbis_Vorbis_vorbis_FRAMEWORKS_FOUND_DEBUG "")
        conan_find_apple_frameworks(vorbis_Vorbis_vorbis_FRAMEWORKS_FOUND_DEBUG "${vorbis_Vorbis_vorbis_FRAMEWORKS_DEBUG}" "${vorbis_Vorbis_vorbis_FRAMEWORK_DIRS_DEBUG}")

        set(vorbis_Vorbis_vorbis_LIBRARIES_TARGETS "")

        ######## Create an interface target to contain all the dependencies (frameworks, system and conan deps)
        if(NOT TARGET vorbis_Vorbis_vorbis_DEPS_TARGET)
            add_library(vorbis_Vorbis_vorbis_DEPS_TARGET INTERFACE IMPORTED)
        endif()

        set_property(TARGET vorbis_Vorbis_vorbis_DEPS_TARGET
                     APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbis_FRAMEWORKS_FOUND_DEBUG}>
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbis_SYSTEM_LIBS_DEBUG}>
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbis_DEPENDENCIES_DEBUG}>
                     )

        ####### Find the libraries declared in cpp_info.component["xxx"].libs,
        ####### create an IMPORTED target for each one and link the 'vorbis_Vorbis_vorbis_DEPS_TARGET' to all of them
        conan_package_library_targets("${vorbis_Vorbis_vorbis_LIBS_DEBUG}"
                              "${vorbis_Vorbis_vorbis_LIB_DIRS_DEBUG}"
                              "${vorbis_Vorbis_vorbis_BIN_DIRS_DEBUG}" # package_bindir
                              "${vorbis_Vorbis_vorbis_LIBRARY_TYPE_DEBUG}"
                              "${vorbis_Vorbis_vorbis_IS_HOST_WINDOWS_DEBUG}"
                              vorbis_Vorbis_vorbis_DEPS_TARGET
                              vorbis_Vorbis_vorbis_LIBRARIES_TARGETS
                              "_DEBUG"
                              "vorbis_Vorbis_vorbis"
                              "${vorbis_Vorbis_vorbis_NO_SONAME_MODE_DEBUG}")


        ########## TARGET PROPERTIES #####################################
        set_property(TARGET Vorbis::vorbis
                     APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbis_OBJECTS_DEBUG}>
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbis_LIBRARIES_TARGETS}>
                     )

        if("${vorbis_Vorbis_vorbis_LIBS_DEBUG}" STREQUAL "")
            # If the component is not declaring any "cpp_info.components['foo'].libs" the system, frameworks etc are not
            # linked to the imported targets and we need to do it to the global target
            set_property(TARGET Vorbis::vorbis
                         APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                         vorbis_Vorbis_vorbis_DEPS_TARGET)
        endif()

        set_property(TARGET Vorbis::vorbis APPEND PROPERTY INTERFACE_LINK_OPTIONS
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbis_LINKER_FLAGS_DEBUG}>)
        set_property(TARGET Vorbis::vorbis APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbis_INCLUDE_DIRS_DEBUG}>)
        set_property(TARGET Vorbis::vorbis APPEND PROPERTY INTERFACE_LINK_DIRECTORIES
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbis_LIB_DIRS_DEBUG}>)
        set_property(TARGET Vorbis::vorbis APPEND PROPERTY INTERFACE_COMPILE_DEFINITIONS
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbis_COMPILE_DEFINITIONS_DEBUG}>)
        set_property(TARGET Vorbis::vorbis APPEND PROPERTY INTERFACE_COMPILE_OPTIONS
                     $<$<CONFIG:Debug>:${vorbis_Vorbis_vorbis_COMPILE_OPTIONS_DEBUG}>)

    ########## AGGREGATED GLOBAL TARGET WITH THE COMPONENTS #####################
    set_property(TARGET vorbis::vorbis APPEND PROPERTY INTERFACE_LINK_LIBRARIES vorbis::vorbisfile-alias)
    set_property(TARGET vorbis::vorbis APPEND PROPERTY INTERFACE_LINK_LIBRARIES vorbis::vorbisenc-alias)
    set_property(TARGET vorbis::vorbis APPEND PROPERTY INTERFACE_LINK_LIBRARIES Vorbis::vorbisfile)
    set_property(TARGET vorbis::vorbis APPEND PROPERTY INTERFACE_LINK_LIBRARIES Vorbis::vorbisenc)
    set_property(TARGET vorbis::vorbis APPEND PROPERTY INTERFACE_LINK_LIBRARIES Vorbis::vorbis)

########## For the modules (FindXXX)
set(vorbis_LIBRARIES_DEBUG vorbis::vorbis)
