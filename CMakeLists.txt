# ######################################################################### #
# Georgiev Lab (c) 2015-2016                                                #
# ######################################################################### #
# Department of Cybernetics                                                 #
# Faculty of Applied Sciences                                               #
# University of West Bohemia in Pilsen                                      #
# ######################################################################### #
#                                                                           #
# This file is part of CeCe.                                                #
#                                                                           #
# CeCe is free software: you can redistribute it and/or modify              #
# it under the terms of the GNU General Public License as published by      #
# the Free Software Foundation, either version 3 of the License, or         #
# (at your option) any later version.                                       #
#                                                                           #
# CeCe is distributed in the hope that it will be useful,                   #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU General Public License for more details.                              #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with CeCe.  If not, see <http://www.gnu.org/licenses/>.             #
#                                                                           #
# ######################################################################### #

cmake_minimum_required(VERSION 3.1)

# ######################################################################### #

project(cece-plugins
    VERSION 0.6.3
    LANGUAGES CXX
)

# ######################################################################### #

# In case of standalone build, checkout CeCe-core repository into subdirectory
if (NOT TARGET cece)
    if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/CeCe-core/CMakeLists.txt")
        find_package(Git REQUIRED)
        execute_process(
            COMMAND ${GIT_EXECUTABLE} submodule update --init CeCe-core
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
    endif()

    add_subdirectory(CeCe-core)

    find_program(CCACHE_FOUND ccache)
    if (CCACHE_FOUND)
        set_property(DIRECTORY CeCe-core PROPERTY RULE_LAUNCH_COMPILE ccache)
        set_property(DIRECTORY CeCe-core PROPERTY RULE_LAUNCH_LINK ccache)
    endif()
endif ()

# ######################################################################### #

option(CECE_PLUGINS_DEV "Enable special mode for developers that symlink/copy plugins into a single directory" Off)

# A list of custom plugins
set(CECE_PLUGINS_CUSTOM "" CACHE String "A list of additional plugins")

# ######################################################################### #

if (CECE_TESTS_BUILD)
    enable_testing()
endif()

# ######################################################################### #

# Standard plugins
set(CECE_PLUGINS_STANDARD
    agglutination
    background
    cell
    cell-python
    chemostat
    cytometry
    diffusion
    diffusion-python
    diffusion-streamlines
    diffusion-streamlines-python
    fluorescence
    object-generator
    object-streamlines-generator
    obstacle
    obstacles-svg
    picture
    python
    stochastic-reactions
    streamlines
#    streamlines-channel
    streamlines-python
    xml
)

# ######################################################################### #

set(CECE_PLUGINS ${CECE_PLUGINS_STANDARD} ${CECE_PLUGINS_CUSTOM})

# ######################################################################### #

##
## Get plugin project name.
##
## param: PLUGIN_NAME  Plugin
## param: PROJECT_NAME Name of output variable.
##
function(plugin_project_name PLUGIN_NAME PROJECT_NAME)
    set(${PROJECT_NAME} "cece-${PLUGIN_NAME}" PARENT_SCOPE)
endfunction ()

# ######################################################################### #

##
## Build plugin helper function
##
## param: NAME             Plugin name.
## param: SOURCES          List of sources.
## param: LIBRARIES        A list of required libraries.
## param: PLUGINS_REQUIRED A list of required plugins.
##
## Function uses PLUGINS_BUILTIN_${NAME} variable.
##
function(build_plugin NAME)

    # Get project name
    plugin_project_name(${NAME} FULLNAME)

    include(CMakeParseArguments)
    cmake_parse_arguments(ARG "EXPORT_HEADER" "" "SOURCES;LIBRARIES;DEPENDENCIES;PLUGINS_REQUIRED" ${ARGN})

    # Modify list of plugins (add prefixes)
    set(PLUGINS_DEPENDENCIES "")

    foreach (PLUGIN ${ARG_PLUGINS_REQUIRED})
        plugin_project_name(${PLUGIN} PLUGIN_NAME)
        list(APPEND PLUGINS_DEPENDENCIES "${PLUGIN_NAME}")
    endforeach ()

    # Shared module
    add_library(${FULLNAME} SHARED ${ARG_SOURCES})

    # Setup dependencies
    if (${ARG_DEPENDENCIES})
        add_dependencies(${FULLNAME}
            ${ARG_DEPENDENCIES}
        )
    endif ()

    # Link libraries
    target_link_libraries(${FULLNAME}
        ${PLUGINS_DEPENDENCIES}
        ${ARG_LIBRARIES}
        cece
    )

    set_target_properties(${FULLNAME} PROPERTIES
        INSTALL_RPATH "$ORIGIN"
        MACOSX_RPATH On
        CXX_STANDARD 11
        CXX_EXTENSIONS Off
        CXX_STANDARD_REQUIRED On
        WINDOWS_EXPORT_ALL_SYMBOLS On
    )

    # Generate export header
    if (ARG_EXPORT_HEADER)
        include(GenerateExportHeader)
        target_include_directories(${FULLNAME}
            PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
        )
        generate_export_header(${FULLNAME}
            EXPORT_FILE_NAME "export.hpp"
            BASE_NAME "CECE_PLUGIN_${NAME}"
        )
    endif ()

    if (CECE_PLUGINS_DEV)
        if (UNIX)
            add_custom_command(
                TARGET ${FULLNAME} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E
                    create_symlink $<TARGET_FILE:${FULLNAME}> ${CMAKE_CURRENT_BINARY_DIR}/../$<TARGET_FILE_NAME:${FULLNAME}>
            )
        elseif (WIN32)
            add_custom_command(
                TARGET ${FULLNAME} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E
                    copy $<TARGET_FILE:${FULLNAME}> ${CMAKE_CURRENT_BINARY_DIR}/../$<TARGET_FILE_NAME:${FULLNAME}>
            )
        endif ()
    endif ()

    install(TARGETS ${FULLNAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
    )
endfunction()

# ######################################################################### #

##
## Build tests
##
## param: PROJECT_NAME Project name.
##
function(build_test PROJECT_NAME)

    if (NOT CECE_TESTS_BUILD)
        message(WARNING "Trying to build test '${PROJECT_NAME}' without CECE_TESTS_BUILD")
    endif ()

    # Parse arguments
    include(CMakeParseArguments)
    cmake_parse_arguments(ARG "" "" "SOURCES;LIBRARIES;DEPENDENCIES;PLUGINS_REQUIRED" ${ARGN})

    # Modify list of plugins (add prefixes)
    set(PLUGINS_DEPENDENCIES "")

    foreach (PLUGIN ${ARG_PLUGINS_REQUIRED})
        plugin_project_name(${PLUGIN} PLUGIN_NAME)
        list(APPEND PLUGINS_DEPENDENCIES "${PLUGIN_NAME}")
    endforeach ()

    # Test name
    set(TEST_NAME "${PROJECT_NAME}_test")

    # Create executable
    add_executable(${TEST_NAME}
        ${ARG_SOURCES}
    )

    set_target_properties(${TEST_NAME} PROPERTIES
        CXX_STANDARD 11
        CXX_EXTENSIONS Off
        CXX_STANDARD_REQUIRED On
    )

    # Link libraries
    target_link_libraries(${TEST_NAME}
        ${ARG_LIBRARIES}
        ${PLUGINS_DEPENDENCIES}
        gtest_main
    )

    # Add test
    add_test(${PROJECT_NAME} ${TEST_NAME})
endfunction ()

# ######################################################################### #

set(CECE_PLUGINS_BUILD "")

# Build plugins option
foreach(PLUGIN_NAME ${CECE_PLUGINS})
    if (IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${PLUGIN_NAME}")
        option(CECE_PLUGINS_BUILD_${PLUGIN_NAME} "Build plugin '${PLUGIN_NAME}'" On)

        if (CECE_PLUGINS_BUILD_${PLUGIN_NAME})
            list(APPEND CECE_PLUGINS_BUILD ${PLUGIN_NAME})
            add_subdirectory(${PLUGIN_NAME})
            if (CECE_PLUGINS_BUILD_${PLUGIN_NAME}_DISABLE)
                # Disable plugin
                set(CECE_PLUGINS_BUILD_${PLUGIN_NAME} Off CACHE BOOL "" FORCE)
            endif ()
        endif ()
    endif()
endforeach()

mark_as_advanced(CECE_PLUGINS_BUILD)

# ######################################################################### #

add_custom_target(${PROJECT_NAME}
    DEPENDS ${CECE_PLUGINS_BUILD}
)

set_target_properties(${PROJECT_NAME} PROPERTIES
    OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
)

# ######################################################################### #

if (${CMAKE_SOURCE_DIR} MATCHES ${CMAKE_CURRENT_SOURCE_DIR})
    install(FILES LICENSE.txt Copyright.txt
        DESTINATION .
    )

    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "CeCe standard plugins")
    set(CPACK_PACKAGE_VENDOR "Georgiev Lab")
    set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
    set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
    set(CPACK_PACKAGE_NAME CeCe-plugins)
    set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")
    set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
    set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
    set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
    set(CPACK_PACKAGE_CONTACT "Jiří Fatka <fatkaj@ntis.zcu.cz>")

    find_package(Git)
    if (GIT_FOUND)
        execute_process(
            COMMAND ${GIT_EXECUTABLE} describe --exact-match --tags
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_TAG
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )

        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_COMMIT
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )

        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_BRANCH
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    endif ()

    message(STATUS "TAG: ${GIT_TAG}")
    message(STATUS "COMMIT: ${GIT_COMMIT}")
    message(STATUS "BRANCH: ${GIT_BRANCH}")

    if (NOT GIT_TAG AND GIT_COMMIT AND GIT_BRANCH)
        set(SYSTEM_NAME ${CMAKE_SYSTEM_NAME})
        if (${SYSTEM_NAME} MATCHES "Windows")
            if (CMAKE_CL_64)
                set(SYSTEM_NAME "Win64")
            else ()
                set(SYSTEM_NAME "Win32")
            endif ()
        endif ()

        string(SUBSTRING "${GIT_COMMIT}" 0 7 COMMIT)
        set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${SYSTEM_NAME}-${GIT_BRANCH}-${COMMIT}")
    endif ()

    # Configure file for travis-ci -> bintray deployment
    if (UNIX AND NOT APPLE)
        set(BINTRAY_PLATFORM "linux-x64")
    elseif (APPLE)
        set(BINTRAY_PLATFORM "macOS")
    endif ()

    if (NOT GIT_TAG AND GIT_COMMIT AND GIT_BRANCH)
        set(BINTRAY_VERSION ${PROJECT_VERSION}-${COMMIT})
    else ()
        set(BINTRAY_VERSION ${PROJECT_VERSION})
    endif ()

    configure_file(resources/descriptor.json.in descriptor.json @ONLY)

    set(CPACK_STRIP_FILES "${INSTALL_DIR_RUNTIME}/${CMAKE_PROJECT_NAME}")

    include(CPack)
endif ()

# ######################################################################### #
