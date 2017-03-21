# ######################################################################### #
# Georgiev Lab (c) 2015-2017                                                #
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

# Try to find OpenCL
find_package(OpenCL)

# OpenCL found, use static linking
if (OpenCL_FOUND)
    include_directories(${OpenCL_INCLUDE_DIRS})
    add_definitions(-DOPENCL_STATIC)
    set(LIBRARIES ${OpenCL_LIBRARIES})
    set(SOURCES)
else ()
    set(LIBRARIES)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR})
    set(SOURCES library.hpp library.cpp)
endif ()

# ######################################################################### #

# Build plugin
build_plugin(streamlines-gpu
    SOURCES
        Plugin.cpp
        Lattice.cpp
        Lattice.hpp
        ${SOURCES}

    LIBRARIES
        ${LIBRARIES}

    PLUGINS_REQUIRED
        streamlines
)

# ######################################################################### #
