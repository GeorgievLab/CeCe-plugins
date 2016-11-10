/* ************************************************************************ */
/* Georgiev Lab (c) 2015-2016                                               */
/* ************************************************************************ */
/* Department of Cybernetics                                                */
/* Faculty of Applied Sciences                                              */
/* University of West Bohemia in Pilsen                                     */
/* ************************************************************************ */
/*                                                                          */
/* This file is part of CeCe.                                               */
/*                                                                          */
/* CeCe is free software: you can redistribute it and/or modify             */
/* it under the terms of the GNU General Public License as published by     */
/* the Free Software Foundation, either version 3 of the License, or        */
/* (at your option) any later version.                                      */
/*                                                                          */
/* CeCe is distributed in the hope that it will be useful,                  */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU General Public License for more details.                             */
/*                                                                          */
/* You should have received a copy of the GNU General Public License        */
/* along with CeCe.  If not, see <http://www.gnu.org/licenses/>.            */
/*                                                                          */
/* ************************************************************************ */

#pragma once

/* ************************************************************************ */

// CeCe
#include "cece/core/Real.hpp"
#include "cece/core/ViewPtr.hpp"
#include "cece/core/Map.hpp"
#include "cece/core/Parameters.hpp"
#include "cece/simulator/Simulation.hpp"

#include "../cell/CellBase.hpp"

// Plugin
#include "Diffusion.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace stochastic_reactions {

/* ************************************************************************ */

/**
 * @brief Container for important pointers: current Cell and Diffusion.
 */
struct Context
{
    simulator::Simulation* simulation = nullptr;
    plugin::diffusion::Module* diffusion = nullptr;
    plugin::cell::CellBase* cell = nullptr;
    const DynamicArray<plugin::diffusion::Module::Coordinate>* coords = nullptr;
    const core::Map<String, RealType>* arguments = nullptr;

    Context() = default;

    Context(
        simulator::Simulation* simulation,
        plugin::diffusion::Module* d,
        plugin::cell::CellBase* c,
        const DynamicArray<plugin::diffusion::Module::Coordinate>* cs,
        const core::Map<String, RealType>* args)
        : simulation(simulation)
        , diffusion(d)
        , cell(c)
        , coords(cs)
        , arguments(args)
    {
        // Nothing to do
    }
};

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
