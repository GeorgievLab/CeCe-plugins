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

// Declaration
#include "ExportModule.hpp"

// CeCe
#include "cece/core/Assert.hpp"
#include "cece/core/IteratorRange.hpp"
#include "cece/core/VectorRange.hpp"
#include "cece/core/Exception.hpp"
#include "cece/config/Configuration.hpp"
#include "cece/simulator/Simulation.hpp"

// Plugin
#include "Module.hpp"
#include "Descriptor.hpp"
#include "NoDynamics.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines {

/* ************************************************************************ */

void ExportModule::loadConfig(const config::Configuration& config)
{
    module::ExportModule::loadConfig(config);

    setDensityExported(config.get("density", isDensityExported()));
    setPopulationsExported(config.get("populations", isPopulationsExported()));
}

/* ************************************************************************ */

void ExportModule::storeConfig(config::Configuration& config) const
{
    module::ExportModule::storeConfig(config);

    config.set("density", isDensityExported());
    config.set("populations", isPopulationsExported());
}

/* ************************************************************************ */

void ExportModule::init()
{
    // Get diffusion module
    m_module = getSimulation().getModule("streamlines");

    if (!m_module)
        throw RuntimeException("Streamlines module required!");

    module::ExportModule::init();

    // Write output header
    if (isPopulationsExported())
    {
        writeHeader(
            "iteration", "totalTime", "x", "y", "vx", "vy", "rho",
            "d0", "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8"
        );
    }
    else
    {
        writeHeader(
            "iteration", "totalTime", "x", "y", "vx", "vy", "rho"
        );
    }
}

/* ************************************************************************ */

void ExportModule::update()
{
    Assert(m_module);

    const auto& sim = getSimulation();
    const auto& lattice = m_module->getLattice();
    const auto& conv = m_module->getConverter();

    for (auto&& c : range(lattice.getSize()))
    {
        const auto& data = lattice[c];

        // Do not save data with no dynamics
        if (data.getDynamics() == NoDynamics::getInstance())
            continue;

        const auto vel = conv.convertVelocity(data.computeVelocity());

        if (isPopulationsExported())
        {
            writeRecord(
                sim.getIteration(),
                sim.getTotalTime().value(),
                c.getX(),
                c.getY(),
                vel.getX().value(),
                vel.getY().value(),
                data.computeDensity(),
                data[0],
                data[1],
                data[2],
                data[3],
                data[4],
                data[5],
                data[6],
                data[7],
                data[8]
            );
        }
        else
        {
            writeRecord(
                sim.getIteration(),
                sim.getTotalTime().value(),
                c.getX(),
                c.getY(),
                vel.getX().value(),
                vel.getY().value(),
                data.computeDensity()
            );
        }
    }

    flush();
}

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
