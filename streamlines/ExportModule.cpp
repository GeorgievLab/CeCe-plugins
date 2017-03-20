/* ************************************************************************ */
/* Georgiev Lab (c) 2015-2017                                               */
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
            "d0", "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8", "dynamics"
        );
    }
    else
    {
        writeHeader(
            "iteration", "totalTime", "x", "y", "vx", "vy", "rho", "dynamics"
        );
    }
}

/* ************************************************************************ */

void ExportModule::update()
{
    Assert(m_module);

    const auto& sim = getSimulation();

    for (auto&& c : range(m_module->getLatticeSize()))
    {
        // Do not save data with no dynamics
        if (m_module->getDynamics(c) == Dynamics::None)
            continue;

        const auto vel = m_module->getVelocity(c);
        const auto p = m_module->getPressure(c);
        const int dynamics = static_cast<int>(m_module->getDynamics(c));

        if (isPopulationsExported())
        {
            const auto& data = m_module->getDistributions(c);

            writeRecord(
                sim.getIteration(),
                sim.getTotalTime().value(),
                c.getX(),
                c.getY(),
                vel.getX().value(),
                vel.getY().value(),
                p.value(),
                data[0],
                data[1],
                data[2],
                data[3],
                data[4],
                data[5],
                data[6],
                data[7],
                data[8],
                dynamics
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
                p.value(),
                dynamics
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
