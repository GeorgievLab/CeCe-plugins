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
#include "Module.hpp"

// C++
#include <random>
#include <string>
#include <cmath>

// CeCe
#include "cece/core/Log.hpp"
#include "cece/core/Assert.hpp"
#include "cece/core/Exception.hpp"
#include "cece/core/StringStream.hpp"
#include "cece/core/TimeMeasurement.hpp"
#include "cece/core/UnitIo.hpp"
#include "cece/object/Object.hpp"
#include "cece/simulator/Simulation.hpp"
#include "cece/simulator/TimeMeasurement.hpp"

// Plugins
#include "../streamlines/Module.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace object_streamlines_generator {

/* ************************************************************************ */

namespace {

/* ************************************************************************ */

// Random engine
std::random_device g_rd;

/* ************************************************************************ */

std::default_random_engine g_gen(g_rd());

/* ************************************************************************ */

/**
 * @brief Parse list of active iterations.
 *
 * @param str
 *
 * @return
 */
DynamicArray<IterationRange> parseActive(String str)
{
    DynamicArray<IterationRange> res;

    InStringStream iss(std::move(str));

    while (true)
    {
        IterationType it;

        if (!(iss >> it))
            break;

        if (iss.peek() == '-')
        {
            IterationType itEnd;
            iss.ignore();
            iss >> itEnd;

            res.emplace_back(it, itEnd);
        }
        else
        {
            // Single item range
            res.emplace_back(it);
        }
    }

    return res;
}

/* ************************************************************************ */

/**
 * @brief Check if iteration is in range.
 *
 * @param list
 * @param it
 *
 * @return
 */
bool inRange(const DynamicArray<IterationRange>& list, IterationType it)
{
    // No limitation
    if (list.empty())
        return true;

    for (const auto& p : list)
    {
        if (p.inRange(it))
            return true;
    }

    return false;
}

/* ************************************************************************ */

}

/* ************************************************************************ */

void Module::init()
{
    m_module = getSimulation().getModule("streamlines");

    if (!m_module)
        throw RuntimeException("[object-streamlines-generator] 'streamlines' module not found");

    // Get boundaries
    const auto& boundaries = m_module->getBoundaries();

    // Check generator objects
    for (const auto& desc : m_objects)
    {
        // Find boundary
        auto boundary = boundaries.find(desc.boundary);

        if (!boundary)
        {
            Log::warning("[object-streamlines-generator] Boundary '", desc.boundary, "' not found, will be ignored");
            continue;
        }

        if (boundary->getType() != streamlines::Boundary::Type::Inlet)
        {
            Log::warning("[object-streamlines-generator] Boundary '", desc.boundary, "' is not inlet");
            continue;
        }

        if (boundary->getInletFlow() == Zero)
        {
            Log::warning("[object-streamlines-generator] Boundary '", desc.boundary, "' have zero flow rate");
            continue;
        }

        if (boundary->getBlocks().empty())
        {
            Log::warning("[object-streamlines-generator] Boundary '" + boundary->getName() + "' with no blocks");
            continue;
        }
    }
}

/* ************************************************************************ */

void Module::update()
{
    auto& simulation = getSimulation();

    auto _ = measure_time("object-streamlines-generator", simulator::TimeMeasurement(simulation));

    // Get current iteration number and world size
    const auto iteration = simulation.getIteration();
    const auto worldSize = simulation.getWorldSize();
    const auto worldSizeH = worldSize * 0.5;

    // Get boundaries
    CECE_ASSERT(m_module);
    const auto& boundaries = m_module->getBoundaries();
    const auto& converter = m_module->getConverter();
    const auto& lattice = m_module->getLattice();

    const units::PositionVector start = worldSize * -0.5;
    const auto step = getSimulation().getWorldSize() / lattice.getSize();

    // Foreach generated objects
    for (const auto& desc : m_objects)
    {
        // Skip inactive generators
        if (!inRange(desc.active, iteration))
            continue;

        // Find boundary
        auto boundary = boundaries.find(desc.boundary);

        if (!boundary)
            continue;

        // Calculate spawn rate
        const auto rate = desc.concentration * boundary->getInletFlow();

        // Create object number + probability
        const auto number = rate * simulation.getTimeStep();

        // Number of created with 100% probability
        const auto baseCount = std::floor(number);

        // Probability of the remaining object
        const auto probability = number - baseCount;
        Assert(probability >= 0);
        Assert(probability <= 1);
        std::bernoulli_distribution distSpawn(probability);

        // Total number of spawned objects
        int count = baseCount + (distSpawn(g_gen) ? 1 : 0);

        // Get inlet blocks
        const auto& blocks = boundary->getBlocks();

        if (blocks.empty())
            continue;

        // Distrbution for block selection
        std::uniform_int_distribution<> blockDis(0, blocks.size() - 1);

        // Boundary position
        const auto position = boundary->getPosition();

        // Generate
        while (count--)
        {
            // Select block
            const auto& block = blocks[blockDis(g_gen)];

            units::Length xMin = Zero;
            units::Length xMax = Zero;
            units::Length yMin = Zero;
            units::Length yMax = Zero;

            // Block size
            const auto blockSize = converter.convertLength(block.first());
            const auto fix = blockSize * 0.1;

            switch (position)
            {
            case streamlines::Boundary::Position::Top:
                yMin = yMax = worldSizeH.getY() - desc.offset;
                xMin = converter.convertLength(block.first()) + fix - worldSizeH.getX();
                xMax = converter.convertLength(block.last()) - fix - worldSizeH.getX();
                break;

            case streamlines::Boundary::Position::Bottom:
                yMin = yMax = -(worldSizeH.getY() - desc.offset);
                xMin = converter.convertLength(block.first()) + fix - worldSizeH.getX();
                xMax = converter.convertLength(block.last()) - fix - worldSizeH.getX();
                break;

            case streamlines::Boundary::Position::Right:
                xMin = xMax = worldSizeH.getX() - desc.offset;
                yMin = converter.convertLength(block.first()) + fix - worldSizeH.getY();
                yMax = converter.convertLength(block.last()) - fix - worldSizeH.getY();
                break;

            case streamlines::Boundary::Position::Left:
                xMin = xMax = -(worldSizeH.getX() - desc.offset);
                yMin = converter.convertLength(block.first()) + fix - worldSizeH.getY();
                yMax = converter.convertLength(block.last()) - fix - worldSizeH.getY();
                break;

            default:
                throw InvalidArgumentException("Unknown boundary position");
            }

            // Uniform distribution
            std::uniform_real_distribution<RealType> xDist(xMin.value(), xMax.value());
            std::uniform_real_distribution<RealType> yDist(yMin.value(), yMax.value());

            // Create object
            auto object = simulation.createObject(desc.className);
            CECE_ASSERT(object);

            // Generate position vector
            const units::PositionVector pos{
                units::Length(xDist(g_gen)),
                units::Length(yDist(g_gen))
            };

            object->configure(desc.config, simulation);
            object->setPosition(pos);

            // Get coordinate to lattice
            const auto coord = Coordinate((pos - start) / step);

            if (!lattice.inRange(coord))
                continue;

            // Extract velocity from LB to match the flow
            const auto velLB = lattice[coord].computeVelocity();

            // Obtain physical velocity.
            auto vel = converter.convertVelocity(velLB);

            // Set object velocity
            object->setVelocity(vel);
        }
    }

}

/* ************************************************************************ */

void Module::loadConfig(const config::Configuration& config)
{
    // Configure parent
    module::Module::loadConfig(config);

    for (auto&& cfg : config.getConfigurations("object"))
    {
        ObjectDesc desc;
        desc.className = cfg.get("class");
        desc.boundary = cfg.get("boundary");
        desc.offset = cfg.get("offset", desc.offset);
        desc.concentration = cfg.get<units::NumberConcentration>("concentration");
        desc.active = parseActive(cfg.get("active", String{}));
        desc.config = cfg.toMemory();

        add(std::move(desc));
    }
}

/* ************************************************************************ */

void Module::storeConfig(config::Configuration& config) const
{
    module::Module::storeConfig(config);

    for (const auto& object : m_objects)
    {
        auto objectConfig = config.addConfiguration("object");

        objectConfig.set("class", object.className);
        objectConfig.set("boundary", object.boundary);
        objectConfig.set("concentration", object.concentration);
        //objectConfig.set("active", object.active);

        // TODO: store remaining config
        // desc.config
    }
}

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
