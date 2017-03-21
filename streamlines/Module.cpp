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
#include "Module.hpp"

// C++
#include <algorithm>
#include <limits>
#include <cmath>

// CeCe
#include "cece/core/Assert.hpp"
#include "cece/core/constants.hpp"
#include "cece/core/StaticArray.hpp"
#include "cece/core/DynamicArray.hpp"
#include "cece/core/Exception.hpp"
#include "cece/core/VectorRange.hpp"
#include "cece/core/constants.hpp"
#include "cece/core/FileStream.hpp"
#include "cece/core/CsvFile.hpp"
#include "cece/core/TimeMeasurement.hpp"
#include "cece/core/BinaryInput.hpp"
#include "cece/core/BinaryOutput.hpp"
#include "cece/core/ShapeToGrid.hpp"
#include "cece/core/Log.hpp"
#include "cece/core/UnitIo.hpp"
#include "cece/object/Object.hpp"
#include "cece/config/Configuration.hpp"
#include "cece/simulator/TimeMeasurement.hpp"
#include "cece/simulator/Simulation.hpp"
#include "cece/simulator/Visualization.hpp"

// Plugin
#include "FactoryManager.hpp"
#include "cpu/Lattice.hpp"
#include "gpu/Lattice.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines {

/* ************************************************************************ */

namespace {

/* ************************************************************************ */

constexpr StaticArray<char, 5> FILE_GUARD{{'C', 'E', 'S', 'L', '\0'}};

/* ************************************************************************ */

/**
 * @brief      Parse array of values from string.
 *
 * @param      str   Source string.
 *
 * @tparam     T     Element type.
 * @tparam     N     Number of result elements.
 *
 * @return     Array of result elements.
 */
template<typename T, size_t N>
StaticArray<T, N> split(String str)
{
    StaticArray<T, N> array;
    InStringStream is(std::move(str));

    // Read values
    for (size_t i = 0; i < N; ++i)
        is >> std::skipws >> array[i];

    return array;
}

/* ************************************************************************ */

/**
 * @brief      Calculare edge block.
 *
 * @param      offset  Center offset.
 * @param      size    Boundary size.
 * @param      count   Total lattice size.
 *
 * @return     The range of indices for boundary.
 */
Range<int> getBoundaryBlock(int offset, int size, int count) noexcept
{
    if (size == 0)
        size = count;

    const int off = count / 2 + offset;
    const int siz = size / 2;

    return {off - siz, off + siz};
}

/* ************************************************************************ */

/**
 * @brief      Split boundary block to multiple blocks according to fluid
 *             dynamics in the lattice.
 *
 * @param      lattice        The lattice
 * @param      range          The range
 * @param      fluidDynamics  The fluid dynamics
 * @param      coordFn        The coordinate function
 *
 * @tparam     CoordFn        Coordinate function.
 *
 * @return     List of index ranges.
 */
template<typename CoordFn>
DynamicArray<Range<int>> detectBoundaryBlocks(
    const Lattice& lattice, Range<int> range, CoordFn coordFn
) noexcept
{
    DynamicArray<Range<int>> blocks;

    // Foreach boundary edge
    for (int i = range.first(); i < range.last(); ++i)
    {
        const auto dynamics1 = lattice.getDynamics(coordFn(i));

        if (dynamics1 != Dynamics::Fluid)
            continue;

        const int start = i++;

        // Inner block
        for (; i < range.last(); ++i)
        {
            const auto dynamics2 = lattice.getDynamics(coordFn(i));

            if (dynamics2 != Dynamics::Fluid)
                break;
        }

        // Insert block
        blocks.push_back({start, i});
    }

    return blocks;
}

/* ************************************************************************ */

/**
 * @brief      Set boundary blocks to lattice.
 *
 * @param      lattice  The lattice
 * @param      blocks   The blocks
 * @param      coordFn  The coordinate function
 * @param      setFn    The set function
 *
 * @tparam     CoordFn  Coordinate function.
 * @tparam     SetFn    Apply function.
 */
template<typename CoordFn, typename SetFn>
void setBoundaryBlocks(
    Lattice& lattice, const DynamicArray<Range<int>>& blocks,
    CoordFn coordFn, SetFn setFn
) noexcept
{
    // Foreach blocks
    for (const auto& block : blocks)
    {
        for (auto i = block.first(); i < block.last(); ++i)
        {
            const auto c = coordFn(i);

            // Set boundary dynamics
            if (lattice.getDynamics(c) == Dynamics::Fluid)
                setFn(c, block.getSize());
        }
    }
}

/* ************************************************************************ */

}

/* ************************************************************************ */

Module::Module(simulator::Simulation& simulation)
    : module::Module(simulation)
{
    // Register CPU builder
    FactoryManager::getInstance().createForLattice<cpu::Lattice>("cpu");

    // Register GPU builder
    FactoryManager::getInstance().createForLattice<gpu::Lattice>("gpu");
}

/* ************************************************************************ */

Module::~Module() = default;

/* ************************************************************************ */

ViewPtr<const Boundary> Module::findBoundary(StringView name) const noexcept
{
    auto it = std::find_if(m_config.boundaries.begin(), m_config.boundaries.end(), [&](const Boundary& b) {
        return b.getName() == name;
    });

    if (it != m_config.boundaries.end())
        return &*it;

    return nullptr;
}

/* ************************************************************************ */

DynamicArray<Pair<units::PositionVector, units::PositionVector>>
Module::getBoundaryBlocks(StringView name) const noexcept
{
    auto boundary = findBoundary(name);

    // No boundary with that name
    if (boundary == nullptr)
        return {};

    DynamicArray<Pair<units::PositionVector, units::PositionVector>> res;

    // No valid boundary type
    if (boundary->getType() == Boundary::Type::None)
        return {};

    const auto gridSize = m_lattice->getSize();
    const auto offset = m_converter.convertLength(boundary->getOffset());
    const auto size = m_converter.convertLength(boundary->getSize());

    const units::PositionVector start = getSimulation().getWorldSize() * -0.5;

    switch (boundary->getPosition())
    {
    case Boundary::Position::Top:
    {
        // Boundary block
        const auto block = getBoundaryBlock(offset, size, gridSize.getWidth());

        auto coordFn = [&](int i) {
            return Lattice::CoordinateType(i, gridSize.getHeight() - 1);
        };

        // Top edge
        const auto blocks = detectBoundaryBlocks(*m_lattice, block, coordFn);

        for (auto block : blocks)
        {
            const auto pos1 = start + m_converter.convertLength(Vector<RealType>(coordFn(block.first())));
            const auto pos2 = start + m_converter.convertLength(Vector<RealType>(coordFn(block.last())));

            res.push_back({pos1, pos2});
        }

        break;
    }

    case Boundary::Position::Bottom:
    {
        // Boundary block
        const auto block = getBoundaryBlock(offset, size, gridSize.getWidth());

        auto coordFn = [&](int i) {
            return Lattice::CoordinateType(i, 0);
        };

        // Bottom edge
        const auto blocks = detectBoundaryBlocks(*m_lattice, block, coordFn);

        for (auto block : blocks)
        {
            const auto pos1 = start + m_converter.convertLength(Vector<RealType>(coordFn(block.first())));
            const auto pos2 = start + m_converter.convertLength(Vector<RealType>(coordFn(block.last())));

            res.push_back({pos1, pos2});
        }

        break;
    }

    case Boundary::Position::Left:
    {
        // Boundary block
        const auto block = getBoundaryBlock(offset, size, gridSize.getHeight());

        auto coordFn = [&](int i) {
            return Lattice::CoordinateType(0, i);
        };

        // Left edge
        const auto blocks = detectBoundaryBlocks(*m_lattice, block, coordFn);

        for (auto block : blocks)
        {
            const auto pos1 = start + m_converter.convertLength(Vector<RealType>(coordFn(block.first())));
            const auto pos2 = start + m_converter.convertLength(Vector<RealType>(coordFn(block.last())));

            res.push_back({pos1, pos2});
        }

        break;
    }

    case Boundary::Position::Right:
    {
        // Boundary block
        const auto block = getBoundaryBlock(offset, size, gridSize.getHeight());

        auto coordFn = [&](int i) {
            return Lattice::CoordinateType(gridSize.getWidth() - 1, i);
        };

        // Right edge
        const auto blocks = detectBoundaryBlocks(*m_lattice, block, coordFn);

        for (auto block : blocks)
        {
            const auto pos1 = start + m_converter.convertLength(Vector<RealType>(coordFn(block.first())));
            const auto pos2 = start + m_converter.convertLength(Vector<RealType>(coordFn(block.last())));

            res.push_back({pos1, pos2});
        }

        break;
    }
    }

    return res;
}

/* ************************************************************************ */

Lattice::SizeType Module::getLatticeSize() const noexcept
{
    return Lattice::SizeType(
        getSimulation().getWorldSize() / m_converter.getCharLength() * m_converter.getNumberNodes()
    ) + Lattice::SizeType(1, 1);
}

/* ************************************************************************ */

void Module::init(AtomicBool& flag)
{
    // Get calculated lattice size
    const auto size = getLatticeSize();

    // Create a lattice
    createLattice(size);

    if (m_lattice->getSize() == Zero)
        throw InvalidArgumentException("[streamlines] Zero size lattice");

    // Allocate moving obstacle map
    m_movingObstacleMap.resize(size);

    Log::info("[streamlines] Viscosity: ", m_converter.getKinematicViscosity(), " um2/s");
    Log::info("[streamlines] Max object speed: ", getSimulation().getMaxObjectTranslation(), " um/it");
    Log::info("[streamlines] Char. length: ", m_converter.getCharLength(), " um");
    Log::info("[streamlines] Char. time: ", m_converter.getCharTime(), " s");
    Log::info("[streamlines] Char. speed: ", m_converter.getCharVelocity(), " um/s");
    Log::info("[streamlines] Number of nodes: ", m_converter.getNumberNodes());
    Log::info("[streamlines] Number of time steps: ", m_converter.getNumberSteps());
    Log::info("[streamlines] Re: ", m_converter.getRe());
    Log::info("[streamlines] ## Lattice ##");
    Log::info("[streamlines] Tau: ", m_converter.getTau());
    Log::info("[streamlines] Omega: ", m_converter.getOmega());
    Log::info("[streamlines] Size: (", size.getWidth(), "; ", size.getHeight(), ")");
    Log::info("[streamlines] Viscosity: ", m_converter.getViscosity());

    if (m_converter.getTau() > 5.0)
    {
        Log::warning("[streamlines] Relaxation parameter Tau is too large, unwanted behaviour can be observed. "
            "Increase number of time steps or decrease viscosity."
        );
    }

    if (m_converter.getViscosity() == 0.0)
        Log::warning("[streamlines] Zero viscosity means unstable simulation");

    Log::info("[streamlines] Boundaries: ", m_config.boundaries.size());

    // Initialize boundary positions
    for (auto& boundary : m_config.boundaries)
    {
        Log::info("[streamlines] Boundary: ", boundary.getName());

        // TODO: define BC for lattice
    }

    // Init lattice from to simulation objects.
    updateDynamics();

    Log::info("[streamlines] Initialization...");

    bool initialized = false;

    // Loading precomputed streamlines
    if (!m_config.initFile.empty() && pathExists(m_config.initFile))
    {
        try
        {
            Log::info("[streamlines] Loading from external file: ", m_config.initFile);

            loadFromFile(m_config.initFile);
            initialized = true;
        }
        catch (const Exception& e)
        {
            Log::warning("[streamlines] ", e.what());
        }
    }

    if (!initialized)
    {
        // Initialize with zero velocity and default density
        m_lattice->initDefault();

        if (getInitIterations() == 0 && !isDynamic())
            Log::warning("[streamlines] Static simulation without precomputed streamlines!");

        // Initialization iterations
        for (IterationType it = 1; it <= getInitIterations(); it++)
        {
            if ((it % 100) == 0)
                Log::info("[streamlines] Initialization ", it, "/", getInitIterations());

            if (!flag)
            {
                Log::info("[streamlines] Initialization interrupted.");
                return;
            }

            // Update lattice
            m_lattice->update();
        }

        // Store precomputed streamlines
        if (!m_config.initFile.empty())
        {
            Log::info("[streamlines] Store streamlines into file: ", m_config.initFile);
            storeToFile(m_config.initFile);
            Log::info("[streamlines] Streamlines stored in file: ", m_config.initFile);
        }
    }

    Log::info("[streamlines] Initialization done.");
}

/* ************************************************************************ */

void Module::update()
{
    auto _ = measure_time("streamlines", simulator::TimeMeasurement(getSimulation()));

    // No recalculation
    if (isDynamic())
    {
        // Update dynamics
        updateDynamics();

        // Update lattice
        m_lattice->update(getInnerIterations());
    }

    // Update simulation objects
    updateObjects();
}

/* ************************************************************************ */

void Module::loadConfig(const config::Configuration& config)
{
    // Configure parent
    module::Module::loadConfig(config);

    // Lattice mode
    m_config.mode = config.get("mode", m_config.mode);

    // Set streamlines dynamicity
    setDynamic(config.get("dynamic", isDynamic()));

    // Number of init iterations
    setInitIterations(config.get("init-iterations", getInitIterations()));

    // Number of inner iterations
    setInnerIterations(config.get("inner-iterations", getInnerIterations()));

    // Load boundaries config
    {
        auto initDefault = [&] {
            if (!m_config.boundaries.empty())
                return;

            m_config.boundaries.resize(4);

            for (int i = 0; i < 4; ++i)
                m_config.boundaries[i].setPosition(static_cast<Boundary::Position>(i));
        };

        // Layout
        if (config.has("layout"))
        {
            // Create default
            initDefault();

            const auto layout = split<String, 4>(config.get("layout"));

            for (size_t i = 0; i < layout.size(); ++i)
            {
                if (layout[i] == "inlet")
                    m_config.boundaries[i].setType(Boundary::Type::Inlet);
                else if (layout[i] == "outlet")
                    m_config.boundaries[i].setType(Boundary::Type::Outlet);
                else
                    m_config.boundaries[i].setType(Boundary::Type::None);
            }
        }

        // Inlet velocities
        if (config.has("inlet-velocity"))
        {
            initDefault();

            const auto velocity = config.get<units::Velocity>("inlet-velocity");

            for (auto& boundary : m_config.boundaries)
                boundary.setInletVelocity(velocity);
        }
        else if (config.has("inlet-velocities"))
        {
            initDefault();

            const auto velocities = split<units::Velocity, 4>(config.get("inlet-velocities"));

            for (size_t i = 0; i < velocities.size(); ++i)
                m_config.boundaries[i].setInletVelocity(velocities[i]);
        }

        if (config.has("inlet-type"))
        {
            initDefault();

            const auto typeStr = config.get<String>("inlet-type");
            Boundary::InletProfileType type = Boundary::InletProfileType::Auto;

            if (typeStr == "constant")
                type = Boundary::InletProfileType::Constant;

            for (auto& boundary : m_config.boundaries)
                boundary.setInletProfileType(type);
        }
        else if (config.has("inlet-types"))
        {
            initDefault();

            const auto types = split<String, 4>(config.get("inlet-types"));

            for (size_t i = 0; i < types.size(); ++i)
            {
                Boundary::InletProfileType type = Boundary::InletProfileType::Auto;

                if (types[i] == "constant")
                    type = Boundary::InletProfileType::Constant;

                m_config.boundaries[i].setInletProfileType(type);
            }
        }

        // Inlets
        for (auto&& boundaryConfig : config.getConfigurations("boundary"))
        {
            Boundary boundary;
            boundary.loadConfig(boundaryConfig);
            m_config.boundaries.push_back(std::move(boundary));
        }
    }

    // Units converter
    {
        // Viscosity
        m_converter.setKinematicViscosity(config.get("kinematic-viscosity", m_converter.getKinematicViscosity()));

        if (m_converter.getKinematicViscosity() == Zero)
            throw config::Exception("Kinematic viscosity ('kinematic-viscosity') cannot be zero");

        // Characteristic length & time
        m_converter.setCharLength(config.get("char-length", m_converter.getCharLength()));
        m_converter.setNumberNodes(config.get("number-nodes", m_converter.getNumberNodes()));

        // Set number of time steps
        m_converter.setCharTime(getSimulation().getTimeStep());
        m_converter.setNumberSteps(getInnerIterations());

        // Set fluid density
        m_converter.setCharDensity(config.get("density", m_converter.getCharDensity()));
    }

    // Enable dynamic object obstacles
    setDynamicObjectsObstacles(config.get("dynamic-object-obstacles", isDynamicObjectsObstacles()));

    m_config.circleObstacleScale = config.get("dynamic-object-circle-scale", m_config.circleObstacleScale);

#ifdef CECE_RENDER
    // Visualization
    {
        m_render.layerDynamics = config.get("layer-dynamics", m_render.layerDynamics);
        m_render.layerMagnitude = config.get("layer-magnitude", m_render.layerMagnitude);
        m_render.layerDensity = config.get("layer-density", m_render.layerDensity);
    }
#endif

    // Maximal force
    m_config.maxForce = config.get("max-force", m_config.maxForce);

    // Get initialization file
    if (config.has("init-file"))
    {
        // Get file name
        auto file = config.get("init-file");

        // In case of %temp%
        if (file.substr(0, 6) == "%temp%")
            m_config.initFile = tempDirectory() / file.substr(6);
        else
            m_config.initFile = file;
    }
}

/* ************************************************************************ */

void Module::storeConfig(config::Configuration& config) const
{
    module::Module::storeConfig(config);

    // TODO: implement
}

/* ************************************************************************ */

#ifdef CECE_RENDER
void Module::draw(const simulator::Visualization& visualization, render::Context& context)
{
    const auto size = m_lattice->getSize();
    const bool drawDynamics = visualization.isEnabled(m_render.layerDynamics);
    const bool drawMagnitude = visualization.isEnabled(m_render.layerMagnitude);
    const bool drawDensity = visualization.isEnabled(m_render.layerDensity);

    if (drawDynamics && !m_render.dynamics)
        m_render.dynamics.create(context, size);

    if (drawMagnitude && !m_render.magnitude)
        m_render.magnitude.create(context, size);

    if (drawDensity && !m_render.density)
        m_render.density.create(context, size);

    const RenderState& state = m_render.state.getFront();

    if (drawDynamics && m_render.dynamics)
        m_render.dynamics->setImage(state.imageDynamics);

    if (drawMagnitude && m_render.magnitude)
        m_render.magnitude->setImage(state.imageMagnitude);

    if (drawDensity && m_render.density)
        m_render.density->setImage(state.imageDensity);

    // Draw color grid
    context.matrixPush();
    context.matrixScale(state.scale);

    if (drawDynamics && m_render.dynamics)
        m_render.dynamics->draw(context);

    if (drawMagnitude && m_render.magnitude)
        m_render.magnitude->draw(context);

    if (drawDensity && m_render.density)
        m_render.density->draw(context);

    context.matrixPop();
}
#endif

/* ************************************************************************ */

#ifdef CECE_RENDER
void Module::drawStoreState(const simulator::Visualization& visualization)
{
    const auto size = m_lattice->getSize();
    const bool drawDynamics = visualization.isEnabled(m_render.layerDynamics);
    const bool drawMagnitude = visualization.isEnabled(m_render.layerMagnitude);
    const bool drawDensity = visualization.isEnabled(m_render.layerDensity);

    // Render state
    RenderState& state = m_render.state.getBack();

    state.scale = getSimulation().getWorldSize() / units::Length(1);

    // Resize image
    state.imageDynamics.resize(size);
    state.imageMagnitude.resize(size);
    state.imageDensity.resize(size);

    // Find min/max density
    Descriptor::DensityType rhoMin = std::numeric_limits<Descriptor::DensityType>::max();
    Descriptor::DensityType rhoMax = 0.0;
    RealType maxVel = 0.0;

    if (drawDynamics || drawMagnitude || drawDensity)
    {
        for (auto&& c : range(size))
        {
            if (m_lattice->isFluidDynamics(c))
            {
                const auto velocity = m_lattice->getVelocity(c);
                const auto density = m_lattice->getDensity(c);

                maxVel = std::max(maxVel, velocity.getLength());
                rhoMin = std::min(density, rhoMin);
                rhoMax = std::max(density, rhoMax);
            }
        }

        // Update texture
        for (auto&& c : range(size))
        {
            // Store dynamics type
            if (drawDynamics)
            {
                const auto dynamics = m_lattice->getDynamics(c);

                render::Color color;

                if (dynamics == Dynamics::Fluid)
                {
                    color = render::colors::BLACK;
                }
                else if (dynamics == Dynamics::Inlet)
                {
                    color = render::colors::RED;
                }
                else if (dynamics == Dynamics::Outlet)
                {
                    color = render::colors::BLUE;
                }
                else if (dynamics == Dynamics::None)
                {
                    color = render::colors::WHITE;
                }
                else
                {
                    color = render::colors::GREEN;
                }

                // Store color
                state.imageDynamics.set(c, color);
            }

            if (drawMagnitude)
            {
                const auto velocity = m_lattice->getVelocity(c);
                const auto dynamics = m_lattice->getDynamics(c);

                if (dynamics == Dynamics::Fluid)
                {
                    state.imageMagnitude.set(
                        c,
                        render::Color::fromGray(velocity.getLength() / maxVel)
                    );
                }
                else
                {
                    render::Color color = render::colors::BLACK;
                    color.setAlpha(0);

                    state.imageMagnitude.set(c, color);
                }
            }

            if (drawDensity)
            {
                const auto density = m_lattice->getDensity(c);
                const auto dynamics = m_lattice->getDynamics(c);

                if (dynamics == Dynamics::Fluid)
                {
                    state.imageDensity.set(
                        c,
                        render::Color::fromGray((density - rhoMin) / (rhoMax - rhoMin))
                    );
                }
                else
                {
                    render::Color color = render::colors::BLACK;
                    color.setAlpha(0);

                    state.imageDensity.set(c, color);
                }
            }
        }
    }
}
#endif

/* ************************************************************************ */

#ifdef CECE_RENDER
void Module::drawSwapState()
{
    m_render.state.swap();
}
#endif

/* ************************************************************************ */

void Module::updateDynamics()
{
    CECE_ASSERT(m_lattice);

    // Set whole scene as fluid
    for (auto&& c : range(m_lattice->getSize()))
        m_lattice->setFluidDynamics(c);

    const units::PositionVector start = getSimulation().getWorldSize() * -0.5;
    const auto step = getSimulation().getWorldSize() / m_lattice->getSize();

    Grid<char> movingObstacleMap;
    movingObstacleMap.resize(m_movingObstacleMap.getSize());

    // Foreach all cells
    for (auto& obj : getSimulation().getObjects())
    {
        const bool isDynamic = obj->getType() == object::Object::Type::Dynamic;

        // Ignore dynamic objects
        if (!isDynamicObjectsObstacles() && isDynamic)
            continue;

        // Get object position
        const auto pos = obj->getPosition() - start;

        // Check if position is in range
        if (!pos.inRange(Zero, getSimulation().getWorldSize()))
            continue;

        // Get grid position
        const auto center = Coordinate(pos / step);

        // Calculate object velocity in LB
        const auto velocity = m_converter.convertVelocity(obj->getVelocity());

        // In this case duplicate coordinates doesn't harm and calling
        // operation multiple times on same coordinate is faster than
        // sorting and erasing non-unique coordinates.

        // Map shapes to grid
        for (const auto& shape : obj->getShapes())
        {
            auto shp = shape;

            if (shp.getType() == ShapeType::Circle)
                shp.getCircle().radius *= m_config.circleObstacleScale;

            mapShapeToGrid(
                [&, this] (Coordinate&& coord) {
                    CECE_ASSERT(inLatticeRange(coord));
                    if (isDynamic && this->isDynamic())
                    {
                        m_lattice->setObjectDynamics(coord, velocity);
                        movingObstacleMap[coord] = true;
                        m_movingObstacleMap[coord] = false;
                    }
                    else
                    {
                        m_lattice->setWallDynamics(coord);
                    }
                },
                [] (Coordinate&& coord) {},
                shp, step, center, obj->getRotation(), m_lattice->getSize()
            );
        }
    }

    // Only for dynamic obstacles
    if (isDynamicObjectsObstacles())
    {
        // Swap maps
        std::swap(m_movingObstacleMap, movingObstacleMap);

        // Nodes with 'true' in movingObstacleMap have changed from obstacle to fluid
        for (const auto& c : range(m_lattice->getSize()))
        {
            if (!movingObstacleMap[c])
                continue;

            const Coordinate cs[4] = {
                Coordinate(c.getX() + 1, c.getY() + 0),
                Coordinate(c.getX() - 1, c.getY() + 0),
                Coordinate(c.getX() + 0, c.getY() + 1),
                Coordinate(c.getX() + 0, c.getY() - 1)
            };

            int count = 0;
            RealType density = 0;
            Vector<RealType> velocity = Zero;

            for (const auto c2 : cs)
            {
                if (!m_lattice->inRange(c2))
                    continue;

                if (!m_lattice->isFluidDynamics(c2))
                    continue;

                ++count;
                density += m_lattice->getDensity(c2);
                velocity += m_lattice->getVelocity(c2);
            }

            if (count)
            {
                m_lattice->setVelocityDensity(c, velocity / count, density / count);
            }
            else
            {
                m_lattice->setVelocityDensity(c, Zero, Descriptor::DEFAULT_DENSITY);
            }
        }
    }

    // Update boundaries blocks
    updateBoundaries();

    //if (isDynamicObjectsObstacles())
    //    m_lattice.fixupObstacles(Node::Dynamics::DynamicObstacle);
}

/* ************************************************************************ */

void Module::updateObjects()
{
    // Foreach objects
    for (auto& obj : getSimulation().getObjects())
        updateObject(*obj);
}

/* ************************************************************************ */

void Module::updateObject(object::Object& object)
{
    // Ignore static objects
    if (object.getType() != object::Object::Type::Dynamic)
        return;

    if (isDynamicObjectsObstacles() && isDynamic())
    {
        updateObjectDynamic(object);
    }
    else
    {
        updateObjectStatic(object);
    }
}

/* ************************************************************************ */

void Module::updateObjectStatic(object::Object& object)
{
    // Maximum object speed that is allowed by physical engine
    const auto maxSpeed = getSimulation().getMaxObjectTranslation() / getSimulation().getTimeStep();
    const auto maxSpeedSq = maxSpeed * maxSpeed;

    const units::PositionVector start = getSimulation().getWorldSize() * -0.5;
    const auto step = getSimulation().getWorldSize() / m_lattice->getSize();

    // Get mass local center
    const auto center = object.getMassCenterOffset();

    // Coefficient used in force calculation
    const auto forceCoefficient = 6 * constants::PI * m_converter.getKinematicViscosity() * object.getDensity();

    auto force = units::ForceVector{Zero};
    auto velocityObjEnv = units::VelocityVector{Zero};

    // Map shapes border to grid
    for (const auto& shape : object.getShapes())
    {
        // Only circle shapes are supported
        if (shape.getType() != ShapeType::Circle)
            continue;

        // Shape alias
        const auto& circle = shape.getCircle();

        // Transform from [-size / 2, size / 2] to [0, size] space
        const auto pos = object.getWorldPosition(circle.center) - start;

        // Check if position is in range
        if (!pos.inRange(Zero, getSimulation().getWorldSize()))
            continue;

        // Get coordinate to lattice
        const auto coord = Coordinate(pos / step);

        Vector<RealType> velocityLB = Zero;
        unsigned long count = 0;

        // Store velocity for each coordinate
        mapShapeBorderToGrid(
            [this, &velocityLB, &count] (Coordinate&& coord) {
                velocityLB += m_lattice->getVelocity(coord);
                ++count;
            },
            [] (Coordinate&& coord) { },
            shape, step, coord, m_lattice->getSize(), {}, 1
        );

        if (count == 0)
            continue;

        // Average
        velocityLB /= count;
        const units::VelocityVector velocityEnv = m_converter.convertVelocity(velocityLB);
        //const VelocityVector velocityEnv{m_inletVelocities[0], Zero};

        if (velocityEnv.getLengthSquared() > maxSpeedSq)
        {
            OutStringStream oss;
            oss <<
                "[streamlines] Physical engine can't handle environment "
                "velocity (" << velocityEnv.getLength() << " um/s < " << maxSpeed << " um/s). "
                "Decrease inlet velocities or change topology."
            ;

            Log::warning(oss.str());
            //throw RuntimeException(oss.str());
        }

        velocityObjEnv += velocityEnv;

        // Shape radius
        const auto radius = circle.radius;
        // Distance from mass center
        const auto offset = circle.center - center;

        // Angular velocity
        const auto omega = object.getAngularVelocity();

        // Calculate shape global velocity
        const auto velocity = object.getVelocity() + cross(omega, offset);

        // Difference between environment velocity and shape velocity
        const auto dv = velocityEnv - velocity;

        // Same velocity
        if (dv == Zero)
            continue;

        // Add force from shape
        force += forceCoefficient * radius * dv;
    }

    Assert(object.getShapes().size() > 0);
    velocityObjEnv /= object.getShapes().size();

    // Difference between velocities
    const auto dv = velocityObjEnv - object.getVelocity();

    // Same velocity
    if (dv == Zero)
        return;

    // Maximum impulse
    const auto impulseMax = object.getMass() * dv;

    // Calculate linear impulse from shapes
    auto impulse = force * getSimulation().getTimeStep();

    // Impulse is to big
    if (impulse.getLengthSquared() > impulseMax.getLengthSquared())
    {
        const RealType ratio = impulseMax.getLength() / impulse.getLength();
        impulse *= ratio;
    }

    // Apply impulse
    object.applyLinearImpulse(impulse);
}

/* ************************************************************************ */

void Module::updateObjectDynamic(object::Object& object)
{
    // Maximum object speed that is allowed by physical engine
    const auto maxSpeed = getSimulation().getMaxObjectTranslation() / getSimulation().getTimeStep();
    const auto maxSpeedSq = maxSpeed * maxSpeed;

    const units::PositionVector start = getSimulation().getWorldSize() * -0.5f;
    const auto step = getSimulation().getWorldSize() / m_lattice->getSize();

    // Get mass local center
    const auto center = object.getMassCenterOffset();

    // Coefficient used in force calculation
    const auto forceCoefficient = 6 * constants::PI * m_converter.getKinematicViscosity() * object.getDensity();

    auto force = units::ForceVector{Zero};

    // Map shapes border to grid
    for (const auto& shape : object.getShapes())
    {
        // Only circle shapes are supported
        if (shape.getType() != ShapeType::Circle)
            continue;

        // Shape alias
        const auto& circle = shape.getCircle();

        // Transform from [-size / 2, size / 2] to [0, size] space
        const auto pos = object.getWorldPosition(circle.center) - start;

        // Check if position is in range
        if (!pos.inRange(Zero, getSimulation().getWorldSize()))
            continue;

        // Get coordinate to lattice
        const auto c = Coordinate(pos / step);

        Vector<RealType> forceLB = Zero;

        // Force function
        auto forceEval = [this, &forceLB] (Coordinate&& coord) {

            Vector<RealType> force = Zero;

            // Foreach all directions
            for (Descriptor::PopIndexType iPop = 1; iPop < Descriptor::SIZE; ++iPop)
            {
                const auto iOpp = Descriptor::DIRECTION_OPPOSITES[iPop];
                const auto ei = Descriptor::DIRECTION_VELOCITIES[iOpp];

                // Neighbour coordinate & node
                const auto nc = coord + Coordinate(ei);

                // Skip
                if (!m_lattice->inRange(nc))
                    continue;

                // Not fluid dynamics
                if (!m_lattice->isFluidDynamics(nc))
                    continue;

                // Get distribution functions
                auto nn = m_lattice->getDistributions(nc);

                force += -ei * (nn[iOpp] + nn[iPop]);
            }

            forceLB += force;
        };

        auto shp = shape;

        if (shp.getType() == ShapeType::Circle)
            shp.getCircle().radius *= m_config.circleObstacleScale;

        // Store force for each coordinate
        mapShapeToGrid(
            forceEval,
            [] (Coordinate&& coord) { },
            shp, step, c, object.getRotation(), m_lattice->getSize()
        );

        // No force
        if (forceLB == Zero)
            continue;

        // Shape force
        const auto forceShape = m_converter.convertForce(forceLB);

        // Add force from shape
        force += forceShape;
    }

    // Force coefficient
    // Converter uses calculated number of steps but the simulation works with inner iterations.
    const auto coeff =
        (static_cast<RealType>(getInnerIterations()) * static_cast<RealType>(getInnerIterations())) /
        (static_cast<RealType>(m_converter.getNumberSteps()) * static_cast<RealType>(m_converter.getNumberSteps()))
    ;

    force *= coeff;

    //Log::info("Force: (", force.getX(), ", ", force.getY(), ")");

    // Limit force
    if (force.getLengthSquared() > m_config.maxForce.getLengthSquared())
    {
        const RealType ratio = m_config.maxForce.getLength() / force.getLength();
        force *= ratio;
    }

    // Calculate linear impulse from shapes
    const auto impulse = force * getSimulation().getTimeStep();

    // Apply impulse
    object.applyLinearImpulse(impulse);
}

/* ************************************************************************ */

void Module::updateBoundaries()
{
    for (const auto& boundary : m_config.boundaries)
        updateBoundary(boundary);
}

/* ************************************************************************ */

void Module::updateBoundary(const Boundary& boundary)
{
    // No boundary
    if (boundary.getType() == Boundary::Type::None)
        return;

    const auto gridSize = m_lattice->getSize();
    const auto offset = m_converter.convertLength(boundary.getOffset());
    const auto size = m_converter.convertLength(boundary.getSize());

    auto inletFlowVelocity = [&] (RealType width) -> units::Velocity {
        const auto len = m_converter.convertLength(width);
        const auto area = len * units::Length(1);
        return boundary.getInletFlow() / area;
    };

    // Calculate inlet velocity
    auto inletVelocity = [&] (Lattice::CoordinateType coord, int width) -> units::VelocityVector
    {
        auto type = boundary.getInletProfileType();

        if (type == Boundary::InletProfileType::Auto)
            type = Boundary::InletProfileType::Constant;

        if (type == Boundary::InletProfileType::Constant)
        {
            // Calculate velocity
            const units::Velocity velocity = boundary.getInletFlow() != Zero
                ? inletFlowVelocity(width)
                : boundary.getInletVelocity()
            ;

            switch (boundary.getPosition())
            {
            case Boundary::Position::Top:
                return {Zero, -velocity};

            case Boundary::Position::Bottom:
                return {Zero, velocity};

            case Boundary::Position::Left:
                return {velocity, Zero};

            case Boundary::Position::Right:
                return {-velocity, Zero};

            default:
                Assert(false && "Invalid boundary position");
                std::abort();
            }
        }

        return Zero;
    };

    // Boundary set function
    auto setFnInlet = [&] (Lattice::CoordinateType c, int width) {
        const auto iu = inletVelocity(c, width);
        const auto u = m_converter.convertVelocity(iu);
        m_lattice->setInletDynamics(c, u);
    };

    // Boundary set function
    auto setFnOutlet = [&] (Lattice::CoordinateType c, int width) {
        m_lattice->setOutletDynamics(c, 1.0);
    };


    switch (boundary.getPosition())
    {
    case Boundary::Position::Top:
    {
        // Boundary block
        const auto block = getBoundaryBlock(offset, size, gridSize.getWidth());

        auto coordFn = [&](int i) {
            return Lattice::CoordinateType(i, gridSize.getHeight() - 1);
        };

        // Top edge
        const auto blocks = detectBoundaryBlocks(*m_lattice, block, coordFn);

        switch (boundary.getType())
        {
        case Boundary::Type::Inlet:
            setBoundaryBlocks(*m_lattice, blocks, coordFn, setFnInlet);
            break;

        case Boundary::Type::Outlet:
            setBoundaryBlocks(*m_lattice, blocks, coordFn, setFnOutlet);
            break;

        case Boundary::Type::None:
            break;
        }

        break;
    }

    case Boundary::Position::Bottom:
    {
        // Boundary block
        const auto block = getBoundaryBlock(offset, size, gridSize.getWidth());

        auto coordFn = [&](int i) {
            return Lattice::CoordinateType(i, 0);
        };

        // Bottom edge
        const auto blocks = detectBoundaryBlocks(*m_lattice, block, coordFn);

        switch (boundary.getType())
        {
        case Boundary::Type::Inlet:
            setBoundaryBlocks(*m_lattice, blocks, coordFn, setFnInlet);
            break;

        case Boundary::Type::Outlet:
            setBoundaryBlocks(*m_lattice, blocks, coordFn, setFnOutlet);
            break;

        case Boundary::Type::None:
            break;
        }

        break;
    }

    case Boundary::Position::Left:
    {
        // Boundary block
        const auto block = getBoundaryBlock(offset, size, gridSize.getHeight());

        auto coordFn = [&](int i) {
            return Lattice::CoordinateType(0, i);
        };

        // Left edge
        const auto blocks = detectBoundaryBlocks(*m_lattice, block, coordFn);

        switch (boundary.getType())
        {
        case Boundary::Type::Inlet:
            setBoundaryBlocks(*m_lattice, blocks, coordFn, setFnInlet);
            break;

        case Boundary::Type::Outlet:
            setBoundaryBlocks(*m_lattice, blocks, coordFn, setFnOutlet);
            break;

        case Boundary::Type::None:
            break;
        }

        break;
    }

    case Boundary::Position::Right:
    {
        // Boundary block
        const auto block = getBoundaryBlock(offset, size, gridSize.getHeight());

        auto coordFn = [&](int i) {
            return Lattice::CoordinateType(gridSize.getWidth() - 1, i);
        };

        // Right edge
        const auto blocks = detectBoundaryBlocks(*m_lattice, block, coordFn);

        switch (boundary.getType())
        {
        case Boundary::Type::Inlet:
            setBoundaryBlocks(*m_lattice, blocks, coordFn, setFnInlet);
            break;

        case Boundary::Type::Outlet:
            setBoundaryBlocks(*m_lattice, blocks, coordFn, setFnOutlet);
            break;

        case Boundary::Type::None:
            break;
        }

        break;
    }
    }
}

/* ************************************************************************ */

void Module::storeToFile(const FilePath& filename)
{
    OutFileStream ofs(filename.toString(), OutFileStream::binary);
    BinaryOutput out(ofs);

    // Write file guard
    out.write(FILE_GUARD);

    // Write lattice size
    out.write(m_lattice->getSize());

    // Store lattice hash
    out.write(calculateLatticeHash());

    // Write relaxation time
    out.write(m_converter.getTau());

    // Number of init iterations
    out.write(m_config.initIterations);

    for (auto&& c : range(m_lattice->getSize()))
    {
        // Write cell populations
        out.write(m_lattice->getDistributions(c));
    }
}

/* ************************************************************************ */

void Module::loadFromFile(const FilePath& filename)
{
    InFileStream ifs(filename.toString(), OutFileStream::binary);

    if (!ifs.is_open())
        throw InvalidArgumentException("Cannot load from file: File not found '" + filename.string() + "'");

    BinaryInput in(ifs);

    // Read guard
    StaticArray<char, FILE_GUARD.size()> guard;
    in.read(guard);

    if (guard != FILE_GUARD)
        throw InvalidArgumentException("Cannot load from file: File is not valid");

    // Read lattice size
    Lattice::SizeType size;
    in.read(size);

    if (size != m_lattice->getSize())
        throw InvalidArgumentException("Cannot load from file: different lattice sizes");

    std::size_t hash;
    in.read(hash);

    if (hash != calculateLatticeHash())
        throw InvalidArgumentException("Cannot load from file: different layout");

    RealType tau;
    in.read(tau);

    if (tau != m_converter.getTau())
        throw InvalidArgumentException("Cannot load from file: different relaxation times");

    decltype(m_config.initIterations) iterations;
    in.read(iterations);

    if (iterations != m_config.initIterations)
        throw InvalidArgumentException("Cannot load from file: different init iterations");

    // Foreach lattice
    for (auto&& c : range(m_lattice->getSize()))
    {
        Lattice::DistributionsType df;

        // Read node populations
        in.read(df);

        // Set distribution functions
        m_lattice->setDistributions(c, df);
    }
}

/* ************************************************************************ */

std::size_t Module::calculateLatticeHash() const noexcept
{
    // Based on djb2
    size_t hash = 5381;

    for (auto&& c : range(m_lattice->getSize()))
    {
        const auto dynamics = m_lattice->getDynamics(c);

        // Get value
        size_t val = static_cast<size_t>(dynamics);

        // Can overflow... don't care
        hash = ((hash << 5) + hash) + val; // hash * 33 + val
    }

    return hash;
}

/* ************************************************************************ */

void Module::createLattice(Lattice::SizeType size)
{
    try
    {
        // Try to create user selected lattice type
        m_lattice = FactoryManager::getInstance().createLattice(m_config.mode, size, m_converter.getOmega());
    }
    catch (const Exception& e)
    {
        Log::warning("[streamlines] Unable to create GPU (OpenCL) lattice, falling back to CPU: ", e.what());

        // CPU lattice fallback
        m_lattice = FactoryManager::getInstance().createLattice("cpu", size, m_converter.getOmega());
    }

    CECE_ASSERT(m_lattice);
}

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
