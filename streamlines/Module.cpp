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
#include "NoDynamics.hpp"
#include "BounceBackDynamics.hpp"
#include "BgkDynamics.hpp"
#include "ZouHeDynamics.hpp"
#include "MovingDynamics.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines {

/* ************************************************************************ */

namespace {

/* ************************************************************************ */

constexpr StaticArray<char, 5> FILE_GUARD{{'C', 'E', 'S', 'L', '\0'}};

/* ************************************************************************ */

}

/* ************************************************************************ */

Module::Module(simulator::Simulation& simulation)
    : module::Module(simulation)
    , m_boundaries(simulation)
{
    // Nothing to do
}

/* ************************************************************************ */

Module::~Module() = default;

/* ************************************************************************ */

void Module::init(AtomicBool& flag)
{
    // Prepare pool for moving obstacle object
    MovingDynamics::poolPrepare(m_lattice.getSize().getWidth() * m_lattice.getSize().getHeight());

    // Print simulation info
    printInfo();

    // Set fluid dynamics
    setFluidDynamics(createFluidDynamics());
    setWallDynamics(createWallDynamics());

    Log::info("[streamlines] Boundaries: ", m_boundaries.getCount());

    // Initialize boundary positions
    for (auto& boundary : m_boundaries)
    {
        Log::info("[streamlines] Boundary: ", boundary.getName());
        boundary.setDynamics(createBorderDynamics(boundary.getPosition()));
    }

    // Initialize boundaries
    m_boundaries.init(m_lattice, getFluidDynamics());

    // Obstacles
    updateObstacleMap();

    // Initialize lattice to equilibrium
    m_lattice.initEquilibrium();

    applyBoundaryConditions();

    Log::info("[streamlines] Initialization...");

    bool initialized = false;

    if (!m_initFile.empty() && pathExists(m_initFile))
    {
        try
        {
            Log::info("[streamlines] Loading from external file: ", m_initFile);

            loadFromFile(m_initFile);
            initialized = true;
        }
        catch (const Exception& e)
        {
            Log::warning("[streamlines] ", e.what());
        }
    }

    if (!initialized)
    {
        m_lattice.initEquilibrium();

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

            m_lattice.collideAndStream();

            // Apply boundary conditions
            applyBoundaryConditions();
        }

        // Store initialization
        if (!m_initFile.empty())
        {
            storeToFile(m_initFile);
        }
    }

    Log::info("[streamlines] Initialization done.");

    if (getInitIterations() == 0 && !isDynamic())
        Log::warning("[streamlines] Static simulation without initialization!");
}

/* ************************************************************************ */

void Module::update()
{
    Assert(m_lattice.getSize() != Zero);

    auto _ = measure_time("streamlines", simulator::TimeMeasurement(getSimulation()));

    // No recalculation
    if (isDynamic())
    {
        // Obstacles
        updateObstacleMap();

        // Compute inner iterations
        for (IterationType it = 0; it < getInnerIterations(); it++)
        {
            // Collide and propagate
            m_lattice.collideAndStream();

            // Apply boundary conditions
            applyBoundaryConditions();
        }
    }

    // Apply streamlines to world objects
    applyToObjects();
}

/* ************************************************************************ */

void Module::loadConfig(const config::Configuration& config)
{
    // Configure parent
    module::Module::loadConfig(config);

    // Set streamlines dynamicity
    setDynamic(config.get("dynamic", isDynamic()));

    // Number of init iterations
    setInitIterations(config.get("init-iterations", getInitIterations()));

    // Number of inner iterations
    setInnerIterations(config.get("inner-iterations", getInnerIterations()));

    // Load boundaries config
    m_boundaries.loadConfig(config);

    // Load converter configuration
    m_converter.setCharTime(getSimulation().getTimeStep());
    m_converter.loadConfig(config);

    // Obsolete grid
    auto gridSize = config.get<Lattice::Size>("grid", Lattice::Size{Zero});

    if (gridSize != Zero)
    {
        Log::warning("[streamlines] Config option 'grid' is obsolete!");

        // Grid size
        m_lattice.setSize(gridSize);

        // Compute characteristic length
        m_converter.setNumberNodes(
            gridSize.getWidth() / getSimulation().getWorldSize().getWidth() * m_converter.getCharLength()
        );
    }
    else
    {
        // Calculate lattice size
        const auto size = Lattice::Size(
            getSimulation().getWorldSize() / m_converter.getCharLength() * m_converter.getNumberNodes()
        );

        if (size == Zero)
            throw InvalidArgumentException("Lattice size cannot be zero");

        // Grid size
        m_lattice.setSize(size + Lattice::Size(1, 1));
    }

    // Allocate moving obstacle map
    m_movingObstacleMap.resize(m_lattice.getSize());

    // Enable dynamic object obstacles
    setDynamicObjectsObstacles(config.get("dynamic-object-obstacles", isDynamicObjectsObstacles()));

    m_circleObstacleScale = config.get("dynamic-object-circle-scale", m_circleObstacleScale);

#ifdef CECE_RENDER
    m_visualizationLayerDynamicsType = config.get("layer-dynamics", m_visualizationLayerDynamicsType);
    m_visualizationLayerMagnitude = config.get("layer-magnitude", m_visualizationLayerMagnitude);
    m_visualizationLayerDensity = config.get("layer-density", m_visualizationLayerDensity);
#endif

    // Maximal force
    m_maxForce = config.get("max-force", m_maxForce);

    // Get initialization file
    if (config.has("init-file"))
    {
        // Get file name
        auto file = config.get("init-file");

        // In case of %temp%
        if (file.substr(0, 6) == "%temp%")
            m_initFile = tempDirectory() / file.substr(6);
        else
            m_initFile = file;
    }
}

/* ************************************************************************ */

void Module::storeConfig(config::Configuration& config) const
{
    module::Module::storeConfig(config);

    config.set("init-iterations", getInitIterations());
    config.set("inner-iterations", getInnerIterations());

    // Store converter config
    m_converter.storeConfig(config);

    // TODO: init file
}

/* ************************************************************************ */

#ifdef CECE_RENDER
void Module::draw(const simulator::Visualization& visualization, render::Context& context)
{
    const auto size = m_lattice.getSize();
    const bool drawDynamicsType = visualization.isEnabled(m_visualizationLayerDynamicsType);
    const bool drawMagnitude = visualization.isEnabled(m_visualizationLayerMagnitude);
    const bool drawDensity = visualization.isEnabled(m_visualizationLayerDensity);

    if (drawDynamicsType && !m_drawableDynamicsType)
        m_drawableDynamicsType.create(context, size);

    if (drawMagnitude && !m_drawableMagnitude)
        m_drawableMagnitude.create(context, size);

    if (drawDensity && !m_drawableDensity)
        m_drawableDensity.create(context, size);

    const RenderState& state = m_drawableState.getFront();

    if (drawDynamicsType && m_drawableDynamicsType)
        m_drawableDynamicsType->setImage(state.imageDynamicsType);

    if (drawMagnitude && m_drawableMagnitude)
        m_drawableMagnitude->setImage(state.imageMagnitude);

    if (drawDensity && m_drawableDensity)
        m_drawableDensity->setImage(state.imageDensity);

    // Draw color grid
    context.matrixPush();
    context.matrixScale(state.scale);

    if (drawDynamicsType && m_drawableDynamicsType)
        m_drawableDynamicsType->draw(context);

    if (drawMagnitude && m_drawableMagnitude)
        m_drawableMagnitude->draw(context);

    if (drawDensity && m_drawableDensity)
        m_drawableDensity->draw(context);

    context.matrixPop();
}
#endif

/* ************************************************************************ */

#ifdef CECE_RENDER
void Module::drawStoreState(const simulator::Visualization& visualization)
{
    const auto size = m_lattice.getSize();
    const bool drawDynamicsType = visualization.isEnabled(m_visualizationLayerDynamicsType);
    const bool drawMagnitude = visualization.isEnabled(m_visualizationLayerMagnitude);
    const bool drawDensity = visualization.isEnabled(m_visualizationLayerDensity);

    // Render state
    RenderState& state = m_drawableState.getBack();

    state.scale = getSimulation().getWorldSize() / units::Length(1);

    // Resize image
    state.imageDynamicsType.resize(size);
    state.imageMagnitude.resize(size);
    state.imageDensity.resize(size);

    // Find min/max density
    Descriptor::DensityType rhoMin = std::numeric_limits<Descriptor::DensityType>::max();
    Descriptor::DensityType rhoMax = 0.0;
    RealType maxVel = 0.0;

    for (auto&& c : range(size))
    {
        const auto& node = m_lattice[c];
        const auto dynamics = node.getDynamics();

        if (dynamics == getFluidDynamics())
        {
            const auto velocity = node.computeVelocity();
            const auto density = node.computeDensity();

            maxVel = std::max(maxVel, velocity.getLength());
            rhoMin = std::min(density, rhoMin);
            rhoMax = std::max(density, rhoMax);
        }
    }

    // Update texture
    for (auto&& c : range(size))
    {
        // Cell alias
        const auto& node = m_lattice[c];
        const auto velocity = node.computeVelocity();
        const auto density = node.computeDensity();
        const auto dynamics = node.getDynamics();

        // Store dynamics type
        if (drawDynamicsType)
        {
            render::Color color;

            if (dynamics == getFluidDynamics())
            {
                color = render::colors::BLACK;
            }
            else if (m_boundaries.isBoundaryDynamics(dynamics, Boundary::Type::Inlet))
            {
                color = render::colors::RED;
            }
            else if (m_boundaries.isBoundaryDynamics(dynamics, Boundary::Type::Outlet))
            {
                color = render::colors::BLUE;
            }
            else if (dynamics == NoDynamics::getInstance())
            {
                color = render::colors::WHITE;
            }
            else
            {
                color = render::colors::GREEN;
            }

            // Store color
            state.imageDynamicsType.set(c, color);
        }

        if (drawMagnitude)
        {
            if (dynamics == getFluidDynamics())
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
            if (dynamics == getFluidDynamics())
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
#endif

/* ************************************************************************ */

#ifdef CECE_RENDER
void Module::drawSwapState()
{
    m_drawableState.swap();
}
#endif

/* ************************************************************************ */

UniquePtr<Dynamics> Module::createFluidDynamics() const
{
    return makeUnique<BgkDynamics>(m_converter.calculateOmega());
}

/* ************************************************************************ */

UniquePtr<Dynamics> Module::createWallDynamics() const
{
    return makeUnique<BounceBackDynamics>();
}

/* ************************************************************************ */

UniquePtr<Dynamics> Module::createBorderDynamics(Boundary::Position pos) const
{
    const auto omega = m_converter.calculateOmega();

    switch (pos)
    {
    case Boundary::Position::Top:
        return makeUnique<ZouHeDynamics>(omega, ZouHeDynamics::Position::Top);

    case Boundary::Position::Bottom:
        return makeUnique<ZouHeDynamics>(omega, ZouHeDynamics::Position::Bottom);

    case Boundary::Position::Left:
        return makeUnique<ZouHeDynamics>(omega, ZouHeDynamics::Position::Left);

    case Boundary::Position::Right:
        return makeUnique<ZouHeDynamics>(omega, ZouHeDynamics::Position::Right);

    default:
        throw InvalidArgumentException("Unknown boundary type");
    }

    return nullptr;
}

/* ************************************************************************ */

void Module::updateObstacleMap()
{
    // Clear previous flag
    m_lattice.setDynamics(getFluidDynamics());

    const units::PositionVector start = getSimulation().getWorldSize() * -0.5f;
    const auto step = getSimulation().getWorldSize() / m_lattice.getSize();

    // Clear pool
    MovingDynamics::poolClear();

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
                shp.getCircle().radius *= m_circleObstacleScale;

            mapShapeToGrid(
                [&, this] (Coordinate&& coord) {
                    Assert(m_lattice.inRange(coord));
                    if (isDynamic && this->isDynamic())
                    {
                        m_lattice[coord].setDynamics(MovingDynamics::poolCreate(velocity));
                        movingObstacleMap[coord] = true;
                        m_movingObstacleMap[coord] = false;
                    }
                    else
                    {
                        m_lattice[coord].setDynamics(getWallDynamics());
                    }
                },
                [] (Coordinate&& coord) {},
                shp, step, center, obj->getRotation(), m_lattice.getSize()
            );
        }
    }

    // Only for dynamic obstacles
    if (isDynamicObjectsObstacles())
    {
        // Swap maps
        std::swap(m_movingObstacleMap, movingObstacleMap);

        // Nodes with 'true' in movingObstacleMap have changed from obstacle to fluid
        for (const auto& c : range(m_lattice.getSize()))
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
                if (!m_lattice.inRange(c2))
                    continue;

                if (m_lattice[c2].getDynamics() != getFluidDynamics())
                    continue;

                ++count;
                density += m_lattice[c2].computeDensity();
                velocity += m_lattice[c2].computeVelocity();
            }

            if (count)
            {
                m_lattice[c].initEquilibrium(velocity / count, density / count);
            }
            else
            {
                m_lattice[c].initEquilibrium();
            }
        }
    }

    m_lattice.fixupObstacles(getWallDynamics());

    // Update boundaries blocks
    m_boundaries.updateBlocks(m_lattice, m_converter, getFluidDynamics());

    //if (isDynamicObjectsObstacles())
    //    m_lattice.fixupObstacles(Node::Dynamics::DynamicObstacle);
}

/* ************************************************************************ */

void Module::applyToObjects()
{
    // Foreach objects
    for (auto& obj : getSimulation().getObjects())
        applyToObject(*obj);
}

/* ************************************************************************ */

void Module::applyToObject(object::Object& object)
{
    // Ignore static objects
    if (object.getType() != object::Object::Type::Dynamic)
        return;

    // Maximum object speed that is allowed by physical engine
    const auto maxSpeed = getSimulation().getMaxObjectTranslation() / getSimulation().getTimeStep();
    const auto maxSpeedSq = maxSpeed * maxSpeed;

    const units::PositionVector start = getSimulation().getWorldSize() * -0.5f;
    const auto step = getSimulation().getWorldSize() / m_lattice.getSize();

    // Get mass local center
    const auto center = object.getMassCenterOffset();

    // Coefficient used in force calculation
    const auto forceCoefficient = 6 * constants::PI * m_converter.getKinematicViscosity() * object.getDensity();

    if (isDynamicObjectsObstacles() && isDynamic())
    {
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

                // Get node
                const auto& node = m_lattice[coord];

                Vector<RealType> force = Zero;

                // Foreach all directions
                for (Descriptor::DirectionType iPop = 1; iPop < Descriptor::SIZE; ++iPop)
                {
                    const auto iOpp = Descriptor::DIRECTION_OPPOSITES[iPop];
                    const auto ei = Descriptor::DIRECTION_VELOCITIES[iOpp];

                    // Neighbour coordinate & node
                    const auto nc = coord + Coordinate(ei);

                    // Skip
                    if (!m_lattice.inRange(nc))
                        continue;

                    const auto& nn = m_lattice[nc];

                    // Not fluid dynamics
                    if (nn.getDynamics() != getFluidDynamics())
                        continue;

                    force += -ei * (nn[iOpp] + nn[iPop]);
                }

                forceLB += force;
            };

            auto shp = shape;

            if (shp.getType() == ShapeType::Circle)
                shp.getCircle().radius *= m_circleObstacleScale;

            // Store force for each coordinate
            mapShapeToGrid(
                forceEval,
                [] (Coordinate&& coord) { },
                shp, step, c, object.getRotation(), m_lattice.getSize()
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
        if (force.getLengthSquared() > m_maxForce.getLengthSquared())
        {
            const RealType ratio = m_maxForce.getLength() / force.getLength();
            force *= ratio;
        }

        // Calculate linear impulse from shapes
        const auto impulse = force * getSimulation().getTimeStep();

        // Apply impulse
        object.applyLinearImpulse(impulse);
    }
    else
    {
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
                    velocityLB += m_lattice[coord].computeVelocity();
                    ++count;
                },
                [] (Coordinate&& coord) { },
                shape, step, coord, m_lattice.getSize(), {}, 1
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
}

/* ************************************************************************ */

void Module::applyBoundaryConditions()
{
    m_boundaries.applyConditions(m_lattice, m_converter, getFluidDynamics());
}

/* ************************************************************************ */

void Module::printInfo()
{
    // Get values
    const auto size = m_lattice.getSize();

    Log::info("[streamlines] Viscosity: ", m_converter.getKinematicViscosity(), " um2/s");
    Log::info("[streamlines] Max object speed: ", getSimulation().getMaxObjectTranslation(), " um/it");
    Log::info("[streamlines] Char. length: ", m_converter.getCharLength(), " um");
    Log::info("[streamlines] Char. time: ", m_converter.getCharTime(), " s");
    Log::info("[streamlines] Char. speed: ", m_converter.getCharVelocity(), " um/s");
    Log::info("[streamlines] Number of nodes: ", m_converter.getNumberNodes());
    Log::info("[streamlines] Number of time steps: ", m_converter.getNumberSteps());
    Log::info("[streamlines] Re: ", m_converter.calculateRe());
    Log::info("[streamlines] ## Lattice ##");
    Log::info("[streamlines] Tau: ", m_converter.calculateTau());
    Log::info("[streamlines] Omega: ", m_converter.calculateOmega());
    Log::info("[streamlines] Grid: (", size.getWidth(), "; ", size.getHeight(), ")");
    Log::info("[streamlines] Viscosity: ", m_converter.calculateViscosity());
}

/* ************************************************************************ */

void Module::storeToFile(const FilePath& filename)
{
    OutFileStream ofs(filename.string(), OutFileStream::binary);
    BinaryOutput out(ofs);

    // Write file guard
    out.write(FILE_GUARD);

    // Write lattice size
    out.write(m_lattice.getSize());

    // Store lattice hash
    out.write(calculateLatticeHash());

    // Write relaxation time
    out.write(m_converter.calculateTau());

    // Number of init iterations
    out.write(m_initIterations);

    for (auto&& c : range(m_lattice.getSize()))
    {
        const Node& cell = m_lattice[c];

        // Write cell populations
        out.write(cell.getData());
    }
}

/* ************************************************************************ */

void Module::loadFromFile(const FilePath& filename)
{
    InFileStream ifs(filename.string(), OutFileStream::binary);

    if (!ifs.is_open())
        throw InvalidArgumentException("Cannot load from file: File not found '" + filename.string() + "'");

    BinaryInput in(ifs);

    // Read guard
    StaticArray<char, FILE_GUARD.size()> guard;
    in.read(guard);

    if (guard != FILE_GUARD)
        throw InvalidArgumentException("Cannot load from file: File is not valid");

    // Read lattice size
    Lattice::Size size;
    in.read(size);

    if (size != m_lattice.getSize())
        throw InvalidArgumentException("Cannot load from file: different lattice sizes");

    std::size_t hash;
    in.read(hash);

    if (hash != calculateLatticeHash())
        throw InvalidArgumentException("Cannot load from file: different layout");

    RealType tau;
    in.read(tau);

    if (tau != m_converter.calculateTau())
        throw InvalidArgumentException("Cannot load from file: different relaxation times");

    decltype(m_initIterations) iterations;
    in.read(iterations);

    if (iterations != m_initIterations)
        throw InvalidArgumentException("Cannot load from file: different init iterations");

    for (auto&& c : range(m_lattice.getSize()))
    {
        Node& cell = m_lattice[c];

        // Read cell populations
        in.read(cell.getData());
    }
}

/* ************************************************************************ */

std::size_t Module::calculateLatticeHash() const noexcept
{
    // Based on djb2
    std::size_t hash = 5381;
    const auto noDynamics = NoDynamics::getInstance();
    const auto wallDynamics = getWallDynamics();
    const auto fluidDynamics = getFluidDynamics();

    for (auto&& c : range(m_lattice.getSize()))
    {
        const Node& node = m_lattice[c];
        const auto dynamics = node.getDynamics();

        std::size_t val;

        // Get value based on dynamics
        if (dynamics == noDynamics)
            val = 1;
        else if (dynamics == wallDynamics)
            val = 2;
        else if (dynamics == fluidDynamics)
            val = 3;
        else
            val = 4;

        // Can overflow... don't care
        hash = ((hash << 5) + hash) + val; // hash * 33 + val
    }

    return hash;
}

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
