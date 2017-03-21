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

#pragma once

/* ************************************************************************ */

// CeCe
#include "cece/config.hpp"
#include "cece/core/Units.hpp"
#include "cece/core/UnitsCtors.hpp"
#include "cece/core/Vector.hpp"
#include "cece/core/Pair.hpp"
#include "cece/core/Units.hpp"
#include "cece/core/Grid.hpp"
#include "cece/core/Real.hpp"
#include "cece/core/StaticArray.hpp"
#include "cece/core/UniquePtr.hpp"
#include "cece/core/OutStream.hpp"
#include "cece/core/InStream.hpp"
#include "cece/core/FilePath.hpp"
#include "cece/core/IterationType.hpp"
#include "cece/module/Module.hpp"
#include "cece/object/Object.hpp"

#ifdef CECE_RENDER
#  include "cece/render/Context.hpp"
#  include "cece/render/Object.hpp"
#  include "cece/render/State.hpp"
#  include "cece/render/Image.hpp"
#  include "cece/render/GridColor.hpp"
#  include "cece/render/GridColorColorMap.hpp"
#endif

// Plugin
#include "Boundary.hpp"
#include "Lattice.hpp"
#include "Converter.hpp"
#include "Descriptor.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines {

/* ************************************************************************ */

/**
 * @brief      Pressure unit.
 */
using Pressure = units::Unit<
    units::List<units::BaseMass>,
    units::List<units::BaseLength, units::BaseTime, units::BaseTime>
>;

/* ************************************************************************ */

/**
 * @brief      Module for streamlines.
 */
class Module : public module::Module
{

// Public Ctors & Dtors
public:


    /**
     * @brief      Constructor.
     *
     * @param      simulation  The simulation.
     */
    Module(simulator::Simulation& simulation);


    /**
     * @brief      Destructor.
     */
    virtual ~Module();


// Public Accessors & Mutators
public:


    /**
     * @brief      Returns units converter.
     *
     * @return     The units converter.
     */
    const Converter& getConverter() const noexcept;


    /**
     * @brief      Returns units converter.
     *
     * @return     The units converter.
     */
    Converter& getConverter() noexcept;


    /**
     * @brief      Returns number of init iterations.
     *
     * @return     Init iterations.
     */
    IterationType getInitIterations() const noexcept;


    /**
     * @brief      Set number of init iterations.
     *
     * @param      iterations  The number of init iterations.
     */
    void setInitIterations(IterationType iterations) noexcept;


    /**
     * @brief      Returns init iteration count.
     *
     * @return     Init iterations.
     */
    IterationType getInnerIterations() const noexcept;


    /**
     * @brief      Set number of inner iterations.
     *
     * @param      iterations  The number of iterations.
     */
    void setInnerIterations(IterationType iterations) noexcept;


    /**
     * @brief      Returns boundaries.
     *
     * @return     The boundaries.
     */
    const DynamicArray<Boundary>& getBoundaries() const noexcept;


    /**
     * @brief      Find boundary by name.
     *
     * @param[in]  name  The boundary name.
     *
     * @return     Pointer to boundary or nullptr.
     */
    ViewPtr<const Boundary> findBoundary(StringView name) const noexcept;


    /**
     * @brief      Obtain boundary blocks where real inlet or outlet is.
     *
     * @param[in]  name  The boundary name.
     *
     * @return     A list of blocks.
     */
    DynamicArray<Pair<units::PositionVector, units::PositionVector>> getBoundaryBlocks(StringView name) const noexcept;


    /**
     * @brief      If dynamic objects are used as obstacles.
     *
     * @return     True if dynamic objects obstacles, False otherwise.
     */
    bool isDynamicObjectsObstacles() const noexcept;


    /**
     * @brief      Enable or disable dynamic objects obstacles.
     *
     * @param      flag  The flag.
     */
    void setDynamicObjectsObstacles(bool flag) noexcept;


    /**
     * @brief      Returns if streamlines are dynamic during simulation.
     *
     * @return     True if dynamic, False otherwise.
     */
    bool isDynamic() const noexcept;


    /**
     * @brief      Set if streamlines should be dynamic.
     *
     * @param[in]  flag  If streamlines are dynamic.
     */
    void setDynamic(bool flag);


    /**
     * @brief      Returns the lattice size.
     *
     * @return     The lattice size.
     */
    Lattice::SizeType getLatticeSize() const noexcept;


    /**
     * @brief      Check if given coordinate is in range.
     *
     * @param[in]  coord  The coordinate.
     *
     * @return     If coordinate is in range if lattice.
     */
    bool inLatticeRange(Coordinate coord) const noexcept;


    /**
     * @brief      Obtain physical velocity at given coordinate.
     *
     * @param[in]  coord  The coordinate.
     *
     * @return     The velocity.
     */
    units::VelocityVector getVelocity(Coordinate coord) const;


    /**
     * @brief      Change physical velocity at given coordinate.
     *
     * @param[in]  coord     The coordinate.
     * @param[in]  velocity  The velocity.
     */
    void setVelocity(Coordinate coord, units::VelocityVector velocity);


    /**
     * @brief      Obtain physical pressure at given coordinate.
     *
     * @param[in]  coord  The coordinate.
     *
     * @return     The pressure.
     */
    Pressure getPressure(Coordinate coord) const;


    /**
     * @brief      Obtain dynamics at given coordinate.
     *
     * @param[in]  coord  The coordinate.
     *
     * @return     The dynamics.
     */
    Dynamics getDynamics(Coordinate coord) const;


    /**
     * @brief      Distribution values at given coordinate.
     *
     * @param[in]  coord  The coordinate.
     *
     * @return     The distribution functions.
     */
    Lattice::DistributionsType getDistributions(Coordinate coord) const;


// Public Operations
public:


    /**
     * @brief      Initialize module.
     *
     * @param      flag  The termination variable.
     */
    void init(AtomicBool& flag) override;


    /**
     * @brief      Load module configuration.
     *
     * @param      config  Source configuration.
     */
    void loadConfig(const config::Configuration& config) override;


    /**
     * @brief      Store module configuration.
     *
     * @param      config  Output configuration.
     */
    void storeConfig(config::Configuration& config) const override;


    /**
     * @brief      Update module state.
     */
    void update() override;


#ifdef CECE_RENDER

    /**
     * @brief      Render module.
     *
     * @param      visualization  The visualization.
     * @param      context        Rendering context.
     */
    void draw(const simulator::Visualization& visualization, render::Context& context) override;


    /**
     * @brief      Store current state for drawing.
     *
     * @details    State should be stored in back state because the front state
     *             is used for rendering. Drawing state should contain data that
     *             can be modified during update() call and are used for
     *             rendering.
     *
     * @param      visualization  Visualization context.
     */
    void drawStoreState(const simulator::Visualization& visualization) override;


    /**
     * @brief      Swap render state.
     *
     * @details    Calling this function should be guarded by mutex for all
     *             modules to ensure all modules are in same render state.
     *             Function should be fast because otherwise it will block
     *             rendering.
     */
    void drawSwapState() override;

#endif


// Protected Operations
protected:


    /**
     * @brief      Update lattice dynamics.
     */
    void updateDynamics();


    /**
     * @brief      Apply streamlines to objects.
     */
    void updateObjects();


    /**
     * @brief      Apply streamlines to object.
     *
     * @param      object  The object.
     */
    void updateObject(object::Object& object);


    /**
     * @brief      Apply streamlines to object in static mode.
     *
     * @param      object  The object.
     */
    void updateObjectStatic(object::Object& object);


    /**
     * @brief      Apply streamlines to object in dynamic mode.
     *
     * @param      object  The object.
     */
    void updateObjectDynamic(object::Object& object);


    /**
     * @brief      Set boundaries to lattice.
     */
    void updateBoundaries();


    /**
     * @brief      Set boundary to lattice.
     *
     * @param[in]  boundary  The boundary.
     */
    void updateBoundary(const Boundary& boundary);


    /**
     * @brief      Store streamlines data to file.
     *
     * @param      filename  The destination filename.
     */
    void storeToFile(const FilePath& filename);


    /**
     * @brief      Load streamlines data from file.
     *
     * @param      filename  The source filename.
     */
    void loadFromFile(const FilePath& filename);


    /**
     * @brief      Calculate lattice hash.
     *
     * @return     The lattice hash.
     */
    std::size_t calculateLatticeHash() const noexcept;


    /**
     * @brief      Create lattice.
     *
     * @param[in]  size  The size.
     */
    void createLattice(Lattice::SizeType size);


// Private Enums
private:


    /**
     * @brief      Implementation mode.
     */
    enum class LatticeMode
    {
        GPU,
        CPU
    };


// Private Structures
private:

#ifdef CECE_RENDER
    struct RenderState
    {
        units::ScaleVector scale;
        render::Image imageDynamics;
        render::Image imageMagnitude;
        render::Image imageDensity;
    };
#endif

// Private Data Members
private:

    /// Units converter.
    Converter m_converter;

    /// Lattice implementation.
    UniquePtr<Lattice> m_lattice;

    struct
    {
        LatticeMode mode = LatticeMode::GPU;

        /// Number of init iterations.
        IterationType initIterations = 0;

        /// Number of inner iterations.
        IterationType innerIterations = 1;

        /// If streamlines is updated during simulation iterations.
        bool dynamic = true;

        /// Path to initialization file.
        FilePath initFile;

        /// Use dynamic objects as obstacles
        bool dynamicObjectsObstacles = false;

        /// Scaling for circle obstacle.
        RealType circleObstacleScale = 1.0;

        /// Boundaries.
        DynamicArray<Boundary> boundaries;

        /// Maximum force which can be applied to object.
        units::ForceVector maxForce = {units::N(1), units::N(1)};
    } m_config;

    /// Map of moving obstacles (Grid uses std::vector and bool is an issue).
    Grid<char> m_movingObstacleMap;

#ifdef CECE_RENDER
    struct
    {
        /// Name of layer for flow dynamics type visualization.
        String layerDynamics;

        /// Name of layer for velocity magnitude visualization.
        String layerMagnitude;

        /// Name of layer for density visualization.
        String layerDensity;

        /// Rendering grid with dynamics type.
        render::ObjectPtr<render::GridColor> dynamics;

        /// Rendering grid with magnitude.
        render::ObjectPtr<render::GridColorColorMap> magnitude;

        /// Rendering grid with density.
        render::ObjectPtr<render::GridColorColorMap> density;

        /// Render state.
        render::State<RenderState> state;
    } m_render;
#endif

};

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines {

/* ************************************************************************ */

inline const Converter& Module::getConverter() const noexcept
{
    return m_converter;
}

/* ************************************************************************ */

inline Converter& Module::getConverter() noexcept
{
    return m_converter;
}

/* ************************************************************************ */

inline IterationType Module::getInitIterations() const noexcept
{
    return m_config.initIterations;
}

/* ************************************************************************ */

inline void Module::setInitIterations(IterationType iterations) noexcept
{
    m_config.initIterations = iterations;
}

/* ************************************************************************ */

inline IterationType Module::getInnerIterations() const noexcept
{
    return m_config.innerIterations;
}

/* ************************************************************************ */

inline void Module::setInnerIterations(IterationType iterations) noexcept
{
    m_config.innerIterations = iterations;
}

/* ************************************************************************ */

inline const DynamicArray<Boundary>& Module::getBoundaries() const noexcept
{
    return m_config.boundaries;
}

/* ************************************************************************ */

inline bool Module::isDynamicObjectsObstacles() const noexcept
{
    return m_config.dynamicObjectsObstacles;
}

/* ************************************************************************ */

inline void Module::setDynamicObjectsObstacles(bool flag) noexcept
{
    m_config.dynamicObjectsObstacles = flag;
}

/* ************************************************************************ */

inline bool Module::isDynamic() const noexcept
{
    return m_config.dynamic;
}

/* ************************************************************************ */

inline void Module::setDynamic(bool flag)
{
    m_config.dynamic = flag;
}

/* ************************************************************************ */

inline bool Module::inLatticeRange(Coordinate coord) const noexcept
{
    CECE_ASSERT(m_lattice);
    return m_lattice->inRange(coord);
}

/* ************************************************************************ */

inline units::VelocityVector Module::getVelocity(Coordinate coord) const
{
    CECE_ASSERT(m_lattice);
    return m_converter.convertVelocity(m_lattice->getVelocity(coord));
}

/* ************************************************************************ */

inline void Module::setVelocity(Coordinate coord, units::VelocityVector velocity)
{
    CECE_ASSERT(m_lattice);
    m_lattice->setVelocity(coord, m_converter.convertVelocity(velocity));
}

/* ************************************************************************ */

inline Pressure Module::getPressure(Coordinate coord) const
{
    CECE_ASSERT(m_lattice);
    return Pressure(m_lattice->getDensity(coord));
}

/* ************************************************************************ */

inline Dynamics Module::getDynamics(Coordinate coord) const
{
    return m_lattice->getDynamics(coord);
}

/* ************************************************************************ */

inline Lattice::DistributionsType Module::getDistributions(Coordinate coord) const
{
    CECE_ASSERT(m_lattice);
    return m_lattice->getDistributions(coord);
}

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
