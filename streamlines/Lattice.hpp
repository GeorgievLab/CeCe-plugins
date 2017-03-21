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
#include "cece/core/Vector.hpp"
#include "cece/core/FilePath.hpp"

// Plugin
#include "Dynamics.hpp"
#include "Descriptor.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines {

/* ************************************************************************ */

/**
 * @brief      Interface for Lattice Boltzmann implementation.
 */
class Lattice
{

// Public Types
public:

    /// Density type.
    using DensityType = Descriptor::DensityType;

    /// Momentum type.
    using MomentumType = Descriptor::MomentumType;

    /// Velocity type.
    using VelocityType = Descriptor::VelocityType;

    /// Type of distribution functions
    using DistributionsType = Descriptor::DistributionsType;

    /// Coordinate type.
    using CoordinateType = Vector<unsigned int>;

    /// Size type.
    using SizeType = Vector<unsigned int>;


// Public Ctors & Dtors
public:


    /**
     * @brief      Constructor.
     *
     * @param[in]  size   The lattice size.
     * @param[in]  omega  The relaxation frequency.
     */
    Lattice(SizeType size, RealType omega);


    /**
     * @brief      Destructor.
     */
    virtual ~Lattice() = 0;


// Public Accessors & Mutators
public:


    /**
     * @brief      Returns lattice size.
     *
     * @return     The lattice size.
     */
    SizeType getSize() const noexcept;


    /**
     * @brief      Returns lattice relaxation frequency.
     *
     * @return     The lattice relaxation frequency.
     */
    RealType getOmega() const noexcept;


    /**
     * @brief      Check if coordinates are in range.
     *
     * @param[in]  coord  The coordinate.
     *
     * @return     Coordinate is in range.
     */
    bool inRange(CoordinateType coord) const noexcept;


    /**
     * @brief      Get dynamics at given coordinate.
     *
     * @param[in]  coord  The coordinate.
     *
     * @return     The dynamics.
     */
    virtual Dynamics getDynamics(CoordinateType coord) const = 0;


    /**
     * @brief      Get velocity at given coordinate.
     *
     * @param[in]  coord  The coordinate.
     *
     * @return     The velocity.
     */
    virtual VelocityType getVelocity(CoordinateType coord) const = 0;


    /**
     * @brief      Set velocity at given coordinate.
     *
     * @details    Function update node velocity using equilibrium distribution
     *             with node's density.
     *
     * @param[in]  coord     The coordinate.
     * @param[in]  velocity  The velocity.
     */
    virtual void setVelocity(CoordinateType coord, VelocityType velocity) = 0;


    /**
     * @brief      Get density at given coordinate.
     *
     * @param[in]  coord  The coordinate.
     *
     * @return     The density.
     */
    virtual DensityType getDensity(CoordinateType coord) const = 0;


    /**
     * @brief      Set density at given coordinate.
     *
     * @details    Function update node density using equilibrium distribution
     *             with node's velocity.
     *
     * @param[in]  coord    The coordinate.
     * @param[in]  density  The density.
     */
    virtual void setDensity(CoordinateType coord, DensityType density) = 0;


    /**
     * @brief      Set velocity and density at given coordinate.
     *
     * @details    Function update node velocity and density using equilibrium
     *             distribution.
     *
     * @param[in]  coord     The coordinate.
     * @param[in]  velocity  The velocity.
     * @param[in]  density   The density.
     */
    virtual void setVelocityDensity(CoordinateType coord, VelocityType velocity, DensityType density) = 0;


    /**
     * @brief      Get distribution functions at given coordinate.
     *
     * @param[in]  coord     The coordinate.
     *
     * @return     The distribution functions.
     */
    virtual DistributionsType getDistributions(CoordinateType coord) const = 0;


    /**
     * @brief      Set distribution functions at given coordinate.
     *
     * @param[in]  coord          The coordinate.
     * @param[in]  distributions  The distribution functions.
     */
    virtual void setDistributions(CoordinateType coord, DistributionsType distributions) = 0;


    /**
     * @brief      Determines if dynamics at given coordinate is none.
     *
     * @param[in]  coord  The coordinate.
     *
     * @return     True if is fluid, False otherwise.
     */
    bool isNoneDynamics(CoordinateType coord) const;


    /**
     * @brief      Determines if dynamics at given coordinate is fluid.
     *
     * @param[in]  coord  The coordinate.
     *
     * @return     True if is fluid, False otherwise.
     */
    bool isFluidDynamics(CoordinateType coord) const;


    /**
     * @brief      Determines if dynamics at given coordinate is wall.
     *
     * @param[in]  coord  The coordinate.
     *
     * @return     True if is wall, False otherwise.
     */
    bool isWallDynamics(CoordinateType coord) const;


    /**
     * @brief      Determines if dynamics at given coordinate is inlet.
     *
     * @param[in]  coord  The coordinate.
     *
     * @return     True if is inlet, False otherwise.
     */
    bool isInletDynamics(CoordinateType coord) const;


    /**
     * @brief      Determines if dynamics at given coordinate is outlet.
     *
     * @param[in]  coord  The coordinate.
     *
     * @return     True if is outlet, False otherwise.
     */
    bool isOutletDynamics(CoordinateType coord) const;


    /**
     * @brief      Determines if dynamics at given coordinate is an object.
     *
     * @param[in]  coord  The coordinate.
     *
     * @return     True if is an object, False otherwise.
     */
    bool isObjectDynamics(CoordinateType coord) const;


    /**
     * @brief      Set none dynamics at given coordinate.
     *
     * @param[in]  coord  The coordinate.
     */
    virtual void setNoneDynamics(CoordinateType coord) = 0;


    /**
     * @brief      Set fluid dynamics at given coordinate.
     *
     * @param[in]  coord  The coordinate.
     */
    virtual void setFluidDynamics(CoordinateType coord) = 0;


    /**
     * @brief      Set wall dynamics at given coordinate.
     *
     * @param[in]  coord  The coordinate.
     */
    virtual void setWallDynamics(CoordinateType coord) = 0;


    /**
     * @brief      Set inlet dynamics at given coordinate.
     *
     * @param[in]  coord     The coordinate.
     * @param[in]  velocity  The inlet velocity.
     */
    virtual void setInletDynamics(CoordinateType coord, VelocityType velocity) = 0;


    /**
     * @brief      Set outlet dynamics at given coordinate.
     *
     * @param[in]  coord    The coordinate.
     * @param[in]  density  The outlet density.
     */
    virtual void setOutletDynamics(CoordinateType coord, DensityType density) = 0;


    /**
     * @brief      Set an dynamics at given coordinate object.
     *
     * @param[in]  coord     The coordinate.
     * @param[in]  velocity  The object velocity.
     */
    virtual void setObjectDynamics(CoordinateType coord, VelocityType velocity) = 0;


// Public Operations
public:


    /**
     * @brief      Initialize whole lattice with default values.
     *
     * @details    As default values it means zero velocity with default density
     *             (equilibrium).
     */
    virtual void initDefault() = 0;


    /**
     * @brief      Update lattice according to set parameters and current state.
     *
     * @param[in]  count  The number of inner iterations.
     */
    virtual void update(unsigned int count = 1) = 0;


// Private Data Members
public:

    /// Lattice size.
    SizeType m_size;

    /// Fluid relaxation frequency.
    RealType m_omega;
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

inline Lattice::Lattice(SizeType size, RealType omega)
    : m_size(size)
    , m_omega(omega)
{
    // Nothing to do
}

/* ************************************************************************ */

inline Lattice::SizeType Lattice::getSize() const noexcept
{
    return m_size;
}

/* ************************************************************************ */

inline RealType Lattice::getOmega() const noexcept
{
    return m_omega;
}

/* ************************************************************************ */

inline bool Lattice::inRange(CoordinateType coord) const noexcept
{
    return coord.inRange(Zero, m_size);
}

/* ************************************************************************ */

inline bool Lattice::isNoneDynamics(CoordinateType coord) const
{
    return getDynamics(coord) == Dynamics::None;
}

/* ************************************************************************ */

inline bool Lattice::isFluidDynamics(CoordinateType coord) const
{
    return getDynamics(coord) == Dynamics::Fluid;
}

/* ************************************************************************ */

inline bool Lattice::isWallDynamics(CoordinateType coord) const
{
    return getDynamics(coord) == Dynamics::Wall;
}

/* ************************************************************************ */

inline bool Lattice::isInletDynamics(CoordinateType coord) const
{
    return getDynamics(coord) == Dynamics::Inlet;
}

/* ************************************************************************ */

inline bool Lattice::isOutletDynamics(CoordinateType coord) const
{
    return getDynamics(coord) == Dynamics::Outlet;
}

/* ************************************************************************ */

inline bool Lattice::isObjectDynamics(CoordinateType coord) const
{
    return getDynamics(coord) == Dynamics::Object;
}

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
