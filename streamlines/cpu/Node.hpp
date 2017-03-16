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
#include "cece/core/ViewPtr.hpp"

// Plugin
#include "../Descriptor.hpp"
#include "../Dynamics.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines {
namespace cpu {

/* ************************************************************************ */

/**
 * @brief      Class for storing lattice data and perform LB dynamics.
 */
class Node
{

// Public Types
public:

    /// Population index type.
    using PopIndexType = typename Descriptor::PopIndexType;

    /// Density type.
    using DensityType = typename Descriptor::DensityType;

    /// Momentum type.
    using MomentumType = typename Descriptor::MomentumType;

    /// Velocity type.
    using VelocityType = typename Descriptor::VelocityType;

    /// Distribution function type.
    using DistributionType = typename Descriptor::DistributionType;

    /// Distributions type.
    using DistributionsType = typename Descriptor::DistributionsType;


// Public Enums
public:


    /**
     * @brief      Boundary position.
     */
    enum class BoundaryPosition
    {
        Right  = 0,
        Left   = 1,
        Top    = 2,
        Bottom = 3
    };


// Public Operators
public:


    /**
     * @brief      Returns microscopic density population in required direction.
     *
     * @param      iPop  Direction index.
     *
     * @return     Distribution function.
     */
    DistributionType operator[](PopIndexType iPop) const noexcept;


    /**
     * @brief      Returns microscopic density population in required direction.
     *
     * @param      iPop  Direction index.
     *
     * @return     Distribution function.
     */
    DistributionType& operator[](PopIndexType iPop) noexcept;


// Public Accessors & Mutators
public:


    /**
     * @brief      Returns node dynamics.
     *
     * @return     The dynamics.
     */
    Dynamics getDynamics() const noexcept;


    /**
     * @brief      Set none dynamics at given coordinate.
     */
    void setNoneDynamics() noexcept;


    /**
     * @brief      Set fluid dynamics at given coordinate.
     */
    void setFluidDynamics() noexcept;


    /**
     * @brief      Set wall dynamics at given coordinate.
     */
    void setWallDynamics() noexcept;


    /**
     * @brief      Set inlet dynamics at given coordinate.
     *
     * @param[in]  position  The boundary position.
     * @param[in]  velocity  The inlet velocity.
     */
    void setInletDynamics(BoundaryPosition position, VelocityType velocity) noexcept;


    /**
     * @brief      Set outlet dynamics at given coordinate.
     *
     * @param[in]  position  The boundary position.
     * @param[in]  density  The outlet density.
     */
    void setOutletDynamics(BoundaryPosition position, DensityType density) noexcept;


    /**
     * @brief      Set an dynamics at given coordinate object.
     *
     * @param[in]  velocity  The object velocity.
     */
    void setObjectDynamics(VelocityType velocity) noexcept;


    /**
     * @brief      Returns microscopic density population in required direction.
     *
     * @param      iPop  Direction index.
     *
     * @return     Distribution function.
     */
    DistributionType& get(PopIndexType iPop) noexcept;


    /**
     * @brief      Returns microscopic density population in required direction.
     *
     * @param      iPop  Direction index.
     *
     * @return     Distribution function.
     */
    DistributionType get(PopIndexType iPop) const noexcept;


    /**
     * @brief      Returns node distribution functions.
     *
     * @return     The distributions.
     */
    DistributionsType& getDistributions() noexcept;


    /**
     * @brief      Returns node distribution functions.
     *
     * @return     The distributions.
     */
    const DistributionsType& getDistributions() const noexcept;


    /**
     * @brief      Set distribution functions data.
     *
     * @param[in]  distributions  The distribution functions.
     */
    void setDistributions(DistributionsType distributions) noexcept;


    /**
     * @brief      Define node's velocity.
     *
     * @param      velocity  Required macroscopic velocity.
     */
    void defineVelocity(const VelocityType& velocity) noexcept;


    /**
     * @brief      Define node's density.
     *
     * @param      density  Required macroscopic density.
     * @param      context  Computation context.
     */
    void defineDensity(DensityType density) noexcept;


// Public Operations
public:


    /**
     * @brief      Compute macroscopic density.
     *
     * @return     The node density.
     */
    DensityType computeDensity() const noexcept;


    /**
     * @brief      Compute macroscopic momentum.
     *
     * @return     The node momentum.
     */
    MomentumType computeMomentum() const noexcept;


    /**
     * @brief      Compute macroscopic velocity.
     *
     * @return     The node velocity.
     */
    VelocityType computeVelocity() const noexcept;


    /**
     * @brief Initialize node by equilibrium populations.
     *
     * @param velocity Cell velocity vector.
     * @param density  Cell density.
     */
    void initEquilibrium(VelocityType velocity = Zero, DensityType density = Descriptor::DEFAULT_DENSITY) noexcept;


    /**
     * @brief      Perform node's populations collision.
     *
     * @param[in]  omega  The relaxation frequency.
     */
    void collide(RealType omega) noexcept;


// Private Operations
private:


    /**
     * @brief      Perform node's populations collision for fluid dynamics.
     *
     * @param[in]  omega  The relaxation frequency.
     */
    void collideFluid(RealType omega) noexcept;


    /**
     * @brief      Perform node's populations collision for wall dynamics.
     */
    void collideWall() noexcept;


    /**
     * @brief      Perform node's populations collision for object dynamics.
     *
     * @param[in]  omega  The relaxation frequency.
     */
    void collideObject(RealType omega) noexcept;


    /**
     * @brief      Fill missing distributions at inlet boundary.
     */
    void boundaryInlet() noexcept;


    /**
     * @brief      Fill missing distributions at outlet boundary.
     */
    void boundaryOutlet() noexcept;


// Private Structures
private:


    /// Flow dynamics data
    struct FlowData { };

    /// Wall dynamics data
    struct WallData { };

    /// Inlet dynamics data
    struct InletData
    {
        BoundaryPosition position;
        VelocityType velocity;
    };

    /// Outlet dynamics data
    struct OutletData
    {
        BoundaryPosition position;
        DensityType density;
    };

    /// Object dynamics data
    struct ObjectData
    {
        VelocityType velocity;
    };


// Private Data Members
private:

    /// Node dynamics type.
    Dynamics m_dynamics = Dynamics::None;

    /// Distribution functions.
    DistributionsType m_df;

    /// Additional data
    struct
    {
        FlowData flow;
        WallData wall;
        InletData inlet;
        OutletData outlet;
        ObjectData object;
    } m_data;
};

/* ************************************************************************ */

}
}
}
}

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines {
namespace cpu {

/* ************************************************************************ */

inline Node::DistributionType Node::operator[](PopIndexType iPop) const noexcept
{
    return get(iPop);
}

/* ************************************************************************ */

inline Node::DistributionType& Node::operator[](PopIndexType iPop) noexcept
{
    return get(iPop);
}

/* ************************************************************************ */

inline Dynamics Node::getDynamics() const noexcept
{
    return m_dynamics;
}

/* ************************************************************************ */

inline Node::DistributionType& Node::get(PopIndexType iPop) noexcept
{
    return m_df[iPop];
}

/* ************************************************************************ */

inline Node::DistributionType Node::get(PopIndexType iPop) const noexcept
{
    return m_df[iPop];
}

/* ************************************************************************ */

inline Node::DistributionsType& Node::getDistributions() noexcept
{
    return m_df;
}

/* ************************************************************************ */

inline const Node::DistributionsType& Node::getDistributions() const noexcept
{
    return m_df;
}

/* ************************************************************************ */

inline void Node::setDistributions(DistributionsType distributions) noexcept
{
    m_df = distributions;
}

/* ************************************************************************ */

}
}
}
}

/* ************************************************************************ */
