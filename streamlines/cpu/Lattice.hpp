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
#include "cece/core/Grid.hpp"

// Plugin
#include "../Lattice.hpp"
#include "Node.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines {
namespace cpu {

/* ************************************************************************ */

/**
 * @brief      CPU implementation of Lattice Boltzman.
 *
 * @note       Memory and speedup improvement taken from OpenLB TR1:
 * @link http://optilb.com/openlb/wp-content/uploads/2011/12/olb-tr1.pdf
 */
class Lattice : public streamlines::Lattice
{

// Public Ctors & Dtors
public:


    /**
     * @brief      Constructor.
     *
     * @param[in]  size   The lattice size.
     * @param[in]  omega  The relaxation frequency.
     */
    Lattice(SizeType size, RealType omega);


// Public Operators
public:


    /**
     * @brief      Node access operator.
     *
     * @param      coord  The coordinate.
     *
     * @return     Node at given coordinate.
     */
    Node& operator[](const CoordinateType& coord) noexcept;


    /**
     * @brief      Node access operator.
     *
     * @param      coord  The coordinate.
     *
     * @return     Node at given coordinate.
     */
    const Node& operator[](const CoordinateType& coord) const noexcept;


// Public Accessors & Mutators
public:


    /**
     * @brief      Get node at given coordinate.
     *
     * @param      coord  The coordinate.
     *
     * @return     Node at given coordinate.
     */
    Node& get(const CoordinateType& coord) noexcept;


    /**
     * @brief      Get node at given coordinate.
     *
     * @param      coord  The coordinate.
     *
     * @return     Node at given coordinate.
     */
    const Node& get(const CoordinateType& coord) const noexcept;


    /**
     * @brief      Get dynamics at given coordinate.
     *
     * @param[in]  coord  The coordinate.
     *
     * @return     The dynamics.
     */
    Dynamics getDynamics(CoordinateType coord) const noexcept override;


    /**
     * @brief      Get velocity at given coordinate.
     *
     * @param[in]  coord  The coordinate.
     *
     * @return     The velocity.
     */
    VelocityType getVelocity(CoordinateType coord) const noexcept override;


    /**
     * @brief      Set velocity at given coordinate.
     *
     * @details    Function update node velocity using equilibrium distribution
     *             with node's density.
     *
     * @param[in]  coord     The coordinate.
     * @param[in]  velocity  The velocity.
     */
    void setVelocity(CoordinateType coord, VelocityType velocity) noexcept override;


    /**
     * @brief      Get density at given coordinate.
     *
     * @param[in]  coord  The coordinate.
     *
     * @return     The density.
     */
    DensityType getDensity(CoordinateType coord) const noexcept override;


    /**
     * @brief      Set density at given coordinate.
     *
     * @details    Function update node density using equilibrium distribution
     *             with node's velocity.
     *
     * @param[in]  coord    The coordinate.
     * @param[in]  density  The density.
     */
    void setDensity(CoordinateType coord, DensityType density) noexcept override;


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
    void setVelocityDensity(CoordinateType coord, VelocityType velocity, DensityType density) noexcept override;


    /**
     * @brief      Get distribution functions at given coordinate.
     *
     * @param[in]  coord     The coordinate.
     *
     * @return     The distribution functions.
     */
    DistributionsType getDistributions(CoordinateType coord) const noexcept override;


    /**
     * @brief      Set distribution functions at given coordinate.
     *
     * @param[in]  coord          The coordinate.
     * @param[in]  distributions  The distribution functions.
     */
    void setDistributions(CoordinateType coord, DistributionsType distributions) noexcept override;


    /**
     * @brief      Set none dynamics at given coordinate.
     *
     * @param[in]  coord  The coordinate.
     */
    void setNoneDynamics(CoordinateType coord) noexcept override;


    /**
     * @brief      Set fluid dynamics at given coordinate.
     *
     * @param[in]  coord  The coordinate.
     */
    void setFluidDynamics(CoordinateType coord) noexcept override;


    /**
     * @brief      Set wall dynamics at given coordinate.
     *
     * @param[in]  coord  The coordinate.
     */
    void setWallDynamics(CoordinateType coord) noexcept override;


    /**
     * @brief      Set inlet dynamics at given coordinate.
     *
     * @param[in]  coord     The coordinate.
     * @param[in]  velocity  The inlet velocity.
     */
    void setInletDynamics(CoordinateType coord, VelocityType velocity) noexcept override;


    /**
     * @brief      Set outlet dynamics at given coordinate.
     *
     * @param[in]  coord    The coordinate.
     * @param[in]  density  The outlet density.
     */
    void setOutletDynamics(CoordinateType coord, DensityType density) noexcept override;


    /**
     * @brief      Set an dynamics at given coordinate object.
     *
     * @param[in]  coord     The coordinate.
     * @param[in]  velocity  The object velocity.
     */
    void setObjectDynamics(CoordinateType coord, VelocityType velocity) noexcept override;


// Public Operations
public:


    /**
     * @brief      Initialize whole lattice with default values.
     *
     * @details    As default values it means zero velocity with default density
     *             (equilibrium).
     */
    void initDefault() override;


    /**
     * @brief      Update lattice according to set parameters and current state.
     *
     * @param[in]  count  The number of inner iterations.
     */
    void update(unsigned int count = 1) override;


// Private Operations
private:


    /**
     * @brief      Perform collision.
     */
    void collide();


    /**
     * @brief      Perform streaming.
     */
    void stream();


    /**
     * @brief      Perform collision and streaming.
     */
    void collideAndStream();


    /**
     * @brief      Remove unreachable dynamics.
     *
     * @details    Remove inner wall dynamics.
     *
     * @param[in]  dynamics  The dynamics to remove.
     */
    void removeUnreachableDynamics(Dynamics dynamics);


// Private Data Members
private:

    /// Current lattice data.
    core::Grid<Node> m_data;

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

inline Node& Lattice::operator[](const CoordinateType& coord) noexcept
{
    return get(coord);
}

/* ************************************************************************ */

inline const Node& Lattice::operator[](const CoordinateType& coord) const noexcept
{
    return get(coord);
}

/* ************************************************************************ */

inline Node& Lattice::get(const CoordinateType& coord) noexcept
{
    return m_data[coord];
}

/* ************************************************************************ */

inline const Node& Lattice::get(const CoordinateType& coord) const noexcept
{
    return m_data[coord];
}

/* ************************************************************************ */

}
}
}
}

/* ************************************************************************ */
