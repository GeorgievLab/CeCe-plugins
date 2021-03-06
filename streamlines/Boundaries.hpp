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

// C++
#include <cstddef>

// CeCe
#include "cece/core/Units.hpp"
#include "cece/core/StringView.hpp"
#include "cece/core/ViewPtr.hpp"
#include "cece/core/DynamicArray.hpp"

// Plugin
#include "Boundary.hpp"

/* ************************************************************************ */

namespace cece { namespace config { class Configuration; } }
namespace cece { namespace simulator { class Simulation; } }

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines {

/* ************************************************************************ */

/**
 * @brief A set of streamlines boundaries.
 */
class Boundaries
{


// Public Ctors & Dtors
public:


    /**
     * @brief Constructor.
     *
     * @param simulation
     */
    explicit Boundaries(simulator::Simulation& simulation);


// Public Operators
public:


    /**
     * @brief Get boundary by position.
     *
     * @param position
     *
     * @return
     */
    Boundary& operator[](std::size_t position) noexcept
    {
        return get(position);
    }


    /**
     * @brief Get boundary by position.
     *
     * @param position
     *
     * @return
     */
    const Boundary& operator[](std::size_t position) const noexcept
    {
        return get(position);
    }


// Public Accessors
public:


    /**
     * @brief Returns a number of boundaries.
     *
     * @return
     */
    int getCount() const noexcept
    {
        return m_boundaries.size();
    }


    /**
     * @brief Return simulation.
     *
     * @return
     */
    simulator::Simulation& getSimulation() const noexcept
    {
        return m_simulation;
    }


    /**
     * @brief Get boundary by position.
     *
     * @param position
     *
     * @return
     */
    Boundary& get(std::size_t position) noexcept
    {
        return m_boundaries[position];
    }


    /**
     * @brief Get boundary by position.
     *
     * @param position
     *
     * @return
     */
    const Boundary& get(std::size_t position) const noexcept
    {
        return m_boundaries[position];
    }


    /**
     * @brief Find boundary by name.
     *
     * @param name Boundary name.
     *
     * @return
     */
    ViewPtr<Boundary> find(StringView name) noexcept;


    /**
     * @brief Find boundary by name.
     *
     * @param name Boundary name.
     *
     * @return
     */
    ViewPtr<const Boundary> find(StringView name) const noexcept;


    /**
     * @brief Check if given dynamics is boundary dynamics.
     *
     * @param dynamics
     *
     * @return
     */
    bool isBoundaryDynamics(ViewPtr<Dynamics> dynamics) const noexcept;


    /**
     * @brief Check if given dynamics is boundary dynamics.
     *
     * @param dynamics
     * @param type
     *
     * @return
     */
    bool isBoundaryDynamics(ViewPtr<Dynamics> dynamics, Boundary::Type type) const noexcept;


    /**
     * @brief      Returns begin iterator.
     *
     * @return     Begin iterator.
     */
    DynamicArray<Boundary>::iterator begin() noexcept
    {
        return m_boundaries.begin();
    }


    /**
     * @brief      Returns begin iterator.
     *
     * @return     Begin iterator.
     */
    DynamicArray<Boundary>::const_iterator begin() const noexcept
    {
        return m_boundaries.begin();
    }


    /**
     * @brief      Returns begin iterator.
     *
     * @return     Begin iterator.
     */
    DynamicArray<Boundary>::const_iterator cbegin() const noexcept
    {
        return m_boundaries.cbegin();
    }


    /**
     * @brief      Returns end iterator.
     *
     * @return     End iterator.
     */
    DynamicArray<Boundary>::iterator end() noexcept
    {
        return m_boundaries.end();
    }


    /**
     * @brief      Returns end iterator.
     *
     * @return     End iterator.
     */
    DynamicArray<Boundary>::const_iterator end() const noexcept
    {
        return m_boundaries.end();
    }


    /**
     * @brief      Returns end iterator.
     *
     * @return     End iterator.
     */
    DynamicArray<Boundary>::const_iterator cend() const noexcept
    {
        return m_boundaries.cend();
    }


// Public Operations
public:


    /**
     * @brief Initialize default boundaries.
     */
    void initDefault();


    /**
     * @brief Initialize boundaries.
     * @param lattice
     * @param fluidDynamics
     */
    void init(Lattice& lattice, ViewPtr<Dynamics> fluidDynamics);


    /**
     * @brief Update conditions blocks.
     *
     * @param lattice
     * @param converter
     * @param fluidDynamics
     */
    void updateBlocks(Lattice& lattice, Converter& converter, ViewPtr<Dynamics> fluidDynamics);


    /**
     * @brief Apply boundary conditions.
     *
     * @param lattice
     * @param converter
     * @param fluidDynamics
     */
    void applyConditions(Lattice& lattice, Converter& converter, ViewPtr<Dynamics> fluidDynamics);


    /**
     * @brief Load module configuration.
     *
     * @param config Source configuration.
     */
    void loadConfig(const config::Configuration& config);


    /**
     * @brief Store module configuration.
     *
     * @param config Output configuration.
     */
    void storeConfig(config::Configuration& config) const;


// Private Data Members
private:

    /// Simulation reference.
    simulator::Simulation& m_simulation;

    /// Stored boundaries.
    DynamicArray<Boundary> m_boundaries;

};

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
