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
#include "cece/core/StringView.hpp"
#include "cece/core/ViewPtr.hpp"
#include "cece/core/DynamicArray.hpp"

// Plugin
#include "Boundary.hpp"

/* ************************************************************************ */

namespace cece { namespace config { class Configuration; } }

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines {

/* ************************************************************************ */

/**
 * @brief      A set of streamlines boundaries.
 */
class Boundaries
{

// Public Operators
public:


    /**
     * @brief      Get boundary by position.
     *
     * @param      position  The position
     *
     * @return     Boundary at position.
     */
    Boundary& operator[](int position) noexcept;


    /**
     * @brief      Get boundary by position.
     *
     * @param      position  The position.
     *
     * @return     Boundary at position.
     */
    const Boundary& operator[](int position) const noexcept;


// Public Accessors
public:


    /**
     * @brief      Returns a number of boundaries.
     *
     * @return     The number of boundaries.
     */
    int getCount() const noexcept;


    /**
     * @brief      Get boundary by position.
     *
     * @param      position  The position.
     *
     * @return     Reference to boundary.
     */
    Boundary& get(int position) noexcept;


    /**
     * @brief      Get boundary by position.
     *
     * @param      position  The position.
     *
     * @return     Reference to boundary.
     */
    const Boundary& get(int position) const noexcept;


    /**
     * @brief      Find boundary by name.
     *
     * @param      name  Boundary name.
     *
     * @return     Pointer to boundary or nullptr.
     */
    ViewPtr<Boundary> find(StringView name) noexcept;


    /**
     * @brief      Find boundary by name.
     *
     * @param      name  Boundary name.
     *
     * @return     Pointer to boundary or nullptr.
     */
    ViewPtr<const Boundary> find(StringView name) const noexcept;


    /**
     * @brief      Returns begin iterator.
     *
     * @return     Begin iterator.
     */
    DynamicArray<Boundary>::iterator begin() noexcept;


    /**
     * @brief      Returns begin iterator.
     *
     * @return     Begin iterator.
     */
    DynamicArray<Boundary>::const_iterator begin() const noexcept;


    /**
     * @brief      Returns begin iterator.
     *
     * @return     Begin iterator.
     */
    DynamicArray<Boundary>::const_iterator cbegin() const noexcept;


    /**
     * @brief      Returns end iterator.
     *
     * @return     End iterator.
     */
    DynamicArray<Boundary>::iterator end() noexcept;


    /**
     * @brief      Returns end iterator.
     *
     * @return     End iterator.
     */
    DynamicArray<Boundary>::const_iterator end() const noexcept;


    /**
     * @brief      Returns end iterator.
     *
     * @return     End iterator.
     */
    DynamicArray<Boundary>::const_iterator cend() const noexcept;


// Public Operations
public:


    /**
     * @brief      Load module configuration.
     *
     * @param      config  Source configuration.
     */
    void loadConfig(const config::Configuration& config);


    /**
     * @brief      Store module configuration.
     *
     * @param      config  Output configuration.
     */
    void storeConfig(config::Configuration& config) const;


// Private Operations
private:


    /**
     * @brief      Initialize default boundaries.
     */
    void initDefault();


// Private Data Members
private:

    /// Stored boundaries.
    DynamicArray<Boundary> m_boundaries;

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

inline Boundary& Boundaries::operator[](int position) noexcept
{
    return get(position);
}

/* ************************************************************************ */

inline const Boundary& Boundaries::operator[](int position) const noexcept
{
    return get(position);
}

/* ************************************************************************ */

inline int Boundaries::getCount() const noexcept
{
    return m_boundaries.size();
}

/* ************************************************************************ */

inline Boundary& Boundaries::get(int position) noexcept
{
    return m_boundaries[position];
}

/* ************************************************************************ */

inline const Boundary& Boundaries::get(int position) const noexcept
{
    return m_boundaries[position];
}

/* ************************************************************************ */

inline DynamicArray<Boundary>::iterator Boundaries::begin() noexcept
{
    return m_boundaries.begin();
}

/* ************************************************************************ */

inline DynamicArray<Boundary>::const_iterator Boundaries::begin() const noexcept
{
    return m_boundaries.begin();
}

/* ************************************************************************ */

inline DynamicArray<Boundary>::const_iterator Boundaries::cbegin() const noexcept
{
    return m_boundaries.cbegin();
}

/* ************************************************************************ */

inline DynamicArray<Boundary>::iterator Boundaries::end() noexcept
{
    return m_boundaries.end();
}

/* ************************************************************************ */

inline DynamicArray<Boundary>::const_iterator Boundaries::end() const noexcept
{
    return m_boundaries.end();
}

/* ************************************************************************ */

inline DynamicArray<Boundary>::const_iterator Boundaries::cend() const noexcept
{
    return m_boundaries.cend();
}

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
