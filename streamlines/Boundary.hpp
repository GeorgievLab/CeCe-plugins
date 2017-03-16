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
#include "cece/core/Units.hpp"
#include "cece/core/VectorUnits.hpp"
#include "cece/core/ViewPtr.hpp"
#include "cece/core/String.hpp"
#include "cece/core/Range.hpp"
#include "cece/core/DynamicArray.hpp"

/* ************************************************************************ */

namespace cece { namespace config { class Configuration; } }

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines {

/* ************************************************************************ */

/**
 * @brief      Streamlines boundary specification.
 */
class Boundary
{

// Public Enums
public:


    /**
     * @brief      Boundary type.
     */
    enum class Type
    {
        None,
        Inlet,
        Outlet
    };


    /**
     * @brief      Boundary position.
     */
    enum class Position
    {
        Right  = 1,
        Left   = 3,
        Top    = 0,
        Bottom = 2
    };


    /**
     * @brief      Inlet velocity profile type.
     */
    enum class InletProfileType
    {
        Auto,
        Constant
    };


// Public Accessors
public:


    /**
     * @brief      Returns boundary name.
     *
     * @return     The name.
     */
    const String& getName() const noexcept;


    /**
     * @brief      Returns boundary type.
     *
     * @return     The type.
     */
    Type getType() const noexcept;


    /**
     * @brief      Returns boundary position.
     *
     * @return     The position.
     */
    Position getPosition() const noexcept;


    /**
     * @brief      Returns boundary offset.
     *
     * @return     The offset.
     */
    units::Length getOffset() const noexcept;


    /**
     * @brief      Returns boundary size.
     *
     * @return     The size.
     */
    units::Length getSize() const noexcept;


    /**
     * @brief      Returns boundary type.
     *
     * @return     The inlet profile type.
     */
    InletProfileType getInletProfileType() const noexcept;


    /**
     * @brief      Returns inlet velocity.
     *
     * @return     The inlet velocity.
     */
    units::Velocity getInletVelocity() const noexcept;


    /**
     * @brief      Returns inlet volumeric flow rate.
     *
     * @return     The inlet flow.
     */
    units::VolumericFlow getInletFlow() const noexcept;


// Public Mutators
public:


    /**
     * @brief      Set boundary name.
     *
     * @param      name  The name
     */
    void setName(String name) noexcept;


    /**
     * @brief      Set boundary type.
     *
     * @param      type  The type
     */
    void setType(Type type) noexcept;


    /**
     * @brief      Set boundary position.
     *
     * @param      position  The position
     */
    void setPosition(Position position) noexcept;


    /**
     * @brief      Set boundary offset.
     *
     * @param      offset  The offset
     */
    void setOffset(units::Length offset) noexcept;


    /**
     * @brief      Set boundary size.
     *
     * @param      size  The size
     */
    void setSize(units::Length size) noexcept;


    /**
     * @brief      Set boundary type.
     *
     * @param      type  The type
     */
    void setInletProfileType(InletProfileType type) noexcept;


    /**
     * @brief      Set inlet velocity.
     *
     * @param      velocity  The velocity
     */
    void setInletVelocity(units::Velocity velocity) noexcept;


    /**
     * @brief      Set inlet volumeric flow rate.
     *
     * @param      flow  The flow
     */
    void setInletFlow(units::VolumericFlow flow) noexcept;


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


// Private Data Members
private:

    /// Boundary name.
    String m_name;

    /// Boundary type.
    Type m_type = Type::None;

    /// Boundary position.
    Position m_position;

    /// Boundary position offset.
    units::Length m_offset = Zero;

    /// Boundary size.
    units::Length m_size = Zero;

    /// Inlet velocity profile type.
    InletProfileType m_inletProfileType = InletProfileType::Auto;

    /// Inlet velocity.
    units::Velocity m_inletVelocity = Zero;

    /// Inlet volumeric flow rate.
    units::VolumericFlow m_inletFlow = Zero;

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

inline const String& Boundary::getName() const noexcept
{
    return m_name;
}

/* ************************************************************************ */

inline Boundary::Type Boundary::getType() const noexcept
{
    return m_type;
}

/* ************************************************************************ */

inline Boundary::Position Boundary::getPosition() const noexcept
{
    return m_position;
}

/* ************************************************************************ */

inline units::Length Boundary::getOffset() const noexcept
{
    return m_offset;
}

/* ************************************************************************ */

inline units::Length Boundary::getSize() const noexcept
{
    return m_size;
}

/* ************************************************************************ */

inline Boundary::InletProfileType Boundary::getInletProfileType() const noexcept
{
    return m_inletProfileType;
}

/* ************************************************************************ */

inline units::Velocity Boundary::getInletVelocity() const noexcept
{
    return m_inletVelocity;
}

/* ************************************************************************ */

inline units::VolumericFlow Boundary::getInletFlow() const noexcept
{
    return m_inletFlow;
}

/* ************************************************************************ */

inline void Boundary::setName(String name) noexcept
{
    m_name = std::move(name);
}

/* ************************************************************************ */

inline void Boundary::setType(Type type) noexcept
{
    m_type = type;
}

/* ************************************************************************ */

inline void Boundary::setPosition(Position position) noexcept
{
    m_position = position;
}

/* ************************************************************************ */

inline void Boundary::setOffset(units::Length offset) noexcept
{
    m_offset = offset;
}

/* ************************************************************************ */

inline void Boundary::setSize(units::Length size) noexcept
{
    m_size = size;
}

/* ************************************************************************ */

inline void Boundary::setInletProfileType(InletProfileType type) noexcept
{
    m_inletProfileType = type;
}

/* ************************************************************************ */

inline void Boundary::setInletVelocity(units::Velocity velocity) noexcept
{
    m_inletVelocity = velocity;
}

/* ************************************************************************ */

inline void Boundary::setInletFlow(units::VolumericFlow flow) noexcept
{
    m_inletFlow = flow;
}

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
