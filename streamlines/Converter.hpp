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
#include "cece/core/Real.hpp"
#include "cece/core/Vector.hpp"
#include "cece/core/Units.hpp"
#include "cece/core/UnitsCtors.hpp"
#include "cece/core/VectorUnits.hpp"

/* ************************************************************************ */

namespace cece { namespace config { class Configuration; } }

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines {

/* ************************************************************************ */

/**
 * @brief      LB and physical units converter.
 */
class Converter
{


// Public Accessors & Mutators
public:


    /**
     * @brief      Returns fluid kinematic viscosity.
     *
     * @return     The kinematic viscosity.
     */
    const units::KinematicViscosity& getKinematicViscosity() const noexcept;


    /**
     * @brief      Set fluid kinematic viscosity.
     *
     * @param      viscosity  The kinematic viscosity.
     */
    void setKinematicViscosity(units::KinematicViscosity viscosity) noexcept;


    /**
     * @brief      Returns characteristic length.
     *
     * @return     The characteristic length.
     */
    const units::Length& getCharLength() const noexcept;


    /**
     * @brief      Set characteristic length.
     *
     * @param      length  The characteristic length.
     */
    void setCharLength(units::Length length) noexcept;


    /**
     * @brief      Returns characteristic time.
     *
     * @return     The characteristic time.
     */
    const units::Time& getCharTime() const noexcept;


    /**
     * @brief      Set characteristic time.
     *
     * @param      time  The characteristic time.
     */
    void setCharTime(units::Time time) noexcept;


    /**
     * @brief      Returns characteristic density.
     *
     * @return     The characteristic density.
     */
    const units::Density& getCharDensity() const noexcept;


    /**
     * @brief      Set characteristic density.
     *
     * @param      density  The characteristic density.
     */
    void setCharDensity(units::Density density) noexcept;


    /**
     * @brief      Returns characteristic velocity.
     *
     * @return     The characteristic velocity.
     */
    units::Velocity getCharVelocity() const noexcept;


    /**
     * @brief      Returns number of nodes in LB along characteristic length.
     *
     * @return     The number nodes.
     */
    unsigned int getNumberNodes() const noexcept;


    /**
     * @brief      Set number of nodes in LB for units conversion.
     *
     * @param      nodes  The nodes.
     */
    void setNumberNodes(unsigned int nodes) noexcept;


    /**
     * @brief      Returns number of time steps in LB for units conversion.
     *
     * @return     The number steps.
     */
    unsigned int getNumberSteps() const noexcept;


    /**
     * @brief      Set number of time steps in LB for units conversion.
     *
     * @param      steps  The steps.
     */
    void setNumberSteps(unsigned int steps) noexcept;


    /**
     * @brief      Return length conversion coefficient.
     *
     * @return     The length conversion coefficient.
     */
    units::Length getLengthCoefficient() const noexcept;


    /**
     * @brief      Return time conversion coefficient.
     *
     * @return     The time conversion coefficient.
     */
    units::Time getTimeCoefficient() const noexcept;


    /**
     * @brief      Return mass conversion coefficient.
     *
     * @return     The mass conversion coefficient.
     */
    units::Mass getMassCoefficient() const noexcept;


    /**
     * @brief      Return velocity conversion coefficient.
     *
     * @return     The velocity conversion coefficient.
     */
    units::Velocity getVelocityCoefficient() const noexcept;


    /**
     * @brief      Return force conversion coefficient.
     *
     * @return     The force conversion coefficient.
     */
    units::Force getForceCoefficient() const noexcept;


    /**
     * @brief      Return viscosity conversion coefficient.
     *
     * @return     The viscosity conversion coefficient.
     */
    units::KinematicViscosity getViscosityCoefficient() const noexcept;


// Public Operations
public:


    /**
     * @brief      Convert length from LB to physical.
     *
     * @param      length  The length in LB units.
     *
     * @return     Length in physical units.
     */
    units::Length convertLength(RealType length) const noexcept;


    /**
     * @brief      Convert length from physical to LB.
     *
     * @param      length  The length in physical units.
     *
     * @return     Length in LB units.
     */
    RealType convertLength(units::Length length) const noexcept;


    /**
     * @brief      Convert velocity from LB to physical.
     *
     * @param      vel   The velocity in LB units.
     *
     * @return     Velocity in physical units.
     */
    units::Velocity convertVelocity(RealType vel) const noexcept;


    /**
     * @brief      Convert velocity vector from LB to physical.
     *
     * @param      vel   The velocity vector in LB units.
     *
     * @return     Velocity vector in physical units.
     */
    units::VelocityVector convertVelocity(Vector<RealType> vel) const noexcept;


    /**
     * @brief      Convert velocity from physical to LB.
     *
     * @param      vel   The velocity in physical units.
     *
     * @return     Velocity in physical units.
     */
    RealType convertVelocity(units::Velocity vel) const noexcept;


    /**
     * @brief      Convert velocity vector from physical to LB.
     *
     * @param      vel   The velocity vector in physical units.
     *
     * @return     Velocity vector in physical units.
     */
    Vector<RealType> convertVelocity(units::VelocityVector vel) const noexcept;


    /**
     * @brief      Convert force from LB to physical.
     *
     * @param      force  The force in LB units.
     *
     * @return     Force in physical units.
     */
    units::Force convertForce(RealType force) const noexcept;


    /**
     * @brief      Convert force vector from LB to physical.
     *
     * @param      force  The force vector in LB units.
     *
     * @return     Force vector in physical units.
     */
    units::ForceVector convertForce(Vector<RealType> force) const noexcept;


    /**
     * @brief      Convert force from physical to LB.
     *
     * @param      force  The force in physical units.
     *
     * @return     Force in LB units.
     */
    RealType convertForce(units::Force force) const noexcept;


    /**
     * @brief      Convert force vector from physical to LB.
     *
     * @param      force  The force vector in physical units.
     *
     * @return     Force vector in LB units.
     */
    Vector<RealType> convertForce(units::ForceVector force) const noexcept;


    /**
     * @brief      Convert kinematic viscosity from physical to LB.
     *
     * @param      force  The kinematic viscosity in physical units.
     *
     * @return     Kinematic viscosity in LB units.
     */
    RealType convertViscosity(units::KinematicViscosity viscosity) const noexcept;


    /**
     * @brief      Get LB viscosity.
     *
     * @return     The viscosity.
     */
    RealType getViscosity() const noexcept;


    /**
     * @brief      Returns relaxation time.
     *
     * @return     The relaxation time.
     */
    RealType getTau() const noexcept;


    /**
     * @brief      Returns relaxation frequency (omega).
     *
     * @return     The omega.
     */
    RealType getOmega() const noexcept;


    /**
     * @brief      Returns Reynolds number.
     *
     * @return     The Reynolds number.
     */
    RealType getRe() const noexcept;


// Private Data Members
private:

    /// Fluid viscosity (water).
    units::KinematicViscosity m_kinematicViscosity = units::mm2_s(0.658);

    /// Characteristic length.
    units::Length m_charLength = units::um(1);

    /// Characteristic time.
    units::Time m_charTime = units::s(1);

    /// Characteristic density (water).
    units::Density m_charDensity = units::g(1) / units::m3(1e-6);

    /// Number of LB nodes for units conversions.
    unsigned int m_numberNodes = 1;

    /// Number of LB time steps for units conversions
    unsigned int m_numberSteps = 1;
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

inline const units::KinematicViscosity& Converter::getKinematicViscosity() const noexcept
{
    return m_kinematicViscosity;
}

/* ************************************************************************ */

inline void Converter::setKinematicViscosity(units::KinematicViscosity viscosity) noexcept
{
    m_kinematicViscosity = std::move(viscosity);
}

/* ************************************************************************ */

inline const units::Length& Converter::getCharLength() const noexcept
{
    return m_charLength;
}

/* ************************************************************************ */

inline void Converter::setCharLength(units::Length length) noexcept
{
    m_charLength = std::move(length);
}

/* ************************************************************************ */

inline const units::Time& Converter::getCharTime() const noexcept
{
    return m_charTime;
}

/* ************************************************************************ */

inline void Converter::setCharTime(units::Time time) noexcept
{
    m_charTime = std::move(time);
}

/* ************************************************************************ */

inline const units::Density& Converter::getCharDensity() const noexcept
{
    return m_charDensity;
}

/* ************************************************************************ */

inline void Converter::setCharDensity(units::Density density) noexcept
{
    m_charDensity = std::move(density);
}

/* ************************************************************************ */

inline units::Velocity Converter::getCharVelocity() const noexcept
{
    return getCharLength() / getCharTime();
}

/* ************************************************************************ */

inline unsigned int Converter::getNumberNodes() const noexcept
{
    return m_numberNodes;
}

/* ************************************************************************ */

inline void Converter::setNumberNodes(unsigned int nodes) noexcept
{
    CECE_ASSERT(nodes > 0);
    m_numberNodes = nodes;
}

/* ************************************************************************ */

inline unsigned int Converter::getNumberSteps() const noexcept
{
    return m_numberSteps;
}

/* ************************************************************************ */

inline void Converter::setNumberSteps(unsigned int steps) noexcept
{
    CECE_ASSERT(steps > 0);
    m_numberSteps = steps;
}

/* ************************************************************************ */

inline units::Length Converter::getLengthCoefficient() const noexcept
{
    return getCharLength() / getNumberNodes();
}

/* ************************************************************************ */

inline units::Time Converter::getTimeCoefficient() const noexcept
{
    return getCharTime() / getNumberSteps();
}

/* ************************************************************************ */

inline units::Mass Converter::getMassCoefficient() const noexcept
{
    const auto volume = getLengthCoefficient() * getLengthCoefficient() * getLengthCoefficient();
    return getCharDensity() * volume;
}

/* ************************************************************************ */

inline units::Velocity Converter::getVelocityCoefficient() const noexcept
{
    return getLengthCoefficient() / getTimeCoefficient();
}

/* ************************************************************************ */

inline units::Force Converter::getForceCoefficient() const noexcept
{
    return getMassCoefficient() * getLengthCoefficient() / (getTimeCoefficient() * getTimeCoefficient());
}

/* ************************************************************************ */

inline units::KinematicViscosity Converter::getViscosityCoefficient() const noexcept
{
    return (getLengthCoefficient() * getLengthCoefficient()) / getTimeCoefficient();
}

/* ************************************************************************ */

inline units::Length Converter::convertLength(RealType length) const noexcept
{
    return getLengthCoefficient() * length;
}

/* ************************************************************************ */

inline RealType Converter::convertLength(units::Length length) const noexcept
{
    return length / getLengthCoefficient();
}

/* ************************************************************************ */

inline units::Velocity Converter::convertVelocity(RealType vel) const noexcept
{
    return getVelocityCoefficient() * vel;
}

/* ************************************************************************ */

inline units::VelocityVector Converter::convertVelocity(Vector<RealType> vel) const noexcept
{
    return getVelocityCoefficient() * vel;
}

/* ************************************************************************ */

inline RealType Converter::convertVelocity(units::Velocity vel) const noexcept
{
    return vel / getVelocityCoefficient();
}

/* ************************************************************************ */

inline Vector<RealType> Converter::convertVelocity(units::VelocityVector vel) const noexcept
{
    return vel / getVelocityCoefficient();
}

/* ************************************************************************ */

inline units::Force Converter::convertForce(RealType force) const noexcept
{
    return getForceCoefficient() * force;
}

/* ************************************************************************ */

inline units::ForceVector Converter::convertForce(Vector<RealType> force) const noexcept
{
    return getForceCoefficient() * force;
}

/* ************************************************************************ */

inline RealType Converter::convertForce(units::Force force) const noexcept
{
    return force / getForceCoefficient();
}

/* ************************************************************************ */

inline Vector<RealType> Converter::convertForce(units::ForceVector force) const noexcept
{
    return force / getForceCoefficient();
}

/* ************************************************************************ */

inline RealType Converter::convertViscosity(units::KinematicViscosity viscosity) const noexcept
{
    return viscosity / getViscosityCoefficient();
}

/* ************************************************************************ */

inline RealType Converter::getViscosity() const noexcept
{
    return convertViscosity(getKinematicViscosity());
}

/* ************************************************************************ */

inline RealType Converter::getOmega() const noexcept
{
    return RealType(1.0) / getTau();
}

/* ************************************************************************ */

inline RealType Converter::getRe() const noexcept
{
    return getCharLength() * getCharLength() / getCharTime() / getKinematicViscosity();
}

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
