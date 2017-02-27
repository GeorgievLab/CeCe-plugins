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
#include "Converter.hpp"

// CeCe
#include "cece/core/UnitIo.hpp"
#include "cece/config/Configuration.hpp"
#include "cece/config/Exception.hpp"

// Plugin
#include "Descriptor.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines {

/* ************************************************************************ */

units::Length Converter::convertLength(RealType length) const noexcept
{
    const auto charLength = m_charLength / getNumberNodes();
    return charLength * length;
}

/* ************************************************************************ */

RealType Converter::convertLength(units::Length length) const noexcept
{
    const auto charLength = m_charLength / getNumberNodes();
    return length / charLength;
}

/* ************************************************************************ */

units::Velocity Converter::convertVelocity(RealType vel) const noexcept
{
    const auto charTime = m_charTime / getNumberSteps();
    const auto charLength = m_charLength / getNumberNodes();
    return charLength / charTime * vel;
}

/* ************************************************************************ */

units::VelocityVector Converter::convertVelocity(Vector<RealType> vel) const noexcept
{
    const auto charTime = m_charTime / getNumberSteps();
    const auto charLength = m_charLength / getNumberNodes();
    return charLength / charTime * vel;
}

/* ************************************************************************ */

RealType Converter::convertVelocity(units::Velocity vel) const noexcept
{
    const auto charTime = m_charTime / getNumberSteps();
    const auto charLength = m_charLength / getNumberNodes();
    return charTime / charLength * vel;
}

/* ************************************************************************ */

Vector<RealType> Converter::convertVelocity(units::VelocityVector vel) const noexcept
{
    const auto charTime = m_charTime / getNumberSteps();
    const auto charLength = m_charLength / getNumberNodes();
    return charTime / charLength * vel;
}

/* ************************************************************************ */

units::Force Converter::convertForce(RealType force) const noexcept
{
    const auto charTime = m_charTime / getNumberSteps();
    const auto charLength = m_charLength / getNumberNodes();
    const auto charDensity = m_charDensity;
    const auto charMass = charDensity * (charLength * charLength * charLength);
    return charMass * charLength / (charTime * charTime) * force;
}

/* ************************************************************************ */

units::ForceVector Converter::convertForce(Vector<RealType> force) const noexcept
{
    const auto charTime = m_charTime / getNumberSteps();
    const auto charLength = m_charLength / getNumberNodes();
    const auto charDensity = m_charDensity;
    const auto charMass = charDensity * (charLength * charLength * charLength);
    return charMass * charLength / (charTime * charTime) * force;
}

/* ************************************************************************ */

RealType Converter::convertForce(units::Force force) const noexcept
{
    const auto charTime = m_charTime / getNumberSteps();
    const auto charLength = m_charLength / getNumberNodes();
    const auto charDensity = m_charDensity;
    const auto charMass = charDensity * (charLength * charLength * charLength);
    return (charTime * charTime) / (charMass * charLength) * force;
}

/* ************************************************************************ */

Vector<RealType> Converter::convertForce(units::ForceVector force) const noexcept
{
    const auto charTime = m_charTime / getNumberSteps();
    const auto charLength = m_charLength / getNumberNodes();
    const auto charDensity = m_charDensity;
    const auto charMass = charDensity * (charLength * charLength * charLength);
    return (charTime * charTime) / (charMass * charLength) * force;
}

/* ************************************************************************ */

RealType Converter::calculateViscosity() const noexcept
{
    const auto charTime = getCharTime() / getNumberSteps();
    const auto charLength = getCharLength() / getNumberNodes();

    return charTime / (charLength * charLength) * getKinematicViscosity();
}

/* ************************************************************************ */

RealType Converter::calculateTau() const noexcept
{
    return Descriptor::SPEED_OF_SOUND_SQ_INV * calculateViscosity() + 0.5;
}

/* ************************************************************************ */

unsigned int Converter::calculateNumberSteps(RealType tau) const noexcept
{
    return Descriptor::SPEED_OF_SOUND_SQ_INV *
        (getCharTime() * getNumberNodes() * getNumberNodes() * getKinematicViscosity())
        /
        ((tau - 0.5) * getCharLength() * getCharLength())
    ;
}

/* ************************************************************************ */

void Converter::loadConfig(const config::Configuration& config)
{
    // Viscosity
    setKinematicViscosity(config.get("kinematic-viscosity", getKinematicViscosity()));

    if (getKinematicViscosity() == Zero)
        throw config::Exception("Kinematic viscosity ('kinematic-viscosity') cannot be zero");

    // Characteristic length & time
    setCharLength(config.get("char-length", getCharLength()));
    setNumberNodes(config.get("number-nodes", getNumberNodes()));

    if (config.has("tau"))
        setNumberSteps(calculateNumberSteps(config.get<RealType>("tau")));
    else
        // Set number of time steps
        setNumberSteps(config.get("number-steps", getNumberNodes() * getNumberNodes() * 20));

    // Set channel height
    setHeight(config.get("height", units::Length(1)));
}

/* ************************************************************************ */

void Converter::storeConfig(config::Configuration& config) const
{
    config.set("kinematic-viscosity", getKinematicViscosity());
    config.set("char-length", getCharLength());
    config.set("number-nodes", getNumberNodes());
    config.set("number-steps", getNumberSteps());
}

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
