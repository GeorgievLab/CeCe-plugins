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
#include "Program.hpp"

// CeCe
#include "cece/object/Object.hpp"
#include "cece/config/Configuration.hpp"
#include "cece/simulator/Simulation.hpp"

// Plugin
#include "Utils.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace python {

/* ************************************************************************ */

/// MSVC 2015 Update 3 contains a bug which doesn't include ref template function body into object file.
#if _MSC_FULL_VER >= 190024213
inline simulator::Simulation& ref(simulator::Simulation* ptr)
{
    return *ptr;
}

inline object::Object& ref(object::Object* ptr)
{
    return *ptr;
}
#endif

/* ************************************************************************ */

UniquePtr<program::Program> Program::clone() const
{
    return makeUnique<Program>(*this);
}

/* ************************************************************************ */

void Program::loadConfig(simulator::Simulation& simulation, const config::Configuration& config)
{
    if (config.has("filename"))
    {
        m_source.initFile(config.get("filename"));
    }
    else
    {
        m_source.initSource(config.getContent());
    }

    // Store function pointer
    m_call = m_source.getFunction("__call__");

    // Fallback
    if (!m_call)
        m_call = m_source.getFunction("call");

    if (!m_call)
        throw InvalidArgumentException("Python program doesn't have '__call__' or 'call' function");
}

/* ************************************************************************ */

void Program::call(simulator::Simulation& simulation, object::Object& object, units::Time dt)
{
    // Call function
    python::call(m_call, &object, &simulation, dt.value());
}

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
