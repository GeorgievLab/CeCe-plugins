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

// This must be first
#include "Python.hpp"

// CeCe
#include "cece/config.hpp"
#include "cece/core/String.hpp"
#include "cece/object/Object.hpp"
#ifdef CECE_RENDER
#  include "cece/render/Context.hpp"
#endif

// Plugin
#include "Handle.hpp"
#include "Source.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace python {

/* ************************************************************************ */

/**
 * @brief Python defined simulation object.
 */
class Object : public object::Object
{

// Public Ctors & Dtors
public:


    /**
     * @brief Constructor.
     *
     * @param name Module name. Can be path to python source.
     */
    explicit Object(simulator::Simulation& simulation, const String& name, Type type = Type::Static);


// Public Operations
public:


    /**
     * @brief Configure object.
     *
     * @param config
     * @param simulation
     */
    void configure(const config::Configuration& config, simulator::Simulation& simulation) override;


    /**
     * @brief Update object state.
     *
     * @param dt Simulation time step.
     */
    void update(units::Time dt) override;


#ifdef CECE_RENDER

    /**
     * @brief Render object.
     *
     * @param context
     */
    void draw(render::Context& context) override;

#endif


// Private Data Members
private:


    /// Source.
    Source m_source;

    /// Configure function.
    Handle<PyObject> m_configureFn;

    /// Update function.
    Handle<PyObject> m_updateFn;

#ifdef CECE_RENDER
    /// Draw function.
    Handle<PyObject> m_drawFn;
#endif

};

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
