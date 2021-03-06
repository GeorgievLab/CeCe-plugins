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

// CeCe
#include "cece/core/Real.hpp"
#include "cece/core/String.hpp"
#include "cece/core/Units.hpp"
#include "cece/core/VectorUnits.hpp"
#include "cece/core/DynamicArray.hpp"
#include "cece/core/IterationRange.hpp"
#include "cece/module/Module.hpp"
#include "cece/config/Configuration.hpp"

/* ************************************************************************ */

namespace cece { namespace plugin { namespace streamlines { class Module; } } }

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace object_streamlines_generator {

/* ************************************************************************ */

/**
 * @brief Structure for storing created object parameters.
 */
struct ObjectDesc
{
    /// Object class name.
    String className;

    /// Boundary name.
    String boundary;

    /// Object concentration.
    units::NumberConcentration concentration;

    /// List of iteration ranges when the generator is active.
    DynamicArray<IterationRange> active;

    /// Source offset
    units::Length offset = units::Length(1);

    /// Object configuration
    config::Configuration config;
};

/* ************************************************************************ */

/**
 * @brief Object generator module.
 */
class Module : public module::Module
{

// Public Ctors & Dtors
public:


    using module::Module::Module;


// Public Mutators
public:


    /**
     * @brief Register object.
     *
     * @param desc Object description.
     */
    void add(ObjectDesc desc)
    {
        m_objects.push_back(std::move(desc));
    }


// Public Operations
public:


    /**
     * @brief Load module configuration.
     *
     * @param config Source configuration.
     */
    void loadConfig(const config::Configuration& config) override;


    /**
     * @brief Store module configuration.
     *
     * @param config Output configuration.
     */
    void storeConfig(config::Configuration& config) const override;


    /**
     * @brief Initialize module.
     */
    void init() override;


    /**
     * @brief Update module state.
     */
    void update() override;


// Private Data Members
private:

    /// List of generated objects.
    DynamicArray<ObjectDesc> m_objects;

    /// Streamlines module.
    ViewPtr<streamlines::Module> m_module;

};

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
