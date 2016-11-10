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
#include "cece/core/DynamicArray.hpp"

// Plugin
#include "Function.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace stochastic_reactions {

/* ************************************************************************ */

/**
 * @brief      User defined function.
 */
class UserFunction
{

// Public Ctors & Dtors
public:


    /**
     * @brief      Default ctor.
     */
    UserFunction() = default;


    /**
     * @brief      Constructor.
     *
     * @param[in]  name    The function name.
     * @param[in]  params  The function parameter names.
     * @param[in]  body    The function body.
     */
    explicit UserFunction(String name, DynamicArray<String> params, SharedPtr<Node<RealType>> body) noexcept
        : m_name(std::move(name))
        , m_parameters(std::move(params))
        , m_body(std::move(body))
    {
        // Nothing to do
    }


// Public Accessors & Mutators
public:


    /**
     * @brief      Returns function name.
     *
     * @return     The function name.
     */
    const String& getName() const noexcept
    {
        return m_name;
    }


    /**
     * @brief      Returns function parameters.
     *
     * @return     The parameters.
     */
    const DynamicArray<String>& getParameters() const noexcept
    {
        return m_parameters;
    }


// Public Operations
public:


    /**
     * @brief      Call user function.
     *
     * @param[in]  args  Call arguments.
     *
     * @return     Result value.
     */
    RealType call(const Context& context, const DynamicArray<RealType>& args = {}) const;


// Public Operations
public:

    /// Function name.
    String m_name;

    /// Parameter names.
    DynamicArray<String> m_parameters;

    /// Function body.
    SharedPtr<Node<RealType>> m_body;
};

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
