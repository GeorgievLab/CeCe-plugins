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
#include "Function.hpp"

// CeCe
#include "cece/core/Assert.hpp"

// Plugin
#include "UserFunction.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace stochastic_reactions {

/* ************************************************************************ */

RealType Function::eval(const Context& context) const
{
    DynamicArray<RealType> args;

    for (int i = 0; i < arguments.size(); ++i)
        args.push_back(arguments[i]->eval(context));

    CECE_ASSERT(function);
    return function->call(args);
}

/* ************************************************************************ */

RealType IdentifierCell::eval(const Context& context) const
{
    if (context.cell)
    {
        // Cell context
        return context.cell->getMoleculeCount(m_identifier);
    }

    // Function context
    auto it = context.arguments.find(m_identifier);
    if (it == context.arguments.end())
        throw InvalidArgumentException("Argument `" + m_identifier + "` not found");

    return it->second;
}

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
