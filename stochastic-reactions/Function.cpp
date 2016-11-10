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
    return function->call(context, args);
}

/* ************************************************************************ */

RealType IdentifierCell::eval(const Context& context) const
{
    if (context.arguments != nullptr)
    {
        // Function context
        auto it = context.arguments->find(m_identifier);
        if (it != context.arguments->end())
            return it->second;
    }

    if (context.cell == nullptr)
        return 0;

    // Cell context
    return context.cell->getMoleculeCount(m_identifier);
}

/* ************************************************************************ */

RealType IdentifierEnv::eval(const Context& context) const
{
    if (context.diffusion == nullptr)
        return 0;

    const auto id = context.diffusion->getSignalId(m_identifier);

    if (id == plugin::diffusion::Module::INVALID_SIGNAL_ID)
        return 0;

    return getMolarConcentration(*context.diffusion, *context.coords, id).value();
}

/* ************************************************************************ */

RealType IdentifierEnvNo::eval(const Context& context) const
{
    if (context.diffusion == nullptr || context.simulation == nullptr)
        return 0;

    const auto id = context.diffusion->getSignalId(m_identifier);

    if (id == plugin::diffusion::Module::INVALID_SIGNAL_ID)
        return 0;

    return getAmountOfMolecules(*context.simulation, *context.diffusion, *context.coords, id);
}

/* ************************************************************************ */

RealType IdentifierPar::eval(const Context& context) const
{
    if (context.simulation == nullptr)
        return 0;

    return units::parse(context.simulation->getParameters().get(m_identifier));
}

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
