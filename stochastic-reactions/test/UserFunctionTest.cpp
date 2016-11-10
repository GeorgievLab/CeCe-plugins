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

// C++
#include <iostream>

// GTest
#include "gtest/gtest.h"

// Plugin
#include "../UserFunction.hpp"
#include "../Function.hpp"

/* ************************************************************************ */

using namespace cece;
using namespace cece::plugin::stochastic_reactions;

/* ************************************************************************ */

TEST(UserFunction, ctor)
{
    {
        UserFunction fn("function", {}, {});
        EXPECT_EQ("function", fn.getName());
        EXPECT_EQ(0, fn.getParameters().size());
    }

    {
        UserFunction fn("function", {"x", "y"}, {});
        EXPECT_EQ("function", fn.getName());
        EXPECT_EQ(2, fn.getParameters().size());
    }
}

/* ************************************************************************ */

TEST(UserFunction, call)
{
    {
        UserFunction fn("function", {}, makeUnique<Amount>(RealType(10)));

        EXPECT_FLOAT_EQ(10, fn.call(Context{}));
    }

    {
        UserFunction fn("function", {"x", "y"},
            // x + y
            makeUnique<OperatorTwo<std::plus<RealType>>>(
                makeUnique<IdentifierCell>("x"),
                makeUnique<IdentifierCell>("y")
            )
        );

        EXPECT_FLOAT_EQ(3, fn.call(Context{}, {RealType(1), RealType(2)}));
        EXPECT_FLOAT_EQ(25, fn.call(Context{}, {RealType(15), RealType(10)}));
    }
}

/* ************************************************************************ */
