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
#include "cece/core/UniquePtr.hpp"
#include "cece/core/UnitIo.hpp"

// Plugins
#include "../cell/CellBase.hpp"

// Plugin
#include "Diffusion.hpp"
#include "Types.hpp"
#include "Context.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace stochastic_reactions {

/* ************************************************************************ */

class UserFunction;

/* ************************************************************************ */

/**
 * @brief Parent of all function graph nodes.
 *
 */
template <typename Return>
struct Node
{
    virtual Return eval(const Context& context) const = 0;
};

/* ************************************************************************ */

/**
 * @brief Operator with 2 children.
 *
 */
template<typename OperatorType>
struct OperatorThree : public Node<typename OperatorType::result_type>
{
    using NodeOne = UniquePtr<Node<typename OperatorType::first_argument_type>>;
    using NodeTwo = UniquePtr<Node<typename OperatorType::second_argument_type>>;
    using NodeThree = UniquePtr<Node<typename OperatorType::third_argument_type>>;

private:
    NodeOne m_one;
    NodeTwo m_two;
    NodeThree m_three;

public:

    typename OperatorType::result_type eval(const Context& context) const override
    {
        return OperatorType{}(m_one->eval(context), m_two->eval(context), m_three->eval(context));
    }

public:

    OperatorThree(NodeOne one, NodeTwo two, NodeThree three):
    m_one(std::move(one)), m_two(std::move(two)), m_three(std::move(three))
    {
        // Nothing to do.
    }
};

/**
 * @brief Operator with 2 children.
 *
 */
template<typename OperatorType>
struct OperatorTwo : public Node<typename OperatorType::result_type>
{
    using NodeLeft = UniquePtr<Node<typename OperatorType::first_argument_type>>;
    using NodeRight = UniquePtr<Node<typename OperatorType::second_argument_type>>;

private:
    NodeLeft m_left;
    NodeRight m_right;

public:

    typename OperatorType::result_type eval(const Context& context) const override
    {
        return OperatorType{}(m_left->eval(context), m_right->eval(context));
    }

public:

    OperatorTwo(NodeLeft l, NodeRight r):
    m_left(std::move(l)), m_right(std::move(r))
    {
        // Nothing to do.
    }
};

/**
 * @brief Operator with 1 children.
 *
 */
template <typename OperatorType>
struct OperatorOne : public Node<typename OperatorType::result_type>
{
    using NodeRoot = UniquePtr<Node<typename OperatorType::argument_type>>;

private:
    NodeRoot root;

public:

    typename OperatorType::result_type eval(const Context& context) const override
    {
        return OperatorType{}(root->eval(context));
    }

public:

    OperatorOne(NodeRoot r):
    root(std::move(r))
    {
        // Nothing to do.
    }
};

/* ************************************************************************ */

/**
 * @brief Node for inner function.
 *
 */
struct Function : public Node<RealType>
{
private:
    SharedPtr<const UserFunction> function;
    DynamicArray<UniquePtr<Node<RealType>>> arguments;

public:

    RealType eval(const Context& context) const override;

public:

    Function(SharedPtr<const UserFunction> fn, DynamicArray<UniquePtr<Node<RealType>>> arguments)
        : function(fn)
        , arguments(std::move(arguments))
    {
        // Nothing to do.
    }
};

/* ************************************************************************ */

/**
 * @brief Leaf which uses context from cell.
 *
 */
struct IdentifierCell : public Node<RealType>
{
private:
    String m_identifier;

public:

    RealType eval(const Context& context) const override;

public:

    IdentifierCell(const String& identifier):
    m_identifier(identifier)
    {
        // Nothing to do.
    }
};

/**
 * @brief Leaf which uses context from diffusion.
 *
 */
struct IdentifierEnv : public Node<RealType>
{
private:
    String m_identifier;

public:

    RealType eval(const Context& context) const override
    {
        if (context.diffusion == nullptr)
            return 0;

        const auto id = context.diffusion->getSignalId(m_identifier);

        if (id == plugin::diffusion::Module::INVALID_SIGNAL_ID)
            return 0;

        return getMolarConcentration(*context.diffusion, *context.coords, id).value();
    }

public:

    IdentifierEnv(const String& identifier):
    m_identifier(identifier)
    {
        // Nothing to do.
    }
};

/**
 * @brief Leaf with fixed value.
 *
 */
struct Amount : public Node<RealType>
{
private:
    RealType m_amount;

public:

    RealType eval(const Context& context) const override
    {
        return m_amount;
    }

public:

    Amount(const RealType amount):
    m_amount(amount)
    {
        // Nothing to do.
    }
};

/**
 * @brief Leaf with value from global parameters.
 *
 */
struct IdentifierPar : public Node<RealType>
{
private:
    String m_identifier;

public:

    RealType eval(const Context& context) const override
    {
        return units::parse(context.parameters.get(m_identifier));
    }

public:

    IdentifierPar(const String& identifier):
    m_identifier(identifier)
    {
        // Nothing to do.
    }
};

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
