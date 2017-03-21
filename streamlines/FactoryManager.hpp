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
#include "cece/core/String.hpp"
#include "cece/core/StringView.hpp"
#include "cece/core/UniquePtr.hpp"
#include "cece/core/FactoryManager.hpp"

// Plugin
#include "Factory.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines {

/* ************************************************************************ */

/**
 * @brief      Lattice factory manager.
 */
class FactoryManager : public core::FactoryManager<Factory>
{

// Public Accessors & Mutators
public:


    /**
     * @brief      Register a new factory for specified lattice implementation.
     *
     * @param      name         Factory name.
     *
     * @tparam     LatticeType  Lattice type.
     */
    template<typename LatticeType>
    void createForLattice(String name) noexcept
    {
        createFor<LatticeType>(std::move(name));
    }


    /**
     * @brief      Gets the singleton instance.
     *
     * @return     The instance.
     */
    static FactoryManager& getInstance() noexcept;


// Public Operations
public:


    /**
     * @brief      Create a lattice by name.
     *
     * @param      name   Factory name.
     * @param[in]  size   The lattice size.
     * @param[in]  omega  The omega.
     *
     * @return     Created lattice.
     *
     * @throws     InitializerFactoryNotFoundException  In case of factory with given
     *                                                  name doesn't exists.
     */
    UniquePtr<Lattice> createLattice(StringView name, Lattice::SizeType size, RealType omega) const;

};

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
