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

// Must be first
#include "../../python/Python.hpp"

// CeCe
#include "cece/core/StringStream.hpp"

// Diffusion
#include "../../streamlines/Module.hpp"

// Plugin
#include "../../python/Utils.hpp"
#include "../../python/Type.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines_python {

/* ************************************************************************ */

namespace {

/* ************************************************************************ */

using namespace plugin::python;

/* ************************************************************************ */

/**
 * @brief Module type.
 */
class ModuleType : public Type<plugin::streamlines::Module*>
{


// Ctors & Dtors
public:


    /**
     * @brief Constructor.
     */
    ModuleType()
        : Type("streamlines.Module")
    {
        tp_getset = m_properties;
        tp_methods = m_methods;
    }


// Public Operations
public:


    /**
     * @brief Returns lattice size.
     *
     * @param self
     *
     * @return
     */
    static PyObject* getLatticeSize(SelfType* self) noexcept
    {
        return makeObject(self->value->getLatticeSize()).release();
    }


    /**
     * @brief Check if there is a obstacle at given position.
     *
     * @param self
     * @param args
     *
     * @return
     */
    static PyObject* isObstacle(SelfType* self, PyObject* args) noexcept
    {
        CECE_ASSERT(self);
        CECE_ASSERT(self->value);
        CECE_ASSERT(args);

        int x;
        int y;

        if (!PyArg_ParseTuple(args, "ii", &x, &y))
            return nullptr;

        // Check if is in range
        const auto* module = self->value;
        const auto size = module->getLatticeSize();
        const auto coord = plugin::streamlines::Lattice::CoordinateType(x, y);

        if (!module->inLatticeRange(coord))
        {
            OutStringStream oss;
            oss << "Coordinates [" << x << ", " << y << "] out of range [" << size.getWidth() << ", " << size.getHeight() << "]";
            PyErr_SetString(PyExc_RuntimeError, oss.str().c_str());
            return nullptr;
        }

        if (module->getDynamics(coord) == plugin::streamlines::Dynamics::Wall)
            Py_RETURN_TRUE;
        else
            Py_RETURN_FALSE;
    }


    /**
     * @brief Change obstacle type at given position.
     *
     * @param self
     * @param args
     *
     * @return
     */
    static PyObject* setObstacle(SelfType* self, PyObject* args) noexcept
    {
        CECE_ASSERT(self);
        CECE_ASSERT(self->value);
        CECE_ASSERT(args);

        int x;
        int y;
        int flag;

        if (!PyArg_ParseTuple(args, "iii", &x, &y, &flag))
            return nullptr;

        // Check if is in range
        auto* module = self->value;
        const auto size = module->getLatticeSize();
        const auto coord = plugin::streamlines::Lattice::CoordinateType(x, y);

        if (!module->inLatticeRange(coord))
        {
            OutStringStream oss;
            oss << "Coordinates [" << x << ", " << y << "] out of range [" << size.getWidth() << ", " << size.getHeight() << "]";
            PyErr_SetString(PyExc_RuntimeError, oss.str().c_str());
            return nullptr;
        }

        // TODO: implement

        return none().release();
    }


    /**
     * @brief Returns velocity.
     *
     * @param self
     * @param args
     *
     * @return
     */
    static PyObject* getVelocity(SelfType* self, PyObject* args) noexcept
    {
        CECE_ASSERT(self);
        CECE_ASSERT(self->value);
        CECE_ASSERT(args);

        int x;
        int y;

        if (!PyArg_ParseTuple(args, "ii", &x, &y))
            return nullptr;

        // Check if is in range
        const auto* module = self->value;
        const auto size = module->getLatticeSize();
        const auto coord = plugin::streamlines::Lattice::CoordinateType(x, y);

        if (!module->inLatticeRange(coord))
        {
            OutStringStream oss;
            oss << "Coordinates [" << x << ", " << y << "] out of range [" << size.getWidth() << ", " << size.getHeight() << "]";
            PyErr_SetString(PyExc_RuntimeError, oss.str().c_str());
            return nullptr;
        }

        // Obtain velocity
        const auto vel = module->getVelocity(coord);

        return makeObject(vel).release();
    }


    /**
     * @brief Returns velocity.
     *
     * @param self
     * @param args
     *
     * @return
     */
    static PyObject* setVelocity(SelfType* self, PyObject* args) noexcept
    {
        CECE_ASSERT(self);
        CECE_ASSERT(self->value);
        CECE_ASSERT(args);

        int x;
        int y;
        PyObject* vel;

        if (!PyArg_ParseTuple(args, "iiO", &x, &y, &vel))
            return nullptr;

        // Check if is in range
        auto* module = self->value;
        const auto size = module->getLatticeSize();
        const auto coord = plugin::streamlines::Lattice::CoordinateType(x, y);

        if (!module->inLatticeRange(coord))
        {
            OutStringStream oss;
            oss << "Coordinates [" << x << ", " << y << "] out of range [" << size.getWidth() << ", " << size.getHeight() << "]";
            PyErr_SetString(PyExc_RuntimeError, oss.str().c_str());
            return nullptr;
        }

        // Set velocity
        module->setVelocity(coord, cast<units::VelocityVector>(vel));

        return none().release();
    }


// Private Data Members
private:

    /// Type properties.
    PyGetSetDef m_properties[2] = {
        {const_cast<char*>("latticeSize"), (getter) getLatticeSize, nullptr, nullptr},
        {nullptr}  /* Sentinel */
    };

    /// Type methods
    PyMethodDef m_methods[5] = {
        {"isObstacle", (PyCFunction) isObstacle, METH_VARARGS, nullptr},
        {"setObstacle", (PyCFunction) setObstacle, METH_VARARGS, nullptr},
        {"getVelocity", (PyCFunction) getVelocity, METH_VARARGS, nullptr},
        {"setVelocity", (PyCFunction) setVelocity, METH_VARARGS, nullptr},
        {nullptr}  /* Sentinel */
    };

};

/* ************************************************************************ */

ModuleType g_type;

/* ************************************************************************ */

}

/* ************************************************************************ */

void init_Module(PyObject* module)
{
    g_type.tp_base = g_type.getBaseType("simulator.Module");
    g_type.add(module);
}

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
