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
#include "cece/core/Assert.hpp"
#include "cece/core/StringStream.hpp"

// Diffusion
#include "../../diffusion/Module.hpp"

// Plugin
#include "../../python/Type.hpp"
#include "../../python/Utils.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace diffusion_python {

/* ************************************************************************ */

namespace {

/* ************************************************************************ */

using namespace plugin::python;

/* ************************************************************************ */

/**
 * @brief Module type.
 */
class ModuleType : public Type<plugin::diffusion::Module*>
{


// Ctors & Dtors
public:


    /**
     * @brief Constructor.
     */
    ModuleType()
        : Type("diffusion.Module")
    {
        tp_getset = m_properties;
        tp_methods = m_methods;
    }


// Public Operations
public:


    /**
     * @brief Returns a number of signals.
     *
     * @param self
     *
     * @return
     */
    static PyObject* getSignalCount(SelfType* self) noexcept
    {
        CECE_ASSERT(self);
        CECE_ASSERT(self->value);
        return makeObject(self->value->getSignalCount()).release();
    }


    /**
     * @brief Returns grid size.
     *
     * @param self
     *
     * @return
     */
    static PyObject* getGridSize(SelfType* self) noexcept
    {
        CECE_ASSERT(self);
        CECE_ASSERT(self->value);
        return makeObject(self->value->getGridSize()).release();
    }


    /**
     * @brief Returns signal identifier.
     *
     * @param self
     * @param args
     *
     * @return
     */
    static PyObject* getSignalId(SelfType* self, PyObject* args) noexcept
    {
        CECE_ASSERT(self);
        CECE_ASSERT(self->value);
        CECE_ASSERT(args);

        char* name;

        if (!PyArg_ParseTuple(args, "s", &name))
            return nullptr;

        CECE_ASSERT(name);
        return makeObject(self->value->getSignalId(name)).release();
    }


    /**
     * @brief Returns signal value.
     *
     * @param self
     * @param args
     *
     * @return
     */
    static PyObject* getSignal(SelfType* self, PyObject* args) noexcept
    {
        CECE_ASSERT(self);
        CECE_ASSERT(self->value);
        CECE_ASSERT(args);

        PyObject* id;
        int x;
        int y;

        if (!PyArg_ParseTuple(args, "Oii", &id, &x, &y))
            return nullptr;

        // Check if is in range
        const auto size = self->value->getGridSize();
        const auto coord = plugin::diffusion::Module::Coordinate(x, y);

        if (!self->value->inRange(coord))
        {
            OutStringStream oss;
            oss << "Coordinates [" << x << ", " << y << "] out of range [" << size.getWidth() << ", " << size.getHeight() << "]";
            PyErr_SetString(PyExc_RuntimeError, oss.str().c_str());
            return nullptr;
        }

        if (PyInt_Check(id) || PyLong_Check(id))
        {
            return makeObject(
                self->value->getSignal(
                    plugin::diffusion::Module::SignalId(PyLong_AsLong(id)),
                    coord
                )
            ).release();
        }
        else if (PyString_Check(id))
        {
            return makeObject(
                self->value->getSignal(
                    PyString_AsString(id),
                    coord
                )
            ).release();
        }

        PyErr_SetString(PyExc_RuntimeError, "Invalid ID argument type, only int and string are allowed");

        return nullptr;
    }


    /**
     * @brief Set signal value.
     *
     * @param self
     * @param args
     *
     * @return
     */
    static PyObject* setSignal(SelfType* self, PyObject* args) noexcept
    {
        CECE_ASSERT(self);
        CECE_ASSERT(self->value);
        CECE_ASSERT(args);

        PyObject* id;
        int x;
        int y;
        double value;

        if (!PyArg_ParseTuple(args, "Oiid", &id, &x, &y, &value))
            return nullptr;

        // Check if is in range
        const auto size = self->value->getGridSize();
        const auto coord = plugin::diffusion::Module::Coordinate(x, y);

        if (!self->value->inRange(coord))
        {
            OutStringStream oss;
            oss << "Coordinates [" << x << ", " << y << "] out of range [" << size.getWidth() << ", " << size.getHeight() << "]";
            PyErr_SetString(PyExc_RuntimeError, oss.str().c_str());
            return nullptr;
        }

        if (PyInt_Check(id) || PyLong_Check(id))
        {
            self->value->setSignal(
                plugin::diffusion::Module::SignalId(PyLong_AsLong(id)),
                coord,
                plugin::diffusion::Module::SignalConcentration(value)
            );
        }
        else if (PyString_Check(id))
        {
            self->value->setSignal(
                PyString_AsString(id),
                coord,
                plugin::diffusion::Module::SignalConcentration(value)
            );
        }

        return none().release();
    }


    /**
     * @brief Set signal value.
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
        const auto size = self->value->getGridSize();
        const auto coord = plugin::diffusion::Module::Coordinate(x, y);

        if (!self->value->inRange(coord))
        {
            OutStringStream oss;
            oss << "Coordinates [" << x << ", " << y << "] out of range [" << size.getWidth() << ", " << size.getHeight() << "]";
            PyErr_SetString(PyExc_RuntimeError, oss.str().c_str());
            return nullptr;
        }

        // FIXME: segfault
        //return makeObject((int) self->value->isObstacle(coord)).release();

        if (self->value->isObstacle(coord))
            Py_RETURN_TRUE;
        else
            Py_RETURN_FALSE;
    }


// Private Data Members
private:

    /// Type Properties
    PyGetSetDef m_properties[3] = {
        {const_cast<char*>("signalCount"),  (getter) getSignalCount,    nullptr,  nullptr},
        {const_cast<char*>("gridSize"),     (getter) getGridSize,       nullptr,  nullptr},
        {nullptr}  /* Sentinel */
    };

    /// Type methods
    PyMethodDef m_methods[5] = {
        {"getSignalId", (PyCFunction) getSignalId,  METH_VARARGS, nullptr},
        {"getSignal",   (PyCFunction) getSignal,    METH_VARARGS, nullptr},
        {"setSignal",   (PyCFunction) setSignal,    METH_VARARGS, nullptr},
        {"isObstacle",  (PyCFunction) isObstacle,   METH_VARARGS, nullptr},
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
