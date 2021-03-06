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
#include "DrawableCell.hpp"

// C++
#include <cmath>
#include <array>

// CeCe
#include "cece/render/Context.hpp"
#include "cece/render/VertexFormat.hpp"

// Shaders
#include "vs.cell.hpp"
#include "fs.cell.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace cell {

/* ************************************************************************ */

struct Vertex
{
    float x, y;
    float u, v;
};

/* ************************************************************************ */

static const std::array<Vertex, 4> g_vertices = {{
    { 0.5f,  0.5f, 1.0f, 1.0f},
    { 0.5f, -0.5f, 1.0f, 0.0f},
    {-0.5f, -0.5f, 0.0f, 0.0f},
    {-0.5f,  0.5f, 0.0f, 1.0f}
}};

/* ************************************************************************ */

DrawableCell::DrawableCell(render::Context& context)
    : m_buffer(context, g_vertices.size() * sizeof(decltype(g_vertices)::value_type), g_vertices.data())
{
    m_vertexShader.init(render::Shader::Type::VERTEX, g_vertexShader);
    m_fragmentShader.init(render::Shader::Type::FRAGMENT, g_fragmentShader);
    m_program.init(m_vertexShader, m_fragmentShader);

    // Get size variable
    m_uniformSize  = m_program.getUniformId("g_size");
    m_uniformColor = m_program.getUniformId("g_color");
}

/* ************************************************************************ */

void DrawableCell::draw(render::Context& context, RealType scale, const render::Color& color) noexcept
{
    static render::VertexFormat vformat{
        render::VertexElement(render::VertexElementType::Position, render::DataType::Float, 2),
        render::VertexElement(render::VertexElementType::TexCoord, render::DataType::Float, 2)
    };

    // Set pointers
    context.setProgram(&m_program);
    context.setVertexBuffer(&m_buffer);
    context.setVertexFormat(&vformat);

    // Set parameters
    context.setProgramParam(m_uniformColor, color);

    // Draw circle
    context.draw(render::PrimitiveType::Quads, 4);

    context.setVertexFormat(nullptr);
    context.setVertexBuffer(nullptr);
    context.setProgram(nullptr);
}

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
