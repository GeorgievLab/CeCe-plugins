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

const char PROGRAM[] = R"%(

/* ************************************************************************ */

#define DF_SIZE 9

/* ************************************************************************ */

__constant float SPEED_OF_SOUND_SQ = 1.0f / 3.0f;
__constant float SPEED_OF_SOUND_SQ_INV = 3.0f;
__constant float2 DIRS[DF_SIZE] = {{ 0,  0}, {-1,  1}, {-1,  0}, {-1, -1}, { 0, -1}, { 1, -1}, { 1,  0}, { 1,  1}, { 0,  1}};
__constant float WEIGHTS[DF_SIZE] = {4.0f / 9.0f, 1.0f / 36.0f, 1.0f / 9.0f, 1.0f / 36.0f, 1.0f / 9.0f, 1.0f / 36.0f, 1.0f / 9.0f, 1.0f / 36.0f, 1.0f / 9.0f};
__constant int OPPOSITE[DF_SIZE] = {0, 5, 6, 7, 8, 1, 2, 3, 4};

/* ************************************************************************ */

__constant int ZOUHE_CENTER_RHO[4][3] = {
    {8, 0, 4},
    {8, 0, 4},
    {2, 0, 6},
    {2, 0, 6}
};

__constant int ZOUHE_KNOWN_RHO[4][3] = {
    {7, 6, 5},
    {1, 2, 3},
    {1, 8, 7},
    {3, 4, 5}
};

__constant int ZOUHE_UNKNOWN_RHO[4][3] = {
    {1, 2, 3},
    {7, 6, 5},
    {3, 4, 5},
    {1, 8, 7}
};

__constant float2 ZOUHE_VELOCITIES[4] = {
    {-1,  0},
    { 1,  0},
    { 0, -1},
    { 0,  1}
};

__constant int ZOUHE_BC_CENTER[4] = { 2, 6, 4, 8 };

__constant int ZOUHE_BC_SIDE1[4][2] = {
    {1, 4},
    {7, 4},
    {5, 2},
    {7, 2}
};

__constant int ZOUHE_BC_SIDE2[4][2] = {
    {3, 8},
    {5, 8},
    {3, 6},
    {1, 6}
};

/* ************************************************************************ */

/// Dynamics
enum Dynamics
{
    Dynamics_None    = 0,    // No dynamics, data is not used
    Dynamics_Fluid   = 1,    // Fluid dynamics, BGK collision
    Dynamics_Wall    = 2,    // Wall dynamics, Bounce-Back
    Dynamics_Inlet   = 3,    // Inlet dynamics, ZouHe velocity BC
    Dynamics_Outlet  = 4,    // Outlet dynamics, ZouHe pressure BC
    Dynamics_Object  = 5     // Moving object dynamics
};

/* ************************************************************************ */

/// Boundary position
enum BoundaryPosition
{
    BoundaryPosition_Right  = 0,
    BoundaryPosition_Left   = 1,
    BoundaryPosition_Top    = 2,
    BoundaryPosition_Bottom = 3
};

/* ************************************************************************ */

/// Fluid dynamics data
struct FluidData { };

/* ************************************************************************ */

/// Wall dynamics data
struct WallData { };

/* ************************************************************************ */

/// Inlet dynamics data
struct InletData
{
    int position;
    float2 velocity;
};

/* ************************************************************************ */

/// Outlet dynamics data
struct OutletData
{
    int position;
    float density;
};

/* ************************************************************************ */

/// Object dynamics data
struct ObjectData
{
    float2 velocity;
};

/* ************************************************************************ */

/// Dynamics data
struct Data
{
    int dynamics;
    struct FluidData fluid;
    struct WallData wall;
    struct InletData inlet;
    struct OutletData outlet;
    struct ObjectData object;
};

/* ************************************************************************ */

int calc_offset(const int2 size, int x, int y)
{
    return y * size.x + x;
}

/* ************************************************************************ */

int calc_df_base(const int2 size, int iPop)
{
    return size.x * size.y * iPop;
}

/* ************************************************************************ */

int calc_df_offset(const int2 size, int x, int y, int iPop)
{
    return calc_df_base(size, iPop) + calc_offset(size, x, y);
}

/* ************************************************************************ */

// Calculate equilibrium
float calc_eq(const float2 u, const float rho, const int iPop)
{
    const float vu = dot(DIRS[iPop], u);
    const float uLen = length(u);
    const float uSqLen = uLen * uLen;

    return rho * WEIGHTS[iPop] * (
          1.0f
        + SPEED_OF_SOUND_SQ_INV * vu
        + 0.5f * SPEED_OF_SOUND_SQ_INV * SPEED_OF_SOUND_SQ_INV * vu * vu
        - 0.5f * SPEED_OF_SOUND_SQ_INV * uSqLen
    );
}

/* ************************************************************************ */

// Collide fluid node
void collideFluid(
    const int2 size,
    const int x,
    const int y,
    const float2 velocity,
    const float density,
    const float omega,
    __global float* df
)
{
    for (int iPop = 0; iPop < DF_SIZE; iPop++)
    {
        const int iF = calc_df_offset(size, x, y, iPop);

        // Calculate equilibrium distribution
        const float feq = calc_eq(velocity, density, iPop);

        // Collide
        df[iF] += -omega * (df[iF] - feq);
    }
}

/* ************************************************************************ */

// Collide wall node
void collideWall(
    const int2 size,
    const int x,
    const int y,
    __global float* df
)
{
    float tmp[9];

    for (int iPop = 1; iPop < 9; iPop++)
    {
        const int iF = calc_df_offset(size, x, y, iPop);

        tmp[iPop] = df[iF];
    }

    for (int iPop = 1; iPop < 9; iPop++)
    {
        const int iPopOpp = OPPOSITE[iPop];
        const int iFopp = calc_df_offset(size, x, y, iPopOpp);

        df[iFopp] = tmp[iPop];
    }
}

/* ************************************************************************ */

// Equilibrium difference
float calc_eq_diff(float2 velocity, float density, int iPop)
{
    return
        calc_eq(velocity, density, iPop) -
        calc_eq(velocity, density, OPPOSITE[iPop])
    ;
}

/* ************************************************************************ */

float calc_df_diff(const int2 size, __global const float* df, int iPop, int i)
{
    return df[calc_df_base(size, iPop) + i] - df[calc_df_base(size, OPPOSITE[iPop]) + i];
}

/* ************************************************************************ */

// Sum of 3 df values
float calc_sum_df_3(const int2 size, __global const float* df, int i, int iPop1, int iPop2, int iPop3)
{
    return
        df[calc_df_base(size, iPop1) + i] +
        df[calc_df_base(size, iPop2) + i] +
        df[calc_df_base(size, iPop3) + i]
    ;
}

/* ************************************************************************ */

// Compute missing df as defined by ZouHe BC.
void zou_he_init(
    const int2 size,
    int position,
    __global float* df,
    float2 velocity,
    float density,
    const int i
)
{
    const int center  = ZOUHE_BC_CENTER[position];
    const int side1_0 = ZOUHE_BC_SIDE1[position][0];
    const int side1_1 = ZOUHE_BC_SIDE1[position][1];
    const int side2_0 = ZOUHE_BC_SIDE2[position][0];
    const int side2_1 = ZOUHE_BC_SIDE2[position][1];

    const int center_opp  = OPPOSITE[ZOUHE_BC_CENTER[position]];
    const int side1_0_opp = OPPOSITE[ZOUHE_BC_SIDE1[position][0]];
    const int side1_1_opp = OPPOSITE[ZOUHE_BC_SIDE1[position][1]];
    const int side2_0_opp = OPPOSITE[ZOUHE_BC_SIDE2[position][0]];
    const int side2_1_opp = OPPOSITE[ZOUHE_BC_SIDE2[position][1]];

    const int center_i  = calc_df_base(size, center) + i;
    const int side1_0_i = calc_df_base(size, side1_0) + i;
    const int side1_1_i = calc_df_base(size, side1_1) + i;
    const int side2_0_i = calc_df_base(size, side2_0) + i;
    const int side2_1_i = calc_df_base(size, side2_1) + i;

    const int center_opp_i  = calc_df_base(size, center_opp) + i;
    const int side1_0_opp_i = calc_df_base(size, side1_0_opp) + i;
    const int side1_1_opp_i = calc_df_base(size, side1_1_opp) + i;
    const int side2_0_opp_i = calc_df_base(size, side2_0_opp) + i;
    const int side2_1_opp_i = calc_df_base(size, side2_1_opp) + i;

    // Center
    df[center_i] = df[center_opp_i]
        + calc_eq_diff(velocity, density, center)
    ;

    // Side 1
    df[side1_0_i] = df[side1_0_opp_i]
        + calc_eq_diff(velocity, density, side1_0)
        + 0.5 * calc_df_diff(size, df, side1_1, i)
    ;

    // Side 2
    df[side2_0_i] = df[side2_0_opp_i]
        + calc_eq_diff(velocity, density, side2_0)
        + 0.5 * calc_df_diff(size, df, side2_1, i)
    ;
}

/* ************************************************************************ */

void boundaryInlet(
    const int2 size,
    __global float* df,
    const struct Data data,
    const int i
)
{
    const int position = data.inlet.position;

    const float center = calc_sum_df_3(size, df, i,
        ZOUHE_CENTER_RHO[position][0],
        ZOUHE_CENTER_RHO[position][1],
        ZOUHE_CENTER_RHO[position][2]
    );

    const float known  = calc_sum_df_3(size, df, i,
        ZOUHE_KNOWN_RHO[position][0],
        ZOUHE_KNOWN_RHO[position][1],
        ZOUHE_KNOWN_RHO[position][2]
    );

    const float velP = dot(data.inlet.velocity, ZOUHE_VELOCITIES[position]);

    const float density = 1.0f / (1.0f - velP) * (center + 2.0f * known);

    zou_he_init(size, position, df, data.inlet.velocity, density, i);
}

/* ************************************************************************ */

void boundaryOutlet(
    const int2 size,
    __global float* df,
    const struct Data data,
    const int i
)
{
    const int position = data.inlet.position;

    const float center = calc_sum_df_3(size, df, i,
        ZOUHE_CENTER_RHO[position][0],
        ZOUHE_CENTER_RHO[position][1],
        ZOUHE_CENTER_RHO[position][2]
    );

    const float known  = calc_sum_df_3(size, df, i,
        ZOUHE_KNOWN_RHO[position][0],
        ZOUHE_KNOWN_RHO[position][1],
        ZOUHE_KNOWN_RHO[position][2]
    );

    // Speed
    const float speed = (1.0f - 1.0f / data.outlet.density * (center + 2.0f * known));

    // Velocity vector
    const float2 velocity = speed * ZOUHE_VELOCITIES[position];

    zou_he_init(size, position, df, velocity, data.outlet.density, i);
}

/* ************************************************************************ */

// Init
__kernel void init(const int2 size, const float2 u, const float rho, __global float* df)
{
    const int x = get_global_id(0);
    const int y = get_global_id(1);
    const int iPop = get_global_id(2);

    const int iF = calc_df_offset(size, x, y, iPop);

    df[iF] = calc_eq(u, rho, iPop);
}

/* ************************************************************************ */

// Collide df
__kernel void collide(
    const int2 size,
    __global const float2* u,
    __global const float* rho,
    __global float* df,
    __global const struct Data* data,
    const float omega
)
{
    const int x = get_global_id(0);
    const int y = get_global_id(1);

    const int i = calc_offset(size, x, y);

    switch (data[i].dynamics)
    {
    default:
        break;

    case Dynamics_Fluid:
    case Dynamics_Inlet:
    case Dynamics_Outlet:
        collideFluid(size, x, y, u[i], rho[i], omega, df);
        break;

    case Dynamics_Wall:
        collideWall(size, x, y, df);
        break;
    }
}

/* ************************************************************************ */

// Stream df
__kernel void stream(
    const int2 size,
    __global const float* fin,
    __global float* fout
)
{
    const int x = get_global_id(0);
    const int y = get_global_id(1);
    const int iPop = get_global_id(2);

    const int iF = calc_df_offset(size, x, y, iPop);
    const int iFmin = calc_df_offset(size, 0, 0, iPop);
    const int iFmax = calc_df_offset(size, size.x - 1, size.y - 1, iPop);

    // New coordinates
    const int xNew = x + DIRS[iPop].x;
    const int yNew = y + DIRS[iPop].y;

    if (!(xNew >= 0 && xNew < size.x))
        return;

    if (!(yNew >= 0 && yNew < size.y))
        return;

    const int iFnew = calc_df_offset(size, xNew, yNew, iPop);

    if (iFnew >= iFmin && iFnew <= iFmax)
        fout[iFnew] = fin[iF];
}

/* ************************************************************************ */

// Collide df
__kernel void bc(
    const int2 size,
    __global float* df,
    __global const struct Data* data
)
{
    const int x = get_global_id(0);
    const int y = get_global_id(1);

    const int i = calc_offset(size, x, y);

    switch (data[i].dynamics)
    {
    case Dynamics_Inlet:
        boundaryInlet(size, df, data[i], i);
        break;

    case Dynamics_Outlet:
        boundaryOutlet(size, df, data[i], i);
        break;

    default:
        break;
    }
}

/* ************************************************************************ */

// Synchronize df, velocity and density
__kernel void sync(
    const int2 size,
    __global const float* df,
    __global float2* u,
    __global float* rho,
    __global const struct Data* data
)
{
    const int x = get_global_id(0);
    const int y = get_global_id(1);

    const int i = calc_offset(size, x, y);

    if (data[i].dynamics == Dynamics_Wall)
        return;

    float ux = 0;
    float uy = 0;
    float r = 0;

    for (int iPop = 0; iPop < DF_SIZE; iPop++)
    {
        const int iF = calc_df_offset(size, x, y, iPop);
        const float fVal = df[iF];

        ux = ux + fVal * DIRS[iPop].x;
        uy = uy + fVal * DIRS[iPop].y;

        r = r + fVal;
    }

    rho[i] = r;
    u[i].x = ux;
    u[i].y = uy;
}

/* ************************************************************************ */

)%";

/* ************************************************************************ */
