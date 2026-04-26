# ME 412 Spring 2026 Final Project
# Amanda Yaklin and Christina Li
# 4/20/2026 - 5/10/2026

# Note: numpy arrays are indexed as [y,x]

# --- Objectives
# Model 2D incompressible flow through a block U channel on a cartesian grid.
# Analyze the effect of Reynolds number and grid resolution on the flow field

# --- Import
import numpy as np
import matplotlib.pyplot as plt

# --- Channel Geometry

#  inlet    outlet
#   chw
#   <->
# ^ ---      ---
# | ---      ---
# h ---      ---
# | ------------
# v ------------
#   <----W----->

l_y = 1.0  # m total U channel height
l_w = 0.2  # m individual channel width
l_x = 1.0  # m total width

# --- Grid
nx = 64
ny = 108

dx = l_x / nx  # m per grid cell
dy = l_y / ny
# issue: something is wrong with our grid calculations
# --- Connect channels to grid
h = ny  #  total U channel height
w = nx  #  total U channel width
cw = int(l_w / dx)  #  individual channel width 

# --- Flow params
rho = 1000  # kg/m3
mu = 0.0010005  # Pa*s
Re = 100
U_in = 0  # m/s, keep at zero for U-shaped channel
V_in = -Re * mu / (rho * l_w)  # m/s (toward negative y)

# --- Initialize grids
# RESET TO ZEROS AFTER TESTING Note that U grid is offset 1/2 step to the right of P grid
U = np.ones([ny + 2, nx + 1])
# Note that V grid is offset 1/2 step up from P grid
V = np.ones([ny+1, nx + 2])
P = np.zeros([ny + 2, nx + 2])

U_star = np.zeros_like(U)
V_star = np.zeros_like(V)
rhs = np.zeros_like(P)

# --- Simulation parameters
max_iter = 5000
poisson_tol = 1e-5
poisson_max_iter = 1e+5
max_error = 0

# --- Boundary Conditions
def apply_bcs(U, V):
    # --- For U --- #
    # No slip walls (horizontal boundaries have ghost cells)
    U[0:h,0] = 0  # outer left wall
    U[0:h,w] = 0  # outer right wall
    U[0,0:w] = -U[1,0:w]  # outer bottom wall
    U[cw + 1,cw:w-cw] = -U[cw,cw:w-cw]  # inner bottom wall
    U[cw:h,cw] = 0  # inner left wall
    U[cw:h,w-cw] = 0  # inner right wall

    # --- For V --- #
    # No slip walls (vertical boundaries have ghost cells)
    V[0:h,1] = -V[0:h,0]  # outer left wall
    V[0:h,w] = -V[0:h,w-1]  # outer right wall
    V[0,0:w] = 0  # outer bottom wall
    V[cw,cw:w-cw] = 0  # inner bottom wall
    V[cw:h,cw - 1] = -V[cw:h,cw]  # inner left wall
    V[cw:h,w-cw + 1] = -V[cw:h,w-cw]  # inner right wall

    # --- Inlet/Outlet --- #
    # Inlet
    U[h,0:cw] = 0
    V[h,0:cw] = V_in
    # Outlet
    U[h,w-cw:w] = U[h - 1,w-cw:w]
    V[h,w-cw:w] = V[h - 1,w-cw:w]

def predictor():
    U_up = np.maximum(U,0) # positive
    U_dn = np.minimum(U,0) # negative
    du2dx = (U[1:ny+1,1:nx+1]+U[1:ny+1,2:nx+2])*\
        (U[1:ny+1,1:nx+1]+U[1:ny+1,2:nx+2]) - #wip
    
    
apply_bcs(U, V)
print(np.min(V))
X = np.linspace(0, nx+2, nx+2)
Y = np.linspace(0, ny + 1, ny + 1)
X, Y = np.meshgrid(X, Y)
plt.contourf(X, Y, V)
plt.show()
