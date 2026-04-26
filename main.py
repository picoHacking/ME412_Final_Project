# ME 412 Spring 2026 Final Project
# Amanda Yaklin and Christina Li
# 4/20/2026 - 5/10/2026

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
# | ---<-gw->---
# h ---      ---
# | ------------
# v ------------
#   <----W----->

ch_h = 1  # m total U channel height
ch_w = 2  # m individual channel width
ch_gw = 3  # m width of gap between vertical channels

# --- Grid
nx = 300
ny = 700

dx = (2 * ch_w + ch_gw)/nx # m per grid cell
dy = dx
# issue: something is wrong with our grid calculations
# --- Connect channels to grid
cx = int(1 / dx)  # grid cells per mm x
cy = int(1 / dy)  # grid cells per mm y
h = int(ch_h*cy)  # m total U channel height
w = int(ch_w*cx)  # m individual channel width
gw = int(ch_gw*cy)  # m gap width
total_w = int(2*w + gw)  # total U channel width

# --- Flow params
rho = 1000 # kg/m3
mu = 0.0010005 #Pa*s
U_in = 0  # m/s, keep at zero for U-shaped channel
V_in = -50*mu/(rho*ch_w)  # m/s (toward negative y)
Re = -V_in*rho*ch_w/mu # Required V_in: -50*mu/(rho*ch_w), -100*mu/(rho*ch_w), -200*mu/(rho*ch_w)

# --- Initialize grids
U = np.ones([nx, ny + 1])  # Note that U grid is offset 1/2 step to the right of P grid
V = np.zeros([nx + 1, ny])  # Note that V grid is offset 1/2 step up from P grid
P = np.zeros([nx + 1, ny + 1])

# --- Boundary Conditions
def apply_bcs(U, V, P):
    # --- For U --- #
    # No slip walls (horizontal boundaries have ghost cells)
    U[0][0:h]=0  # outer left wall
    U[2*w+gw][0:h]=0  # outer right wall
    U[0:2*w+gw][0]=-U[0:2*w+gw][1]  # outer bottom wall
    U[w:w+gw][w+1]=-U[w:w+gw][w]  # inner bottom wall
    U[w][w:h]=0 # inner left wall
    U[w+gw][w:h]=0 # inner right wall
    # Inlet
    U[0:w][h]=0
    # Outlet
    U[w+gw:2*w+gw][h]=0
    # --- For V --- #


apply_bcs(U, V, P)
print(np.min(V))
X = np.linspace(0, nx, nx)
Y = np.linspace(0, ny+1, ny+1)
X, Y = np.meshgrid(X, Y)
plt.contourf(X, Y, U.T)
plt.show()


