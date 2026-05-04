#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Nx 60
#define Ny 90
#define Lx 10.0
#define Ly 15.0
#define W 1.0 // width of channel
#define Re 100.0
#define rho 1000.0 // kg/m3
#define mu 0.0010005 // Pa*s
#define dt 0.001
#define max_iter 5000
#define poisson_tol 1e-5
#define poisson_max_iter 100000
double max_error;
int iter;
const double dx = Lx / Nx;
const double dy = Ly / Ny;
const int cw = W/dx; // width of channel

double u[Ny+2][Nx+1] = {0};       // u(i+1/2,j)
double v[Ny+1][Nx+2] = {0};       // v(i,j+1/2)
double p[Ny+2][Nx+2] = {0};       // p(i,j)
double u_star[Ny+2][Nx+1] = {0};
double v_star[Ny+1][Nx+2] = {0};
double rhs[Ny+2][Nx+2] = {0};

void apply_boundary_conditions() {
    int w = Nx;
    int h = Ny;
    double V_in = -Re * mu / (rho * W); // compute inlet velocity from Re

    // --- For u --- //

    // Outer left and right walls
    for (int j = 0; j < h; j++) {
        u[j][0] = 0.0;
        u[j][w] = 0.0;
    }

    // Outer bottom wall (ghost cell reflection)
    for (int i = 0; i < w; i++) {
        u[0][i] = -u[1][i];
    }

    // Inner bottom wall
    for (int i = cw; i < w - cw; i++) {
        u[cw][i] = -u[cw - 1][i];
    }

    // Inner left and right walls
    for (int j = cw; j < h; j++) {
        u[j][cw] = 0.0;
        u[j][w - cw] = 0.0;
    }


    // --- For v --- //

    // Outer left wall (ghost)
    for (int j = 0; j < h; j++) {
        v[j][1] = -v[j][0];
    }

    // Outer right wall (ghost)
    for (int j = 0; j < h; j++) {
        v[j][w] = -v[j][w - 1];
    }

    // Outer bottom wall
    for (int i = 0; i < w; i++) {
        v[0][i] = 0.0;
    }

    // Inner bottom wall
    for (int i = cw; i < w - cw; i++) {
        v[cw][i] = 0.0;
    }

    // Inner left wall (ghost)
    for (int j = cw; j < h; j++) {
        v[j][cw - 1] = -v[j][cw];
    }

    // Inner right wall (ghost)
    for (int j = cw; j < h; j++) {
        v[j][w - cw + 1] = -v[j][w - cw];
    }


    // --- Inlet --- //
    for (int i = 0; i < cw; i++) {
        u[h][i] = 0.0;
        v[h][i] = V_in;
    }

    // --- Outlet --- //
    for (int i = w - cw; i < w; i++) {
        u[h][i] = u[h - 1][i];
        v[h][i] = v[h - 1][i];
    }
}

double max(double x, double y) {
    return (x > y) ? x : y;
}

double min(double x, double y) {
    return (x < y) ? x : y;
}

void predictor() {
    for (int j = 1; j < Ny+1; j++) {
        for (int i = 1; i < Nx; i++) {
            double up = 0.5*(u[j][i+1] + u[j][i]);
            double vp = 0.25*(v[j][i] + v[j][i+1] + v[j+1][i] + v[j+1][i+1]);
            double du2dx = max(0,up) * (up - 0.5 * (u[j][i-1] + u[j][i])) / dx + min(0,up) * (0.5 * (u[j][i+2] + u[j][i+1]) - up) / dx;
            double duvdy = max(0,vp) * (u[j+1][i] - u[j-1][i]) / (2 * dy) + min(0,vp) * (u[j+1][i] - u[i][j]) / (2 * dy);

            double d2udx2 = (u[j][i+1] - 2*u[j][i] + u[j][i-1]) / (dx * dx);
            double d2udy2 = (u[j+1][i] - 2*u[j][i] + u[j-1][i]) / (dy * dy);

            u_star[j][i] = u[j][i] + dt * (-du2dx - duvdy + (d2udx2 + d2udy2)/Re);
        }
    }

    for (int j = 1; j < Ny; j++) {
        for (int i = 1; i < Nx+1; i++) {
            double vp = 0.5*(v[j+1][i] + v[j][i]);
            double up = 0.25*(u[j][i] + u[j][i+1] + u[j+1][i] + u[j+1][i+1]);
            double dv2dy = max(0,vp) * (vp - 0.5*(v[j-1][i] + v[j][i]))/dy + min(0,vp)*(0.5*(v[j+2][i] + v[j+1][i]) - vp)/dy;

            double duvdx = max(0,up) * (v[j][i+1] - u[j][i-1]) / (2 * dx) + min(0,up) * (u[j][i+1] - u[i][j]) / (2 * dx);
            double d2vdx2 = (v[j][i+1] - 2*v[j][i] + v[j][i-1]) / (dx * dx);
            double d2vdy2 = (v[j+1][i] - 2*v[j][i] + v[j-1][i]) / (dy * dy);

            v_star[j][i] = v[j][i] + dt * (-duvdx - dv2dy + (d2vdx2 + d2vdy2)/Re);
        }
    }
}

void find_divergence() {
    for (int j = 1; j < Ny+1; j++) {
        for (int i = 1; i < Nx+1; i++) {
            rhs[j][i] = ((u_star[j][i] - u_star[j][i-1]) / dx + (v_star[j][i] - v_star[j-1][i]) / dy) / dt;
        }
    }
}

void solve_poisson() {
    for (iter = 0; iter < poisson_max_iter; iter++) {
        max_error = 0.0;
        for (int j = 1; j < Ny+1; j++) {
            for (int i = 1; i < Nx+1; i++) {
                double p_old = p[j][i];
                p[j][i] = (
                    (p[j][i+1] + p[j][i-1]) * dy * dy +
                    (p[j+1][i] + p[j-1][i]) * dx * dx -
                    rhs[j][i] * dx * dx * dy * dy
                ) / (2.0 * (dx * dx + dy * dy));
                double err = fabs(p[j][i] - p_old);
                if (err > max_error) max_error = err;
            }
        }
        for (int i = 0; i < Nx+2; i++) {
            p[0][i]     = p[1][i];       // bottom
            p[Ny+1][i]  = p[Ny][i];      // top
        }
        for (int j = 0; j < Ny+2; j++) {
            p[j][0]     = p[j][1];       // left
            p[j][Nx+1]  = p[j][Nx];      // right
        }
        // Fix pressure reference
        p[Ny/2][Nx/2] = 0.0;

        if (max_error < poisson_tol) break;
    }
}

void corrector(){
    for (int j = 1; j < Ny+1; j++) {
        for (int i = 1; i < Nx; i++) {
            u[j][i] = u_star[j][i] - dt * (p[j][i+1] - p[j][i]) / dx;
        }
    }

    for (int j = 1; j < Ny; j++) {
        for (int i = 1; i < Nx+1; i++) {
            v[j][i] = v_star[j][i] - dt * (p[j+1][i] - p[j][i]) / dy;
        }
    }
}

void write_vtk(int step) {
    char filename[64];
    snprintf(filename, sizeof(filename), "output_%04d.vtk", step);
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        perror("File opening failed");
        return;
    }

    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "Lid-driven cavity flow\nASCII\n");
    fprintf(fp, "DATASET STRUCTURED_POINTS\n");
    fprintf(fp, "DIMENSIONS %d %d 1\n", Nx, Ny);
    fprintf(fp, "ORIGIN 0 0 0\n");
    fprintf(fp, "SPACING %f %f 1\n", dx, dy);
    fprintf(fp, "POINT_DATA %d\n", Nx * Ny);

    // Pressure field
    fprintf(fp, "SCALARS pressure float 1\nLOOKUP_TABLE default\n");
    for (int j = 1; j < Ny+1; j++) {
        for (int i = 1; i < Nx+1; i++) {
            fprintf(fp, "%f\n", p[j][i]);
        }
    }

    // Velocity magnitude
    fprintf(fp, "SCALARS velocity_magnitude float 1\nLOOKUP_TABLE default\n");
    for (int j = 1; j < Ny+1; j++) {
        for (int i = 1; i < Nx+1; i++) {
            double uc = 0.5 * (u[j][i] + u[j][i-1]);
            double vc = 0.5 * (v[j][i] + v[j-1][i]);
            fprintf(fp, "%f\n", sqrt(uc * uc + vc * vc));
        }
    }

    // Vorticity (z-component of curl)
    fprintf(fp, "SCALARS vorticity float 1\nLOOKUP_TABLE default\n");
    for (int j = 1; j < Ny+1; j++) {
        for (int i = 1; i < Nx+1; i++) {
            double dudy = (u[j+1][i] - u[j-1][i]) / (2.0 * dy);
            double dvdx = (v[j][i+1] - v[j][i-1]) / (2.0 * dx);
            double omega = dvdx - dudy;
            fprintf(fp, "%f\n", omega);
        }
    }

    // u component
    fprintf(fp, "SCALARS u float 1\nLOOKUP_TABLE default\n");
    for (int j = 1; j < Ny+1; j++) {
        for (int i = 1; i < Nx+1; i++) {
            double uc = 0.5 * (u[j][i] + u[j][i-1]);
            fprintf(fp, "%f\n", uc);
        }
    }

    // v component
    fprintf(fp, "SCALARS v float 1\nLOOKUP_TABLE default\n");
    for (int j = 1; j < Ny+1; j++) {
        for (int i = 1; i < Nx+1; i++) {
            double vc = 0.5 * (v[j][i] + v[j-1][i]);
            fprintf(fp, "%f\n", vc);
        }
    }

    // Velocity vector
    fprintf(fp, "VECTORS velocity float\n");
    for (int j = 1; j < Ny+1; j++) {
        for (int i = 1; i < Nx+1; i++) {
            double uc = 0.5 * (u[j][i] + u[j][i-1]);
            double vc = 0.5 * (v[j][i] + v[j-1][i]);
            fprintf(fp, "%f %f 0.0\n", uc, vc);
        }
    }

    fclose(fp);
}


int main() {
    for (int step = 0; step <= max_iter; step++) {
        apply_boundary_conditions();
        predictor();
        find_divergence();
        solve_poisson();
        corrector();

        printf("Step=%d iter=%d maxErr=%e u=%.4f v=%.4f p=%.4f\n",step, iter, max_error, u[Ny/2][Nx/2], v[Ny/2][Nx/2],p[Ny/2][Nx/2]);
        if (step % 500 == 0) write_vtk(step);
    }

    printf("Simulation complete.\n");
    return 0;
}
