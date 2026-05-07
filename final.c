#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Nx 60
#define Ny 90
#define Lx 0.5          //10
#define Ly 0.75         //15
#define W 0.1 // width of channel
#define Re 100.0
#define rho 1.0 // kg/m3
#define mu 0.0010005 // Pa*s
#define dt 0.0002
#define max_iter 5000
#define poisson_tol 1e-20 // trust in the power of FVM
#define poisson_max_iter 10000
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
int fluid[Ny][Nx] = {0}; // where the fluid is = where the barriers aren't

void define_fluid_region() {
    for (int j = 1; j < Ny; j++){
        for (int i = 1; i < cw; i++){
            fluid[j][i] = 1;
            fluid[j][Nx-cw+i] = 1;
        }
    }
    for (int j = 1; j < cw; j++){
        for (int i = 1; i < Nx; i++){
            fluid[j][i] = 1;
        }
    }

}

void apply_boundary_conditions() {
    double V_in = -Re * mu / (rho * W); // compute inlet velocity from Re
    // --- For u --- //
    // assignment of fluid regions is included here
    // Outer left and right walls
    for (int j = 0; j < Ny; j++) {
        u[j][0] = 0.0;
        u[j][Nx+1] = 0.0;
    
    }

    // Outer bottom wall (ghost cell reflection)
    for (int i = 0; i < Nx; i++) {
        u[1][i] = -u[0][i];
    }

    // Inner bottom wall
    for (int i = cw; i < Nx - cw; i++) {
        u[cw][i] = -u[cw+1][i];
    }

    // Inner left and right walls
    for (int j = cw; j < Ny; j++) {
        u[j][cw] = 0.0;
        u[j][Nx - cw+1] = 0.0;
    }


    // --- For v --- //

    // Outer left wall (ghost)
    for (int j = 0; j < Ny; j++) {
        v[j][1] = -v[j][0];
    }

    // Outer right wall (ghost)
    for (int j = 0; j < Ny; j++) {
        v[j][Nx] = -v[j][Nx+1];
    }

    // Outer bottom wall
    for (int i = 0; i < Nx; i++) {
        v[0][i] = 0.0;
    }

    // Inner bottom wall
    for (int i = cw; i < Nx - cw; i++) {
        v[cw][i] = 0.0;
    }

    // Inner left wall (ghost)
    for (int j = cw; j < Ny; j++) {
        v[j][cw] = -v[j][cw+1];
    }

    // Inner right wall (ghost)
    for (int j = cw; j < Ny; j++) {
        v[j][Nx-cw+1] = -v[j][Nx-cw];
    }


    // --- Inlet --- //
    for (int i = 1; i < cw+1; i++) {
        u[Ny][i] = 0.0;
        v[Ny][i] = V_in;
    }

    // --- Outlet --- //
    for (int i = Nx - cw; i < Nx; i++) {
        u[Ny][i] = u[Ny - 1][i];
        v[Ny][i] = v[Ny - 1][i];
    }
}

double max(double x, double y) {
    return (x > y) ? x : y;
}

double min(double x, double y) {
    return (x < y) ? x : y;
}

void apply_pressure_bc(void)
{
    int i, j;
    /* bottom and top: dp/dn = 0    Nx=w */
    /*for (i = 0; i < Nx + 2; i++)
    {
        p[0][i]    = p[1][i];
        p[Ny+1][i] = p[Ny][i]; 
    }

    left and right: dp/dn = 0 
    for (j = 0; j < Ny + 2; j++)
    {
        p[j][0]    = p[j][1];
        p[j][Nx+1] = p[j][Nx];
    }*/

    /* outlet: p = 0 */ 
    /* for U-channel outlet at top of right branch */
    for (i = Nx - cw; i <= Nx; i++)
    {
        p[Ny][i] = 0.0;
    }

    // top face
    for (i = 1; i < Nx-cw; i++)
    {
        p[Ny-1][i] = p[Ny][i];
    }

    // Outer left wall (ghost)
    for (int j = 0; j < Ny; j++) {
        p[j][0] = p[j][1];
    }

    // Outer right wall (ghost)
    for (int j = 0; j < Ny; j++) {
        p[j][Nx + 1] = p[j][Nx];
    }

    // Outer bottom wall (ghost)
    for (int i = 0; i < Nx; i++) {
        p[0][i] = p[1][i];
    }

    // Inner bottom wall (ghost)
    for (int i = cw; i < Nx - cw; i++) {
        p[cw][i] = p[cw+1][i];
    }

    // Inner left wall (ghost)
    for (int j = cw; j < Ny; j++) {
        p[j][cw+1] = p[j][cw];
    }

    // Inner right wall (ghost)
    for (int j = cw; j < Ny; j++) {
        p[j][Nx - cw] = p[j][Nx - cw+1];
    }
}

void predictor() {
    for (int j = 1; j < Ny+1; j++) {
        for (int i = 1; i < Nx; i++) {
            double ue_i = 0.5*(u[j][i]+u[j][i+1]);
            double uw_i = 0.5*(u[j][i]+u[j][i-1]);
            double vn_i = 0.25*(v[j][i]+v[j+1][i]+v[j][i+1]+v[j+1][i+1]);
            double vs_i = 0.25*(v[j][i]+v[j-1][i]+v[j][i+1]+v[j-1][i+1]);

            double du2dx = (max(0,ue_i)*u[j][i] + min(0,ue_i)*u[j][i+1] - max(0,uw_i)*u[j][i-1] - min(0,uw_i)*u[j][i])/dx;
            double duvdy = (max(0,vn_i)*u[j][i] + min(0,vn_i)*u[j+1][i] - max(0,vs_i)*u[j-1][i] - min(0,vs_i)*u[j][i])/dy;
            
            double d2udx2 = (u[j][i+1] - 2*u[j][i] + u[j][i-1]) / (dx * dx);
            double d2udy2 = (u[j+1][i] - 2*u[j][i] + u[j-1][i]) / (dy * dy);
            
            u_star[j][i] = u[j][i] + dt * (-du2dx - duvdy + (d2udx2 + d2udy2)/Re);
        }
    }

    for (int j = 1; j < Ny; j++) {
        for (int i = 1; i < Nx+1; i++) {
            double vn_j = 0.5*(v[j][i]+v[j+1][i]);
            double vs_j = 0.5*(v[j][i]+v[j-1][i]);
            double ue_j = 0.25*(u[j][i]+u[j][i+1]+u[j+1][i]+u[j+1][i+1]);
            double uw_j = 0.25*(u[j][i]+u[j][i-1]+u[j+1][i]+u[j+1][i-1]);

            double dv2dy = (max(0,vn_j)*v[j][i] + min(0,vn_j)*v[j+1][i] - max(0,vs_j)*v[j-1][i] - min(0,vs_j)*v[j][i])/dy;
            double duvdx = (max(0,ue_j)*u[j][i] + min(0,ue_j)*u[j][i+1] - max(0,uw_j)*u[j][i-1] - min(0,uw_j)*u[j][i])/dx;
            
            double d2vdx2 = (v[j][i+1] - 2*v[j][i] + v[j][i-1]) / (dx * dx);
            double d2vdy2 = (v[j+1][i] - 2*v[j][i] + v[j-1][i]) / (dy * dy);

            v_star[j][i] = v[j][i] + dt * (-dv2dy - duvdx + (d2vdx2 + d2vdy2)/Re);

        }
    }
}

void find_divergence() {
    double mean_rhs = 0.0;
    int count = 0;
    for (int j = 1; j < Ny+1; j++) {
        for (int i = 1; i < Nx+1; i++) {
            rhs[j][i] = ((u_star[j][i] - u_star[j][i-1]) / dx + (v_star[j][i] - v_star[j-1][i]) / dy) / dt;
            mean_rhs += rhs[j][i];
            count++;
        }
    }
    mean_rhs /= count;
    for (int j=1; j<Ny+1; j++){
        for (int i=1; i<Nx+1; i++){
            rhs[j][i] -= mean_rhs;
        }
    }
}

void solve_poisson() { // GS 
    for (iter = 0; iter < poisson_max_iter; iter++) {
        max_error = 0.0;
        for (int j = 1; j < Ny+1; j++) {
            for (int i = 1; i < Nx+1; i++) {
                if (!fluid[j][i]) continue;
                double pE = fluid[j][i+1] ? p[j][i+1] : p[j][i]; // if there is fluid to the east, pE is p[j][i+1]; if there is not, set the gradient to 0
                double pW = fluid[j][i-1] ? p[j][i-1] : p[j][i]; // similarly for west
                double pN = fluid[j+1][i] ? p[j+1][i] : p[j][i]; // south
                double pS = fluid[j-1][i] ? p[j-1][i] : p[j][i]; // and north
                double residual = (p[j][i-1]-2*p[j][i]+p[j][i+1])*(dy*dy) + (p[j-1][i]-2*p[j][i]+p[j+1][i])*(dx*dx) - rhs[j][i]*dy*dy*dx*dx;
                double p_old = p[j][i];
                p[j][i] = (
                    (pE + pW) * dy * dy +
                    (pN + pS) * dx * dx -
                    rhs[j][i] * dx * dx * dy * dy
                ) / (2.0 * (dx * dx + dy * dy));
                double err = fabs(residual);
                if (err > max_error) max_error = err;
            }
        }
        apply_pressure_bc();

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
    snprintf(filename, sizeof(filename), "output_%04d.csv", step);
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        perror("File opening failed");
        return;
    }

    //fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "U-channel incompressible flow\n");
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
    define_fluid_region();
    for (int step = 0; step <= max_iter; step++) {
        apply_boundary_conditions();
        predictor();
        find_divergence();
        solve_poisson();
        corrector();
        apply_boundary_conditions();

        printf("Step=%d iter=%d maxErr=%e u=%.4f v=%.4f p=%.4f\n",step, iter, max_error, u[cw/2][Nx/2], v[cw/2][Nx/2],p[cw/2][Nx/2]);
        if (step % 500 == 0) write_vtk(step);
    }

    printf("Simulation complete.\n");
    return 0;
}