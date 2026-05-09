#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>

#define Nx 30
#define Ny 30
#define Lx 2.0
#define Ly 2.0
#define W 0.5 // channel width
#define Re 100.0
#define rho 1000 // kg/m3
#define mu 0.0010005 // Pa*s
#define dt 5e-4
#define max_iter 5000
#define poisson_tol 1e-5
#define poisson_max_iter 50000
#define save_every 500 // step
double max_error;
double V_in = -3.0; //-Re * mu / (rho * W); // compute inlet velocity from Re
int iter;
const double dx = Lx / Nx;
const double dy = Ly / Ny;
int cw = W/dx;
int fluid[Ny+2][Nx+2] = {0};
double u[Ny+2][Nx+1] = {0};       // u(i+1/2,j)
double ufluid[Ny+2][Nx+1] = {0};
double v[Ny+1][Nx+2] = {0};       // v(i,j+1/2)
double vfluid[Ny+1][Nx+2] = {0}; 
double p[Ny+2][Nx+2] = {0};       // p(i,j)
double u_star[Ny+2][Nx+1] = {0};
double v_star[Ny+1][Nx+2] = {0};
double rhs[Ny+2][Nx+2] = {0};

void init_fluid_p()
{   
   for (int i = 1; i <= Nx; i++){
        for (int j = 1; j <= Ny; j++){
            fluid[j][i] = 1;
        }
    }
    // mask (middle)
    for (int i = cw+1; i <= Nx-cw; i++){ 
        for (int j = cw+1; j <= Ny+1; j++){
            fluid[j][i] = 0;
        }
    }
    // mask (top and bottom ghost nodes are not fluid, except outlet)
    for (int i = 0; i <= Nx+1; i++){
        //if (i <= Nx-cw) fluid[Ny+1][i] = 0;
        fluid[0][i] = 0;
    }
    // mask (left and right ghost nodes are not fluid)
    for (int j = 0; j <= Ny+1; j++){
        fluid[j][0] = 0;
        fluid[j][Nx+1] = 0;
    }
    for (int i = 0; i<=Nx+1; i++){
        for (int j=0; j<=Ny+1; j++){

            printf("%d %d %d %d\n",i, j, fluid[j][i],cw);
        }
    }
    
}
void init_fluid_u()
{
   for (int i = 0; i <= Nx; i++){
        for (int j=0; j<=Ny+1; j++){
            if (fluid[j][i+1] == 1 && fluid[j][i] == 1) ufluid[j][i] = 1;
        }
    }
}
void init_fluid_v()
{
   for (int i = 0; i <= Nx+1; i++){
        for (int j=0; j<=Ny; j++){
            if (fluid[j+1][i] == 1 && fluid[j][i] == 1) vfluid[j][i] = 1;
        }
    }
}


void apply_star_bcs(){
 // --- For u --- //

    // // Outer left and right walls
    // for (int j = 1; j <= Ny; j++) {
    //     u_star[j][1] = 0.0;
    //     u_star[j][Nx] = 0.0;
    // }

    // Inner left and right walls
    for (int j = cw; j <= Ny; j++) {
        u_star[j][cw] = 0.0;
        u_star[j][Nx-cw+1] = 0.0;
    }

    // // --- For v --- //

    // // Outer bottom wall
    // for (int i = 1; i <= Nx; i++) {
    //     v_star[1][i] = 0.0;
    // }

    // Inner bottom wall
    for (int i = cw; i <= Nx+1-cw; i++) {
        v_star[cw][i] = 0.0;
    }
    // inlet
    for (int i = 1; i <= cw; i++) {
        v_star[Ny][i] = V_in;
    }
}


void apply_boundary_conditions() {

    // --- For u --- //

    // Outer left and right walls
    for (int j = 0; j <= Ny+1; j++) {
        u[j][0] = 0.0;
        u[j][Nx] = 0.0;
    }

    // Outer bottom wall (ghost cell reflection)
    for (int i = 1; i <= Nx-1; i++) {
        u[0][i] = -u[1][i];
    }

    // Inner bottom wall
    for (int i = cw; i <= Nx+1-cw; i++) {
        u[cw+1][i] = -u[cw][i];
    }

    // Inner left and right walls
    for (int j = cw+1; j < Ny+1; j++) {
        u[j][cw] = 0.0;
        u[j][Nx-cw+1] = 0.0;
    }


    // --- For v --- //

    // Outer left wall (ghost)
    for (int j = 0; j <= Ny+1; j++) {
        v[j][0] = -v[j][1];

    // Outer right wall (ghost)
        v[j][Nx+1] = -v[j][Nx];
    }

    // Outer bottom wall
    for (int i = 0; i <= Nx+1; i++) {
        v[0][i] = 0.0;
    }

    // Inner bottom wall
    for (int i = cw+1; i <= Nx-cw+1; i++) {
        v[cw][i] = 0.0;
    }

    // Inner left wall (ghost)
    for (int j = cw; j <= Ny; j++) { 
        v[j][cw+1] = -v[j][cw];

    // Inner right wall (ghost) 
        v[j][Nx+1-cw] = -v[j][Nx-cw+1];
    }


    // --- Inlet --- //
    for (int i = 1; i <= cw; i++) {
        u[Ny+1][i] = u[Ny][i];
        v[Ny][i] = V_in;
    }

    // --- Outlet --- //
    for (int i = Nx - cw + 1; i <= Nx+1; i++) {
        u[Ny+1][i] = u[Ny][i];
    }
    for (int i = Nx - cw; i <= Nx; i++){
        v[Ny][i] = v[Ny+1][i];
    }
}

void apply_pressure_bc(){

    /* --- Inlet --- */
    for (int i = 1; i <= cw; i++) {
        p[Ny+1][i] = p[Ny][i];
    }
    /* --- Outlet --- */ 
    // applied pressure pin to single point at outlet in the poisson function
    // /* for U-channel outlet at top of right branch */
    // for (int i = Nx - cw + 1; i <= Nx; i++)
    // {
    //     p[Ny][i] = 0.0;
    // }
    p[Ny][Nx-cw/2] = 0.0;
}

double max(double x, double y) {
    return (x > y) ? x : y;
}

double min(double x, double y) {
    return (x < y) ? x : y;
}

void predictor() {
    for (int j = 0; j <= Ny+1; j++) {
        for (int i = 0; i <= Nx; i++) {
            if (ufluid[j][i]==0) continue;
            double ue_i = u[j][i];
            if (ufluid[j][i+1]) ue_i = 0.5*(u[j][i]+u[j][i+1]);
            double uw_i = u[j][i];
            if (ufluid[j][i-1]) uw_i = 0.5*(u[j][i]+u[j][i-1]);
            double vn_i = 0.0;
            double vs_i = 0.0;
            if ((vfluid[j][i] 
                + vfluid[j+1][i] 
                + vfluid[j][i+1] 
                + vfluid[j+1][i+1]) == 0) {vn_i =0.0;}
            else {
                vn_i = (
                    vfluid[j][i]*v[j][i] 
                    + vfluid[j+1][i]*v[j+1][i]
                    + vfluid[j][i+1]*v[j][i+1]
                    + vfluid[j+1][i+1]*v[j+1][i+1]) / 
                    (vfluid[j][i] 
                    + vfluid[j+1][i]
                    + vfluid[j][i+1]
                    + vfluid[j+1][i+1]);
                }
            
            if ((vfluid[j][i] 
                + vfluid[j-1][i] 
                + vfluid[j][i+1] 
                + vfluid[j-1][i+1]) == 0) {vs_i =0.0;}
            else {
                vs_i = (
                vfluid[j][i]*v[j][i] 
                + vfluid[j-1][i]*v[j-1][i]
                + vfluid[j][i+1]*v[j][i+1]
                + vfluid[j-1][i+1]*v[j-1][i+1]) / 
                (vfluid[j][i] 
                + vfluid[j-1][i]
                + vfluid[j][i+1]
                + vfluid[j-1][i+1]);
            }
            
            double du2dx = (max(0,ue_i)*u[j][i] + min(0,ue_i)*u[j][i+1] - max(0,uw_i)*u[j][i-1] - min(0,uw_i)*u[j][i])/dx;
            double duvdy = (max(0,vn_i)*u[j][i] + min(0,vn_i)*u[j+1][i] - max(0,vs_i)*u[j-1][i] - min(0,vs_i)*u[j][i])/dy;
            double d2udx2 = 0.0;
            double d2udy2 = 0.0;
            if (ufluid[j][i+1] && ufluid[j][i-1]) d2udx2 = (u[j][i+1] - 2*u[j][i] + u[j][i-1]) / (dx * dx);
            if (ufluid[j+1][i] && ufluid[j-1][i]) d2udy2 = (u[j+1][i] - 2*u[j][i] + u[j-1][i]) / (dy * dy);
            
            u_star[j][i] = u[j][i] + dt * (-du2dx - duvdy + 2.0*(d2udx2 + d2udy2)); //
        }
    }

    for (int j = 0; j <= Ny; j++) {
        for (int i = 0; i <= Nx+1; i++) {
            if (vfluid[j][i]==0) continue;
            double vn_j = v[j][i];
            if (vfluid[j+1][i]) vn_j = 0.5*(v[j][i]+ v[j+1][i]);
            double vs_j = v[j][i];
            if (vfluid[j-1][i]) vs_j = 0.5*(v[j][i]+v[j-1][i]);
            double ue_j = 0;
            double uw_j = 0;
            if ((ufluid[j][i] 
                + ufluid[j+1][i]
                + ufluid[j][i+1]
                + ufluid[j+1][i+1])==0) {ue_j = 0.0;}
            else {
                ue_j = (
                    ufluid[j][i]*u[j][i] 
                    + ufluid[j][i+1]*u[j][i+1]
                    + ufluid[j+1][i]*u[j+1][i]
                    + ufluid[j+1][i+1]*u[j+1][i+1]) / 
                    (ufluid[j][i] 
                    + ufluid[j+1][i]
                    + ufluid[j][i+1]
                    + ufluid[j+1][i+1]);
            }
            if ((ufluid[j][i] 
                + ufluid[j][i-1]
                + ufluid[j+1][i]
                + ufluid[j+1][i-1])==0) {uw_j = 0;}
            else {
                uw_j = 
                (ufluid[j][i]*u[j][i] 
                + ufluid[j][i-1]*u[j][i-1]
                + ufluid[j+1][i]*u[j+1][i]
                + ufluid[j+1][i-1]*u[j+1][i-1]) / 
                (ufluid[j][i] 
                + ufluid[j][i-1]
                + ufluid[j+1][i]
                + ufluid[j+1][i-1]);
            }
                
            double dv2dy = (max(0,vn_j)*v[j][i] + min(0,vn_j)*v[j+1][i] - max(0,vs_j)*v[j-1][i] - min(0,vs_j)*v[j][i])/dy;
            double duvdx = (max(0,ue_j)*v[j][i] + min(0,ue_j)*v[j][i+1] - max(0,uw_j)*v[j][i-1] - min(0,uw_j)*v[j][i])/dx;
            
            double d2vdx2 = 0.0;
            double d2vdy2 = 0.0;
            if (vfluid[j+1][i] && vfluid[j-1][i]) d2vdx2 = (v[j][i+1] - 2*v[j][i] + v[j][i-1]) / (dx * dx);
            if (vfluid[j][i+1] && vfluid[j][i-1]) d2vdy2 = (v[j+1][i] - 2*v[j][i] + v[j-1][i]) / (dy * dy);
            
         
            v_star[j][i] = v[j][i] + dt * (-dv2dy - duvdx + 2.0*(d2vdx2 + d2vdy2)); // 
            
        } 
    }
}

void find_divergence() { 
    int count = 0;
    double sum = 0.0;
    
    for (int j = 1; j <= Ny; j++) {
        for (int i = 1; i <= Nx; i++) {
            // skip inlet/outlet faces since v_star and u_star do not get bcs for them
            //if (j == Ny && i <= cw) continue;        // inlet column
            //if (j == Ny && i >= Nx-cw) continue;     // outlet column
            if (fluid[j][i]==0) continue;
            rhs[j][i] = ((u_star[j][i] - u_star[j][i-1]) / dx + (v_star[j][i] - v_star[j-1][i]) / dy) / dt;
            sum += rhs[j][i]; 
            count++;
            //printf("RHS %f\n",rhs[j][i]);
        }
    }
    // double rhs_sum = 0.0;
    // double rhs_max = 0.0;
    // for (int j=1; j<=Ny; j++)
    //     for (int i=1; i<=Nx; i++)
    //         if (fluid[j][i]) {
    //             rhs_sum += rhs[j][i];
    //             if (fabs(rhs[j][i]) > rhs_max) rhs_max = fabs(rhs[j][i]);
    //         }
    printf("RHS sum=%e \n", sum);
    // for (int j=1; j<=Ny; j++) { // subtract mean of rhs from rhs
    //     for (int i=1; i<=Nx; i++){
    //         if (fluid[j][i]) rhs[j][i] -= sum/count;
    //     }
    // }
}

void solve_poisson() {
    double omega = 1.0;
    for (iter = 0; iter < poisson_max_iter; iter++) {
        max_error = 0.0;
        for (int j = 1; j < Ny+1; j++) {
            for (int i = 1; i < Nx+1; i++) {
                if (fluid[j][i]==0) continue;
                double pE = fluid[j][i+1] ? p[j][i+1] : p[j][i]; // if there is fluid to the east, pE is p[j][i+1]; if there is not, set the gradient to 0
                double pW = fluid[j][i-1] ? p[j][i-1] : p[j][i]; // similarly for west
                double pN = fluid[j+1][i] ? p[j+1][i] : p[j][i]; // south
                double pS = fluid[j-1][i] ? p[j-1][i] : p[j][i]; // and north
                //double p_old = p[j][i];
                // double fluid_neighbors_x = fluid[j][i+1] + fluid[j][i-1];
                // double fluid_neighbors_y = fluid[j+1][i] + fluid[j-1][i];
                
                // if ((fluid_neighbors_x * (dx * dx) + fluid_neighbors_y * (dy * dy)) < 1e-6) {
                //     printf("error: zero fluid neighbors");
                //     continue;
                // }

                double residual = ((pW-2.0*p[j][i]+pE)*(dy*dy) + (pS-2.0*p[j][i]+pN)*(dx*dx) - rhs[j][i]*dy*dy*dx*dx)
                 / (2.0 * (dx*dx + dy*dy));
                double p_old = p[j][i];
                double p_new = (
                    (pE + pW) * (dy * dy) +
                    (pN + pS) * (dx * dx) -
                    rhs[j][i] * (dy * dy) * (dx * dx)
                ) / (2.0 * ((dx * dx + dy * dy)));
                p[j][i] = (1.0-omega)*p[j][i] + omega*p_new;
                double err = fabs(residual);
                if (err > max_error) max_error = err;
                if (j==3 && i==3 && iter % 10000 == 0) {
                //printf("iter=%d p=%e pE=%e pW=%e pN=%e pS=%e rhs=%e\n",
                //iter, p[j][i], pE, pW, pN, pS, rhs[j][i]);
                
                }
            }
            //p[Ny][Nx-cw+1] = 0.0;
        }
       
        apply_pressure_bc();
        
        if (max_error < poisson_tol) break;
    }
}


void corrector(){
    
    for (int j = 1; j < Ny+1; j++) {
        for (int i = 1; i < Nx; i++) {
            if (!ufluid[j][i]) continue;
            u[j][i] = u_star[j][i] - dt * (p[j][i+1] - p[j][i]) / dx;
        }
    }

    for (int j = 1; j < Ny; j++) {
        for (int i = 1; i < Nx+1; i++) {
            if (!vfluid[j][i]) continue;
            v[j][i] = v_star[j][i] - dt * (p[j+1][i] - p[j][i]) / dy;
        }
    }
}

void find_flux() {
    double inflow = 0.0;
    double outflow = 0.0;

    for (int i=1; i<=cw; i++)
        inflow += v[Ny][i] * dx;

    for (int i=Nx-cw; i<=Nx; i++)
        outflow += v[Ny][i] * dx;

    printf("in=%e out=%e\n", inflow, outflow);
}

void write_csv(int step) {
    char filename[64];
    snprintf(filename, sizeof(filename), "output_%04d.csv", step);
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        perror("File opening failed");
        return;
    }

    fprintf(fp, "Test flow problems\n");
    fprintf(fp, "DIMENSIONS, %d, %d, 1, , \n", Nx, Ny);
    //fprintf(fp, "ORIGIN 0 0 0\n");
    fprintf(fp, "SPACING, %f, %f, 1, , \n", dx, dy);
    fprintf(fp, "POINT_DATA, %d\n, , , , \n", Nx * Ny);

    // Print to columns
    // header
    fprintf(fp, "Pressure, Velocity Magnitude, RHS, U-component, V-component, Fluid\n"); // Vorticity (z-comp of curl)
    for (int j = 1; j <= Ny; j++) {
        for (int i = 1; i <= Nx; i++) {
            double uc = 0.5 * (u[j][i] + u[j][i-1]);
            double vc = 0.5 * (v[j][i] + v[j-1][i]);
            double dudy = (u[j+1][i] - u[j-1][i]) / (2.0 * dy);
            double dvdx = (v[j][i+1] - v[j][i-1]) / (2.0 * dx);
            double omega = dvdx - dudy;
            fprintf(fp, "%f, %f, %f, %f, %f, %d \n", p[j][i], sqrt(uc * uc + vc * vc), rhs[j][i], uc, vc, fluid[j][i]);
        }
    }

    fclose(fp);
}


int main() {
    init_fluid_p();
    init_fluid_u();
    init_fluid_v();
    for (int step = 0; step <= max_iter; step++) {
        // if(V_in > -1.0) V_in -= 0.2;
        // printf("u at inner bottom wall: %f %f\n", u[cw][cw], u[cw+1][cw]);
        // printf("v at inner bottom wall: %f %f\n", v[cw][cw+1], v[cw+1][cw+1]);
        // printf("Corner check:\n");
        // for (int j = cw-1; j <= cw+3; j++)
        //     for (int i = cw-1; i <= cw+3; i++)
        //         printf("u[%d][%d]=%e v[%d][%d]=%e fluid=%d\n",
        //             j,i,u[j][i],j,i,v[j][i],fluid[j][i]);
        apply_boundary_conditions();
        predictor();
        //apply_star_bcs();
        find_divergence();
        solve_poisson();
        corrector();
        //apply_boundary_conditions();
        //apply_star_bcs();
        find_flux();
        int probe_x = cw/2;
        int probe_y = Ny-2;
        printf("Step=%d iter=%d maxErr=%e u=%.4f v=%.4f p=%.4f, u_star=%.4f, v_star=%.4f \n",step, iter,
            max_error, u[probe_y][probe_x], v[probe_y][probe_x],p[probe_y][probe_x],u_star[probe_y][probe_x],v_star[probe_y][probe_x]);
        if (step % save_every == 0) write_csv(step);
    }
        

    printf("Simulation complete.\n");
    return 0;
}
