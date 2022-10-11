
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;

#include <cblas.h>

#include <mpi.h>

#include <chrono>
using namespace std::chrono;

#include "ReactionDiffusion.h"

ReactionDiffusion::ReactionDiffusion(double dtf, double integration_timef, int Nxf, 
            int Nyf, double af, double bf, double mu_1f, 
            double mu_2f, double epsf, int rankf, int sizef) {

    if (sizef >= 49) {
    dt = dtf;
    integration_time = integration_timef;
    Nx = Nxf;
    Ny = Nyf;
    a = af;
    b = bf;
    mu_1 = mu_1f;
    mu_2 = mu_2f;
    eps = epsf;

    rank = rankf;
    size = sizef;


    // make the cartesian communicator
    num_1D = sqrt(size);
    dims = 2;
    sizes = new int[dims];
    sizes[0] = num_1D;
    sizes[1] = num_1D;
    periods = new int[dims];
    periods[0] = 0;
    periods[1] = 1;
    reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, dims, sizes, periods, reorder, &mygrid);

    // get coords for each process
    MPI_Comm_rank(mygrid, &mygrid_rank);
    coords = new int[dims];
    MPI_Cart_coords(mygrid, mygrid_rank, dims, coords);

    // make row and column communicators
    keep = new int[dims];

    keep[0] = 0;
    keep[1] = 1; // keep rows in subgrid
    MPI_Cart_sub(mygrid, keep, &myrowcomm);
    MPI_Comm_rank(myrowcomm, &row_rank);
    MPI_Comm_size(myrowcomm, &row_size);    

    keep[0] = 1; // keep columns in subgrid
    keep[1] = 0;
    MPI_Cart_sub(mygrid, keep, &mycolcomm);
    MPI_Comm_rank(mycolcomm, &col_rank); // first column, second column, etc
    MPI_Comm_size(mycolcomm, &col_size);


    // split up domain in y-direction
    remainder_y = Ny % num_1D;
    Ny_small = (Ny - remainder_y)/num_1D; // if all were this size we would just not have enough
    Ny_big = Ny_small + 1; // if all were this size we would have just too many

    Ny_sizes = new double[num_1D];
    for (int i = 0; i < num_1D; ++i) {
        if (i < remainder_y) {
            Ny_sizes[i] = double(Ny_big);
        }
        else {
            Ny_sizes[i] = double(Ny_small);
        }
    }


    // make Ny_sizes_all
    Ny_sizes_all = new double[size];

    if (rank == 0) {
        for (int i = 0; i < num_1D; ++i) {
            for (int j = 0; j < num_1D; ++j) {
                Ny_sizes_all[i*num_1D+j] = Ny_sizes[j];

            }
        }

    }


    // split up domain in x-direction
    remainder_x = Nx % num_1D;
    Nx_small = (Nx - remainder_x)/num_1D; // if all were this size we would just not have enough
    Nx_big = Nx_small + 1; // if all were this size we would have just too many

    Nx_sizes = new double[num_1D];
    for (int i = 0; i < num_1D; ++i) {
        if (i < remainder_x) {
            Nx_sizes[i] = double(Nx_big);
        }
        else {
            Nx_sizes[i] = double(Nx_small);

        }
    }


    if (coords[1] < remainder_y) {
        Ny_size = Ny_big;
    }
    else {
        Ny_size = Ny_small;
    }

    if (coords[0] < remainder_x) {
        Nx_size = Nx_big;
    }
    else {
        Nx_size = Nx_small;
    }


    // make Nx_sizes_all
    Nx_sizes_all = new double[size];
    if (rank == 0) {
        for (int i = 0; i < num_1D; ++i) {
            for (int j = 0; j < num_1D; ++j) {
                Nx_sizes_all[i*num_1D+j] = Nx_sizes[i];

            }
        }
    }



    // make empty arrays
    uuu = new double[Nx_size*Ny_size];
    vvv = new double[Nx_size*Ny_size];

    uuu_work = new double[Nx_size*2+Ny_size*2];
    vvv_work = new double[Nx_size*2+Ny_size*2];
    uuu_work_send_horiz = new double[Ny_size*2];
    vvv_work_send_horiz = new double[Ny_size*2];


    // initialise more arrays
    uuu_main = new double[Nx*Ny];
    vvv_main = new double[Ny*Nx];
    u_lap_x = new double[Nx_size*Ny_size];
    v_lap_x = new double[Nx_size*Ny_size];
    for (int i = 0; i < Nx_size*Ny_size; ++i) v_lap_x[i] = 0;
    uuu_top_bottom = new double[Nx*2];
    vvv_top_bottom = new double[Nx*2];
    u_lap_y = new double[Nx_size*Ny_size];
    v_lap_y = new double[Nx_size*Ny_size];
    for (int i = 0; i < Nx_size*Ny_size; ++i) v_lap_y[i] = 0;
    f_1 = new double[Nx_size*Ny_size];
    f_2 = new double[Nx_size*Ny_size];

    b_over_a = bf/af;
    over_a = 1/af;

    }

    else {
    dt = dtf;
    integration_time = integration_timef;
    Nx = Nxf;
    Ny = Nyf;
    a = af;
    b = bf;
    mu_1 = mu_1f;
    mu_2 = mu_2f;
    eps = epsf;

    rank = rankf;
    size = sizef;


    // split up domain in y-direction
    int remainder_y = Ny % size;
    int Ny_small = (Ny - remainder_y)/size; // if all were this size we would just not have enough
    int Ny_big = Ny_small + 1; // if all were this size we would have just too many

    // int Ny_size;
    if (rank < remainder_y) {
        Ny_size = Ny_big;
    }
    else {
        Ny_size = Ny_small;
    }


    Ny_sizes = new double[size];
    y_recvcounts = new int[size];
    for (int i = 0; i < size; ++i) {
        if (i < remainder_y) {
            Ny_sizes[i] = double(Ny_big);
            y_recvcounts[i] = Nx*Ny_big;
        }
        else {
            Ny_sizes[i] = double(Ny_small);
            y_recvcounts[i] = Nx*Ny_small;
        }
    }


    // make empty arrays
    uuu_across = new double[Ny_size*Nx];
    vvv_across = new double[Ny_size*Nx];


    // initialise more arrays
    uuu_main = new double[Nx*Ny];
    vvv_main = new double[Ny*Nx];
    u_lap_x = new double[Ny_size*Nx];
    v_lap_x = new double[Ny_size*Nx];
    for (int i = 0; i < Ny_size*Nx; ++i) v_lap_x[i] = 0;
    uuu_top_bottom = new double[Nx*2];
    vvv_top_bottom = new double[Nx*2];
    u_lap_y = new double[Ny_size*Nx];
    v_lap_y = new double[Ny_size*Nx];
    for (int i = 0; i < Ny_size*Nx; ++i) v_lap_y[i] = 0;
    f_1 = new double[Ny_size*Nx];
    f_2 = new double[Ny_size*Nx];

    b_over_a = bf/af;
    over_a = 1/af;
    }



}


void ReactionDiffusion::SetInitialConditions() {

    if (size >= 49) {
    // initialise the uuu and vvv arrays
    count_j = 0;
    count_i = 0;
    for (int j = int(cblas_dasum(row_rank, Ny_sizes, 1)); j < int(cblas_dasum(row_rank+1, Ny_sizes, 1)); ++j) {
        for (int i = int(cblas_dasum(col_rank, Nx_sizes, 1)); i < int(cblas_dasum(col_rank+1, Nx_sizes, 1)); ++i) {
            
            // initialise uuu
            if (j < int(floor(Ny/2)) ) {
                uuu[count_j*int(Nx_sizes[col_rank]) + count_i] = 1;
            }
            else {
                uuu[count_j*int(Nx_sizes[col_rank]) + count_i] = 0;
            }

            // initialise vvv
            if (i < floor(Nx/2) ) {
                vvv[count_j*int(Nx_sizes[col_rank]) + count_i] = a/2;
            }
            else {
                vvv[count_j*int(Nx_sizes[col_rank]) + count_i] = 0;
            }

            ++count_i;


        }
        count_i = 0;
        ++count_j;
    }

    }

    else {
    // initialise uuu_across and vvv_across
    int count = 0;
    for (int i = 0; i < Nx; ++i) {
        for (int j = int(cblas_dasum(rank, Ny_sizes, 1)); j < int(cblas_dasum(rank+1, Ny_sizes, 1)); ++j) {
            // initialise uuu_across
            if (j < int(floor(Ny/2)) ) {
                uuu_across[count*Nx + i] = 1;
                // uuu_across[count*Nx + i] = i*Nx+j;
            }
            else {
                uuu_across[count*Nx + i] = 0;
                // uuu_across[count*Nx + i] = i*Nx+j;
            }

            // initialise vvv_across
            if (i < floor(Nx/2) ) {
                vvv_across[count*Nx + i] = a/2;
            }
            else {
                vvv_across[count*Nx + i] = 0;
            }

            ++count;

        }
        count = 0;
    }
    }

}

void ReactionDiffusion::TimeIntegrate() {

    if (size >= 49) {

    int main_count = 0;
//     int count;
//     int through_rows;
//     int index;
//     int i;
    int num_of_iterations = ceil(integration_time/dt);
    // num_of_iterations = 1;

    while (main_count < num_of_iterations) {
    
        // tags are only needed to distinguish between u and v communication

        // send up
        if (coords[1] > 0) {
            MPI_Isend(uuu, Nx_size, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &u_send_up);
            MPI_Isend(vvv, Nx_size, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &v_send_up);
        }

        // send up being received
        if (coords[1] < num_1D-1) {
            MPI_Irecv(uuu_work, Nx_size, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &u_recv_up);
            MPI_Irecv(vvv_work, Nx_size, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &v_recv_up);
            // MPI_Wait(&u_recv_up, MPI_STATUS_IGNORE);
            // MPI_Wait(&v_recv_up, MPI_STATUS_IGNORE);
        }

        // send down
        if (coords[1] < num_1D-1) {
            MPI_Isend(uuu + (Ny_size-1)*Nx_size, Nx_size, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &u_send_down);
            MPI_Isend(vvv + (Ny_size-1)*Nx_size, Nx_size, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &v_send_down);
        }

        // send down being received
        if (coords[1] > 0) {
            MPI_Irecv(uuu_work + Nx_size, Nx_size, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &u_recv_down);
            MPI_Irecv(vvv_work + Nx_size, Nx_size, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &v_recv_down);
            // MPI_Wait(&u_recv_down, MPI_STATUS_IGNORE);
            // MPI_Wait(&v_recv_down, MPI_STATUS_IGNORE);
        }


        // // find y-direction contribution to laplacian for updating uuu and vvv
        count_j = 0;
        count_i = 0;
        for (int j = int(cblas_dasum(row_rank, Ny_sizes, 1)); j < int(cblas_dasum(row_rank+1, Ny_sizes, 1)); ++j) {
            for (int i = int(cblas_dasum(col_rank, Nx_sizes, 1)); i < int(cblas_dasum(col_rank+1, Nx_sizes, 1)); ++i) {
                
                // deal with top boundary
                if (j == 0) {
                    u_lap_y[count_j*Nx_size + count_i] = -uuu[0*Nx_size + count_i] + uuu[1*Nx_size + count_i];
                    if (mu_2 != 0) v_lap_y[count_j*Nx_size + count_i] = -vvv[0*Nx_size + count_i] + vvv[1*Nx_size + count_i];
                }

                else if (j == Ny-1) {
                    u_lap_y[count_j*Nx_size + count_i] = -uuu[count_j*Nx_size + count_i] + uuu[(count_j-1)*Nx_size + count_i];
                    if (mu_2 != 0) v_lap_y[count_j*Nx_size + count_i] = -vvv[count_j*Nx_size + count_i] + vvv[(count_j-1)*Nx_size + count_i];
                }

                // // main inner part
                // not the top or bottom boundaries of a rank
                else if (count_j != 0 && count_j != Ny_size-1) {
                    u_lap_y[count_j*Nx_size + count_i] = uuu[(count_j+1)*Nx_size + count_i] + uuu[(count_j-1)*Nx_size + count_i];
                    if (mu_2 != 0) v_lap_y[count_j*Nx_size + count_i] = vvv[(count_j+1)*Nx_size + count_i] + vvv[(count_j-1)*Nx_size + count_i];
                }

                // top boundary of a rank
                else if (count_j == 0) {
                    MPI_Wait(&u_recv_down, MPI_STATUS_IGNORE);
                    u_lap_y[count_j*Nx_size + count_i] = uuu[(count_j+1)*Nx_size + count_i] + uuu_work[Nx_size + count_i];

                    if (mu_2 != 0) {
                        MPI_Wait(&v_recv_down, MPI_STATUS_IGNORE);
                        v_lap_y[count_j*Nx_size + count_i] = vvv[(count_j+1)*Nx_size + count_i] + vvv_work[Nx_size + count_i];
                    }

                }

                // bottom boundary of a rank
                else if (count_j == Ny_size-1) {
                    MPI_Wait(&u_recv_up, MPI_STATUS_IGNORE);
                    u_lap_y[count_j*Nx_size + count_i] = uuu_work[count_i] + uuu[(count_j-1)*Nx_size + count_i];

                    if (mu_2 != 0) {
                        MPI_Wait(&v_recv_up, MPI_STATUS_IGNORE);
                        v_lap_y[count_j*Nx_size + count_i] = vvv_work[count_i] + vvv[(count_j-1)*Nx_size + count_i];
                    }

                }

                ++count_i;


            }
            count_i = 0;
            ++count_j;
        }


        // send to the right
        // send and recv   | coords[0] is going from column to column
        if (coords[0] < num_1D-1) {
            for (int i = 0; i < Ny_size; ++i) {
                uuu_work_send_horiz[Ny_size + i] = uuu[i*Nx_size + Nx_size-1];
                vvv_work_send_horiz[Ny_size + i] = vvv[i*Nx_size + Nx_size-1];
            }
            MPI_Isend(uuu_work_send_horiz + Ny_size, Ny_size, MPI_DOUBLE, rank + num_1D, 0, MPI_COMM_WORLD, &u_send_right);
            MPI_Isend(vvv_work_send_horiz + Ny_size, Ny_size, MPI_DOUBLE, rank + num_1D, 1, MPI_COMM_WORLD, &v_send_right);
        }

        // send to the right being received
        if (coords[0] > 0) {
            MPI_Irecv(uuu_work + Nx_size*2 + Ny_size, Ny_size, MPI_DOUBLE, rank - num_1D, 0, MPI_COMM_WORLD, &u_recv_right);
            // MPI_Wait(u_recv_right, MPI_STATUS_IGNORE);

            MPI_Irecv(vvv_work + Nx_size*2 + Ny_size, Ny_size, MPI_DOUBLE, rank - num_1D, 1, MPI_COMM_WORLD, &v_recv_right);
        }

        // send to the left
        // send and recv   | coords[0] is going from column to column
        if (coords[0] > 0) {
            for (int i = 0; i < Ny_size; ++i) {
                uuu_work_send_horiz[i] = uuu[i*Nx_size];
                vvv_work_send_horiz[i] = vvv[i*Nx_size];
            }
            MPI_Isend(uuu_work_send_horiz, Ny_size, MPI_DOUBLE, rank - num_1D, 0, MPI_COMM_WORLD, &u_send_left);
            MPI_Isend(vvv_work_send_horiz, Ny_size, MPI_DOUBLE, rank - num_1D, 1, MPI_COMM_WORLD, &v_send_left);
        }

        // send to the left being received
        if (coords[0] < num_1D-1) {
            MPI_Irecv(uuu_work + Nx_size*2, Ny_size, MPI_DOUBLE, rank + num_1D, 0, MPI_COMM_WORLD, &u_recv_left);
            MPI_Irecv(vvv_work + Nx_size*2, Ny_size, MPI_DOUBLE, rank + num_1D, 1, MPI_COMM_WORLD, &v_recv_left);
        }



        // find laplacian and time integrate
        // // find x-direction contribution to laplacian for updating uuu and vvv
        count_j = 0;
        count_i = 0;
        for (int j = int(cblas_dasum(row_rank, Ny_sizes, 1)); j < int(cblas_dasum(row_rank+1, Ny_sizes, 1)); ++j) {
            for (int i = int(cblas_dasum(col_rank, Nx_sizes, 1)); i < int(cblas_dasum(col_rank+1, Nx_sizes, 1)); ++i) {
                
                // deal with left boundary
                if (i == 0) {
                    u_lap_x[count_j*Nx_size + count_i] = uuu[count_j*Nx_size + 1] - uuu[count_j*Nx_size + 0];
                    if (mu_2 != 0) v_lap_x[count_j*Nx_size + count_i] = vvv[count_j*Nx_size + 1] - vvv[count_j*Nx_size + 0];
                }

                // deal with right boundary
                else if (i == Nx-1) {
                    u_lap_x[count_j*Nx_size + count_i] = -uuu[count_j*Nx_size + count_i] + uuu[count_j*Nx_size + (count_i-1)];
                    if (mu_2 != 0) v_lap_x[count_j*Nx_size + count_i] = -vvv[count_j*Nx_size + count_i] + vvv[count_j*Nx_size + (count_i-1)];
                }

                // // main inner part
                // not the left or right boundaries of a rank
                else if (count_i != 0 && count_i != Nx_size-1) {
                    u_lap_x[count_j*Nx_size + count_i] = uuu[count_j*Nx_size + (count_i+1)] + uuu[count_j*Nx_size + (count_i-1)];
                    if (mu_2 != 0) v_lap_x[count_j*Nx_size + count_i] = vvv[count_j*Nx_size + (count_i+1)] + vvv[count_j*Nx_size + (count_i-1)];
                }

                // left boundary of a rank
                else if (count_i == 0) {
                    MPI_Wait(&u_recv_right, MPI_STATUS_IGNORE);
                    u_lap_x[count_j*Nx_size + count_i] = uuu[count_j*Nx_size + (count_i+1)] + uuu_work[Nx_size*2 + Ny_size + count_j];
                    
                    if (mu_2 != 0) {
                        MPI_Wait(&v_recv_right, MPI_STATUS_IGNORE);
                        v_lap_x[count_j*Nx_size + count_i] = vvv[count_j*Nx_size + (count_i+1)] + vvv_work[Nx_size*2 + Ny_size + count_j];
                    }

                }

                // right boundary of a rank
                else if (count_i == Nx_size-1) {
                    MPI_Wait(&u_recv_left, MPI_STATUS_IGNORE);
                    u_lap_x[count_j*Nx_size + count_i] = uuu_work[Nx_size*2 + count_j] + uuu[count_j*Nx_size + (count_i-1)];
                    
                    if (mu_2 != 0) {
                        MPI_Wait(&v_recv_left, MPI_STATUS_IGNORE);
                        v_lap_x[count_j*Nx_size + count_i] = vvv_work[Nx_size*2 + count_j] + vvv[count_j*Nx_size + (count_i-1)];
                    }

                }

                ++count_i;


            }
            count_i = 0;
            ++count_j;
        }


        // // calculate f_1 and f_2
        count_j = 0;
        count_i = 0;
        for (int j = int(cblas_dasum(row_rank, Ny_sizes, 1)); j < int(cblas_dasum(row_rank+1, Ny_sizes, 1)); ++j) {
            for (int i = int(cblas_dasum(col_rank, Nx_sizes, 1)); i < int(cblas_dasum(col_rank+1, Nx_sizes, 1)); ++i) {
                
                f_1[count_j*Nx_size + count_i] = eps * uuu[count_j*Nx_size + count_i] * (1-uuu[count_j*Nx_size + count_i]) * ( uuu[count_j*Nx_size + count_i] - ( vvv[count_j*Nx_size + count_i]*over_a + b_over_a ) );
                f_2[count_j*Nx_size + count_i] = uuu[count_j*Nx_size + count_i]*uuu[count_j*Nx_size + count_i]*uuu[count_j*Nx_size + count_i] - vvv[count_j*Nx_size + count_i];
                

                ++count_i;
            }
            count_i = 0;
            ++count_j;
        }


        // // update uuu with the value of u^{n+1}
        count_j = 0;
        count_i = 0;
        for (int j = int(cblas_dasum(row_rank, Ny_sizes, 1)); j < int(cblas_dasum(row_rank+1, Ny_sizes, 1)); ++j) {
            for (int i = int(cblas_dasum(col_rank, Nx_sizes, 1)); i < int(cblas_dasum(col_rank+1, Nx_sizes, 1)); ++i) {
                
                // deal with corners i.e. don't minus u_{i, j} at all
                if ((i==0 && j==0) || (i==Nx-1 && j==0) || (i==0 && j==Ny-1) || (i==Nx-1 && j==Ny-1)) {
                    uuu[count_j*Nx_size + count_i] = dt * ( mu_1 * (u_lap_x[count_j*Nx_size + count_i] + u_lap_y[count_j*Nx_size + count_i] ) + f_1[count_j*Nx_size + count_i] ) + uuu[count_j*Nx_size + count_i];
                }

                // the below block runs if we are not at a corner but we are at a boundary and so we need to: minus just 2 times u_{i, j}
                else if (i == 0 || i == Nx-1 || j == 0 || j == Ny-1) {
                    uuu[count_j*Nx_size + count_i] = dt * ( mu_1 * (u_lap_x[count_j*Nx_size + count_i] + u_lap_y[count_j*Nx_size + count_i] - 2 * uuu[count_j*Nx_size + count_i]) + f_1[count_j*Nx_size + count_i] ) + uuu[count_j*Nx_size + count_i];
                }
                else {
                    uuu[count_j*Nx_size + count_i] = dt * ( mu_1 * (u_lap_x[count_j*Nx_size + count_i] + u_lap_y[count_j*Nx_size + count_i] - 4 * uuu[count_j*Nx_size + count_i]) + f_1[count_j*Nx_size + count_i] ) + uuu[count_j*Nx_size + count_i];
                }

                ++count_i;
            }
            count_i = 0;
            ++count_j;
        }


        // // update vvv with the value of v^{n+1}
        count_j = 0;
        count_i = 0;
        for (int j = int(cblas_dasum(row_rank, Ny_sizes, 1)); j < int(cblas_dasum(row_rank+1, Ny_sizes, 1)); ++j) {
            for (int i = int(cblas_dasum(col_rank, Nx_sizes, 1)); i < int(cblas_dasum(col_rank+1, Nx_sizes, 1)); ++i) {
                
                // deal with corners i.e. don't minus u_{i, j} at all
                if ((i==0 && j==0) || (i==Nx-1 && j==0) || (i==0 && j==Ny-1) || (i==Nx-1 && j==Ny-1)) {
                    vvv[count_j*Nx_size + count_i] = dt * ( mu_2 * (v_lap_x[count_j*Nx_size + count_i] + v_lap_y[count_j*Nx_size + count_i] ) + f_2[count_j*Nx_size + count_i] ) + vvv[count_j*Nx_size + count_i];
                }

                // the below block runs if we are not at a corner but we are at a boundary and so we need to: minus just 2 times u_{i, j}
                else if (i == 0 || i == Nx-1 || j == 0 || j == Ny-1) {
                    vvv[count_j*Nx_size + count_i] = dt * ( mu_2 * (v_lap_x[count_j*Nx_size + count_i] + v_lap_y[count_j*Nx_size + count_i] - 2 * vvv[count_j*Nx_size + count_i]) + f_2[count_j*Nx_size + count_i] ) + vvv[count_j*Nx_size + count_i];
                }
                else {
                    vvv[count_j*Nx_size + count_i] = dt * ( mu_2 * (v_lap_x[count_j*Nx_size + count_i] + v_lap_y[count_j*Nx_size + count_i] - 4 * vvv[count_j*Nx_size + count_i]) + f_2[count_j*Nx_size + count_i] ) + vvv[count_j*Nx_size + count_i];
                }    

                ++count_i;
            }
            count_i = 0;
            ++count_j;
        }


        ++main_count;
    }


    // find the indices of where the uuu arrays are located within uuu_main
    int indices[size]; // 0, Nx_sizes[0], Nx_sizes[0] + Nx_sizes[1], etc;
                          // 0+Nx*Ny_sizes[0], Nx_sizes[0]+Nx*Ny_sizes[0], etc;
                          // etc
        // but use column storage

    for (int i = 0; i < num_1D; ++i) {
        for (int j = 0; j < num_1D; ++j) {
            indices[i*num_1D + j] = int( cblas_dasum(i, Nx_sizes, 1) + Nx*cblas_dasum(j, Ny_sizes, 1) );
        }
    }


    int Nx_sizes_all_int[size];
    // int Nx_sizes_all_int_1less[size];
    for (int i = 0; i < size; ++i) {
        Nx_sizes_all_int[i] = int(Nx_sizes_all[i]);
        // Nx_sizes_all_int_1less[i] = int(Nx_sizes_all[i]) - 1;
        
    }

    // send the uuu and vvv arrays from each rank to uuu_main and vvv_main on rank 0
    for (int num_times = 0; num_times < Ny_small; ++num_times) {
        MPI_Gatherv(uuu + Nx_size*num_times, Nx_size, MPI_DOUBLE, uuu_main, Nx_sizes_all_int, indices, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gatherv(vvv + Nx_size*num_times, Nx_size, MPI_DOUBLE, vvv_main, Nx_sizes_all_int, indices, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        for (int i = 0; i < size; ++i) indices[i] += (Nx);        
    }


    for (int i = 0; i < size; ++i) {
        if (int(Ny_sizes_all[i]) == Ny_big) {
            Nx_sizes_all_int[i] = Nx_sizes_all_int[i];
        }
        else {
            Nx_sizes_all_int[i] = 0;
        }
    }

    MPI_Gatherv(uuu + Nx_size*Ny_small - 1*(1-(Ny_size-Ny_small)), (Nx_size)*(Ny_size-Ny_small), MPI_DOUBLE, uuu_main, Nx_sizes_all_int, indices, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(vvv + Nx_size*Ny_small - 1*(1-(Ny_size-Ny_small)), (Nx_size)*(Ny_size-Ny_small), MPI_DOUBLE, vvv_main, Nx_sizes_all_int, indices, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (1 && rank == 0) {
        // write uuu and vvv to file
        ofstream myfile;
        myfile.open ("output.txt");

        for (int j = Ny-1; j >=0; --j) { // going from bottom to top i.e. from y=0 to y=height
            for (int i = 0; i < Nx; ++i) { // fastest changing dimension is x
                myfile << i << " " << (Ny-1-j) << " " << uuu_main[j*Nx + i] << " " << vvv_main[j*Nx + i]  << "\n";
            }
            myfile << "\n";
        }

        myfile.close();
    }



    }

    else {
    int main_count = 0;
    int count;
    int through_rows;
    int index;
    int i;
    int num_of_iterations = ceil(integration_time/dt);

    while (main_count < num_of_iterations) {


        // find laplacian and time integrate
        // // find x-direction contribution to laplacian for updating uuu and vvv
        count = 0;
        for (int j = int(cblas_dasum(rank, Ny_sizes, 1)); j < int(cblas_dasum(rank+1, Ny_sizes, 1)); ++j) {

            // deal with left boundary
            through_rows = count*Nx;
            u_lap_x[through_rows] = uuu_across[through_rows + 1] - uuu_across[through_rows + 0];
            if (mu_2 != 0) v_lap_x[through_rows] = vvv_across[through_rows + 1] - vvv_across[through_rows + 0];

            // go across in i for the main making of lap_x
            for (int i = 1; i < Nx-1; ++i) {
                u_lap_x[through_rows + i] = uuu_across[through_rows + (i+1)] + uuu_across[through_rows + (i-1)];
                if (mu_2 != 0) v_lap_x[through_rows + i] = vvv_across[through_rows + (i+1)] + vvv_across[through_rows + (i-1)];
            }

            // deal with right boundary
            i = Nx-1;
            u_lap_x[through_rows + i] = -uuu_across[through_rows + i] + uuu_across[through_rows + (Nx-2)];
            if (mu_2 != 0) v_lap_x[through_rows + i] = -vvv_across[through_rows + i] + vvv_across[through_rows + (Nx-2)];            


            ++count;
        }


        // // find y-direction contribution to laplacian for updating uuu and vvv

        // all ranks except last: send last referenced row of uuu (i.e. last stored column) to next rank
        for (int i = 0; i < size-1; ++i) {
            // if (rank == 0) cout << "\n in here" << endl;

            // send
            if (rank == i) {
                MPI_Isend(uuu_across + ( (Ny_size-1)*Nx ), Nx, MPI_DOUBLE, i+1, 0, MPI_COMM_WORLD, &u_send_last_req);
                if (mu_2 != 0) MPI_Isend(vvv_across + ( (Ny_size-1)*Nx ), Nx, MPI_DOUBLE, i+1, 1, MPI_COMM_WORLD, &v_send_last_req);
            }

            // recv
            if (rank == i+1) {
                MPI_Irecv(uuu_top_bottom, Nx, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &u_recv_last_req);
                if (mu_2 != 0) MPI_Irecv(vvv_top_bottom, Nx, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &v_recv_last_req);
            }
        }


        // all ranks except first: send first referenced row of uuu (i.e. first stored column) to previous rank
        for (int i = size-1; i > 0; --i) {
            // send
            if (rank == i) {
                MPI_Isend(uuu_across, Nx, MPI_DOUBLE, i-1, 0, MPI_COMM_WORLD, &u_send_first_req);
                if (mu_2 != 0) MPI_Isend(vvv_across, Nx, MPI_DOUBLE, i-1, 1, MPI_COMM_WORLD, &v_send_first_req);
            }

            // recv
            if (rank == i-1) {
                MPI_Irecv(uuu_top_bottom + Nx, Nx, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &u_recv_first_req);
                if (mu_2 != 0) MPI_Irecv(vvv_top_bottom + Nx, Nx, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &v_recv_first_req);
            }
        }


        count = 0;
        for (int j = int(cblas_dasum(rank, Ny_sizes, 1)); j < int(cblas_dasum(rank+1, Ny_sizes, 1)); ++j) {
            through_rows = count*Nx;

            for (int i = 0; i < Nx; ++i) {

                // deal with top boundary
                if (j == 0) {
                    u_lap_y[through_rows + i] = -uuu_across[i] + uuu_across[Nx + i];
                    if (mu_2 != 0) v_lap_y[through_rows + i] = -vvv_across[i] + vvv_across[Nx + i];
                }

                // deal with bottom boundary
                else if (j == Ny-1) {
                    u_lap_y[through_rows + i] = -uuu_across[through_rows + i] + uuu_across[(count-1)*Nx + i];
                    if (mu_2 != 0) v_lap_y[through_rows + i] = -vvv_across[through_rows + i] + vvv_across[(count-1)*Nx + i];
                }

                // // main inner part
                // not the top or bottom of a rank
                else if (count != 0 && count != Ny_size-1) {
                    u_lap_y[through_rows + i] = uuu_across[(count+1)*Nx + i] + uuu_across[(count-1)*Nx + i];
                    if (mu_2 != 0) v_lap_y[through_rows + i] = vvv_across[(count+1)*Nx + i] + vvv_across[(count-1)*Nx + i];
                }

                // top boundary of a rank
                else if (count == 0) {
                    MPI_Wait(&u_recv_last_req, MPI_STATUS_IGNORE);
                    u_lap_y[through_rows + i] = uuu_across[(count+1)*Nx + i] + uuu_top_bottom[i];

                    if (mu_2 != 0) {
                        MPI_Wait(&v_recv_last_req, MPI_STATUS_IGNORE);
                        v_lap_y[through_rows + i] = vvv_across[(count+1)*Nx + i] + vvv_top_bottom[i];
                    }
                }

                // bottom boundary of a rank
                else if (count == Ny_size-1) {
                    MPI_Wait(&u_recv_first_req, MPI_STATUS_IGNORE);
                    u_lap_y[through_rows + i] = uuu_top_bottom[Nx + i] + uuu_across[(count-1)*Nx + i];
                    
                    if (mu_2 != 0) {
                        MPI_Wait(&v_recv_first_req, MPI_STATUS_IGNORE);
                        v_lap_y[through_rows + i] = vvv_top_bottom[Nx + i] + vvv_across[(count-1)*Nx + i];
                    }
                }


            }
            ++count;
        }


        // // calculate f_1 and f_2
        count = 0;
        for (int j = int(cblas_dasum(rank, Ny_sizes, 1)); j < int(cblas_dasum(rank+1, Ny_sizes, 1)); ++j) {
            through_rows = count*Nx;

            // go across in i
            for (int i = 0; i < Nx; ++i) {
                index = through_rows + i;
                f_1[index] = eps * uuu_across[index] * (1-uuu_across[index]) * ( uuu_across[index] - ( vvv_across[index]*over_a + b_over_a ) );
                f_2[index] = uuu_across[index]*uuu_across[index]*uuu_across[index] - vvv_across[index];
            }

            ++count;
        }



        // // update uuu with the value of u^{n+1}
        // we are updating uuu_across
        count = 0;
        for (int j = int(cblas_dasum(rank, Ny_sizes, 1)); j < int(cblas_dasum(rank+1, Ny_sizes, 1)); ++j) {
            through_rows = count*Nx;

            // go across in i
            for (int i = 0; i < Nx; ++i) {
                index = through_rows + i;
                
                // deal with corners i.e. don't minus u_{i, j} at all
                if ((i==0 && j==0) || (i==Nx-1 && j==0) || (i==0 && j==Ny-1) || (i==Nx-1 && j==Ny-1)) {
                    // uuu[i*Ny + j] = dt * ( mu_1 * (u_lap_x[i*Ny + j] + u_lap_y[i*Ny + j] ) + f_1[i*Ny + j] ) + uuu[i*Ny + j];
                    uuu_across[index] = dt * ( mu_1 * (u_lap_x[index] + u_lap_y[index] ) + f_1[index] ) + uuu_across[index];
                }

                // the below block runs if we are not at a corner but we are at a boundary and so we need to: minus just 2 times u_{i, j}
                else if (i == 0 || i == Nx-1 || j == 0 || j == Ny-1) {
                    uuu_across[index] = dt * ( mu_1 * (u_lap_x[index] + u_lap_y[index] - 2 * uuu_across[index]) + f_1[index] ) + uuu_across[index];
                }
                else {
                    uuu_across[index] = dt * ( mu_1 * (u_lap_x[index] + u_lap_y[index] - 4 * uuu_across[index]) + f_1[index] ) + uuu_across[index];
                }


            }

            ++count;
        }


        // // update vvv with the value of v^{n+1}
        // we are updating vvv_across
        if (mu_2 != 0) { // use the full equation for v^{n+1}
            count = 0;
            for (int j = int(cblas_dasum(rank, Ny_sizes, 1)); j < int(cblas_dasum(rank+1, Ny_sizes, 1)); ++j) {
                through_rows = count*Nx;

                // go across in i
                for (int i = 0; i < Nx; ++i) {
                    index = through_rows + i;
                    
                    // deal with corners i.e. don't minus u_{i, j} at all
                    if ((i==0 && j==0) || (i==Nx-1 && j==0) || (i==0 && j==Ny-1) || (i==Nx-1 && j==Ny-1)) {
                        vvv_across[index] = dt * ( mu_2 * (v_lap_x[index] + v_lap_y[index] ) + f_2[index] ) + vvv_across[index];
                    }

                    // the below block runs if we are not at a corner but we are at a boundary and so we need to: minus just 2 times u_{i, j}
                    else if (i == 0 || i == Nx-1 || j == 0 || j == Ny-1) {
                        vvv_across[index] = dt * ( mu_2 * (v_lap_x[index] + v_lap_y[index] - 2 * vvv_across[index]) + f_2[index] ) + vvv_across[index];
                    }
                    else {
                        vvv_across[index] = dt * ( mu_2 * (v_lap_x[index] + v_lap_y[index] - 4 * vvv_across[index]) + f_2[index] ) + vvv_across[index];
                    }
                }

                ++count;
            }
        }
        else {
            count = 0;
            for (int j = int(cblas_dasum(rank, Ny_sizes, 1)); j < int(cblas_dasum(rank+1, Ny_sizes, 1)); ++j) {
                through_rows = count*Nx;

                // go across in i
                for (int i = 0; i < Nx; ++i) {
                    index = through_rows + i;

                    // when mu_2 = 0 the equation for v^{n+1} is the same at the boundaries as it is in the middle
                    vvv_across[index] = dt * ( f_2[index] ) + vvv_across[index];
                }

                ++count;
            }
        }


        ++main_count;
    }


    
    // make y_displacements array
    double y_displacements_double[size];
    int y_displacements[size];
    for (int i = 0; i < size; ++i) {
        y_displacements_double[i] = Nx*cblas_dasum(i, Ny_sizes, 1);
        y_displacements[i] = int(y_displacements_double[i]); // make y_displacements be of type int
    }

    // gather uuu_across arrays to uuu_main on rank 0
    MPI_Gatherv(uuu_across, Nx*Ny_size, MPI_DOUBLE, uuu_main, y_recvcounts, y_displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // gather vvv_across arrays to vvv_main on rank 0
    MPI_Gatherv(vvv_across, Nx*Ny_size, MPI_DOUBLE, vvv_main, y_recvcounts, y_displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    // write to file
    if (1 && rank == 0) {
        // cout << " \n Writing uuu and vvv to file" << endl;
        // write uuu and vvv to file
        ofstream myfile;
        myfile.open ("output.txt");

        for (int j = Ny-1; j >=0; --j) { // going from bottom to top i.e. from y=0 to y=height
            for (int i = 0; i < Nx; ++i) { // fastest changing dimension is x
                myfile << i << " " << (Ny-1-j) << " " << uuu_main[j*Nx + i] << " " << vvv_main[j*Nx + i]  << "\n";
                // cout << i << " " << (Ny-1-j) << " " << uuu_main[j*Nx + i] << " " << vvv_main[j*Nx + i]  << "\n";
            }
            myfile << "\n";
            // cout << "\n";
        }

        myfile.close();
        // cout << " Finished writing uuu and vvv to file" << endl;
    }
    }


}
