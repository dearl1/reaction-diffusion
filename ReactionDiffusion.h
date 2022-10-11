
#ifndef CLASS_REACTIONDIFFUSION
#define CLASS_REACTIONDIFFUSION

class ReactionDiffusion {
    public:
        // user-defined constructor
        ReactionDiffusion(double dtf, double integration_timef, int Nxf, 
            int Nyf, double af, double bf, double mu_1f, 
            double mu_2f, double epsf, int rankf, int sizef);

        void SetInitialConditions();

        void TimeIntegrate();
    
    private:
        double dt;
        double integration_time;
        int Nx;
        int Ny;
        double a;
        double b;
        double mu_1;
        double mu_2;
        double eps;

        int rank;
        int size;

        int Ny_size;
        double* Ny_sizes;
        int* y_recvcounts;

        double* uuu_across;
        double* vvv_across;
        double* uuu_top_bottom;
        double* vvv_top_bottom;

        double* uuu;
        double* vvv;

        double* uuu_main;
        double* vvv_main;

        double* uuu_work;
        double* vvv_work;
        double* uuu_work_send_horiz;
        double* vvv_work_send_horiz;

        double* u_lap_x; // for updating uuu
        double* v_lap_x; // for updating vvv
        double* u_lap_y; // for updating uuu
        double* v_lap_y; // for updating vvv
        double* f_1;
        double* f_2;
    
        MPI_Request u_send_up, u_send_down, u_send_left, u_send_right, v_send_up, v_send_down, v_send_left, v_send_right;
        MPI_Request u_recv_up, u_recv_down, u_recv_left, u_recv_right, v_recv_up, v_recv_down, v_recv_left, v_recv_right;

        MPI_Request u_send_last_req, v_send_last_req, u_recv_last_req, v_recv_last_req;
        MPI_Request u_send_first_req, v_send_first_req, u_recv_first_req, v_recv_first_req;

        double b_over_a;
        double over_a;

        // make the cartesian communicator
        int num_1D;
        MPI_Comm mygrid;
        int dims;
        int* sizes;
        int* periods;
        int reorder;

        // get coords for each process
        int mygrid_rank;
        int* coords;

        // make row and column communicators
        MPI_Comm myrowcomm, mycolcomm;
        int* keep;

        int row_rank, row_size;
        int col_rank, col_size;

        // split up domain in y-direction
        int remainder_y;
        int Ny_small;
        int Ny_big;

        // make Ny_sizes_all
        double* Ny_sizes_all;


        // split up domain in x-direction
        int remainder_x;
        int Nx_small;
        int Nx_big;

        double* Nx_sizes;

        
        int Nx_size;


        // make Nx_sizes_all
        double* Nx_sizes_all;


        // initialise the uuu and vvv arrays
        int count_j;
        int count_i;




   
};

#endif