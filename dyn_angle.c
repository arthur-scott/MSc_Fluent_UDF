#include "udf.h"

/* Dynamic Contact angle UDF */

/* Initilisation of parameters */

// Surface properties

typedef enum{HYDROPHILIC, HYDROPHOBIC, N_HYDRO_TYPES} HydroType;

// Creating arrays for each parameter for the different wall types

static double theta_ed_ar[N_HYDRO_TYPES] = {100.0, 100.0}; // Static contact angle
static double mu_ar[N_HYDRO_TYPES] = {0.001003, 0.001003}; // Viscosity
static double sigma_ar[N_HYDRO_TYPES] = {0.0728, 0.0728}; // Surface tension
static double theta_a_ar[N_HYDRO_TYPES] = {105.0, 105.0}; // Advancing contact angle
static double theta_r_ar[N_HYDRO_TYPES] = {95.0, 95.0}; //Receding contact angle

HydroType get_hydro_type(Thread *ft)
{
    int i;
    int n_hydrophilic = 1;
    int hydrophilic_threads[] = {23};
    int thread_id;

    thread_id = THREAD_ID(ft);

    for (i=0; i<n_hydrophilic; i++)
        if (thread_id == hydrophilic_threads[i])
            return HYDROPHILIC;
    
    return HYDROPHOBIC;
}

// UDM Specification

enum{
    VOF_G_X,
    VOF_G_Y,
    VOF_G_Z,
    VTFS,
    CONTACT_ANGLE,
    REL_VEL_U,
    REL_VEL_V,
    REL_VEL_W,
    C_ucell ,
    C_lambda ,
    C_nw,
    C_tw ,
    CONTACT_N_REQ_UDM
};

#define C_VOF_G_X(C, T)C_UDMI(C, T,VOF_G_X)
#define C_VOF_G_Y(C, T)C_UDMI(C, T,VOF_G_Y)
#define C_VOF_G_Z(C, T)C_UDMI(C, T,VOF_G_Z)
#define C_VTFS(C, T)C_UDMI(C, T, VTFS)
#define C_CONTACT_ANGLE(C, T)C_UDMI(C, T,CONTACT_ANGLE)
#define C_REL_VEL_U(C, T)C_UDMI(C, T, REL_VEL_U)
#define C_REL_VEL_V(C, T)C_UDMI(C, T, REL_VEL_V)
#define C_REL_VEL_W(C, T)C_UDMI(C, T, REL_VEL_W)
#define C_ucell (C, T)C_UDMI(C, T, C_ucell )
#define C_lambda (C, T)C_UDMI(C, T, C_lambda )
#define C_tw(C, T)C_UDMI(C, T, C_tw)
#define C_nw(C, T)C_UDMI(C, T, C_nw)

// Main DEFINE_ADJUST loop

DEFINE_ADJUST(dyn_angle, domain)
{
    Thread *ct, *pct, *t0; // current thread, parent thread, master thread
    Thread *ft; // face thread
    cell_t c, c0; // current cell, an adjacent cell - used to compute gradients or source terms
    face_t f; // current face

    double theta_n;
    real Vtfs, temp[ND_ND], vel_re[ND_ND], vof_n[ND_ND], vec_nw[ND_ND], vec_tw[ND_ND], vec_snw[ND_ND]; // simply declaring variables with ND_ND dimensions
    real vec_ntw[ND_ND], vec_stw[ND_ND], A[ND_ND], Amag, Smag, u_cell[ND_ND], veltan[ND_ND];
    double f_Hoff_inverse, x_hoff, Ca, gapp, xapp, lambda, theta_ed, mu, sigma, theta_a, theta_r;

    int phase_domain_index = 1;
    Domain *pDomain = DOMAIN_SUB_DOMAIN(domain, phase_domain_index); /* selecting the index 1 secondary (water) phase, creates a pointer called pDomain which
    points to water by DOMAIN_SUB_DOMAIN. takes main domain and extracts specific phase. Hence, pDomain == water phase */

    HydroType hydro_type;

    if(first_iteration) // calls only if first iteration, inbuilt ANSYS terminology
    {
        // First checking enough UDM allocation

        if(N_UDM < CONTACT_N_REQ_UDM)
        {
            Message0("\n WARNING: Require at least %d UDMs to be setup. \n", CONTACT_N_REQ_UDM);
            return;
        }

        // Calculating VOF gradient

        Alloc_Storage_Vars(pDomain, SV_VOF_RG, SV_VOF_G, SV_NULL); // primary storage of variables
        Scalar_Reconstruction(pDomain, SV_VOF, -1, SV_VOF_RG, NULL);
        Scalar_Derivatives(pDomain, SV_VOF, -1, SV_VOF_G, SV_VOF_RG, Vof_Deriv_Accumulate);

        // Storing VOF gradient in memory

        thread_loop_c(ct, domain) // loops through all threads
            if (FLUID_THREAD_P(ct)) // checks if current thread is a fluid
            {
                pct = THREAD_SUB_THREAD(ct, phase_domain_index); // accesses water phase thread

                begin_c_loop(c, ct) // loops over all cells in the thread
                {
                    ND_V(C_VOF_G_X(c, ct), C_VOF_G_Y(c, ct), C_VOF_G_Z(c, ct), =, C_VOF_G(c, pct)); // calculates gradient in all cells
                }
                end_c_loop(c, ct)
            }
        
        // Free memory for VOF gradient

        Free_Storage_Vars(pDomain, SV_VOF_RG, SV_VOF_G, SV_NULL);

        // Calculation of wall contact angle

        thread_loop_f(ft, domain) // looping over all threads
        {
            if (THREAD_TYPE(ft) == THREAD_F_WALL) // checking if thread is a wall thread
            {
                t0 = THREAD_T0(ft);

                if (FLUID_THREAD_P(t0)) // fluid zone at wall
                {
                    // Setting contact angle parameters for the particular type of wall surface (hydrophilic or hydrophobic)

                    HydroType hydro_type; // initiating hydro_type enum

                    hydro_type = get_hydro_type(ft); // getting hydro type of face thread current

                    theta_ed = theta_ed_ar[hydro_type]; // setting variables for contact angle to the wall type's parameters
                    mu = mu_ar[hydro_type];
                    sigma = sigma_ar[hydro_type];
                    theta_a = theta_a_ar[hydro_type];
                    theta_r = theta_r_ar[hydro_type];

                    // Looping through all faces in thread

                    begin_f_loop(f, ft)
                    {
                        c0 = F_C0(f, ft); // accessing properties of current cell (face) in loop

                        // Flow velocity relative to face
                        NV_DD(vel_re ,=, C_U(c0, t0), C_V(c0, t0), C_W(c0, t0), -, F_U(f, ft), F_V(f, ft), F_W(f, ft));
                        ND_V(C_REL_VEL_U(c0, t0), C_REL_VEL_V(c0, t0), V_REL_VEL_W(c0, t0) ,=, vel_re);

                        F_AREA(A, f, ft); // area vector
                        Amag = NV_MAG(A);
                        NV_S(A, /=, Amag); // A holds outwards facing normal

                        // Navier slip
                        lambda = (C_VOLUME(c0, t0)/Amag)/2;
                        NV_D(u_cell ,=, lambda*C_DUDX(c0, t0), lambda*C_DVDY(c0, t0), lambda*C_DWDZ(c0, t0));

                        C_lambda(c0, t0) = lambda;
                        C_ucell(c0, t0) = NV_MAG(u_cell);

                        NV_D(vof_n ,=, C_VOF_G_X(c0, t0), C_VOF_G_Y(c0, t0), C_VOF_G_Z(c0, t0));

                        // VoF along the wall
                        
                        NV_VS(vec_nw, =, A, *, NV_DOT(vel_re, A)); // normal relative velocity vector
                        NV_VV(vec_tw, =, vel_re, -, vec_nw); // tangential velocity vector

                        C_nw(c0, t0) = NV_MAG(vec_nw);
                        C_tw(c0, t0) = NV_MAG(vec_tw);

                        NV_VS(vec_snw, =, A, *, NV_DOT(vof_n, A)); // VoF gradient vector normal
                        NV_VV(vec_stw, =, vof_n, -, vec_snw); // VoF gradient vector tangential

                        // Contact line direction

                        NV_CROSS(temp, u_cell, vec_nw); // takes cross product of u_cell and vec_nw, stores in temp
                        NV_CROSS(veltan, vec_nw, temp); // takes cross product of vec_nw and temp, stores in veltan

                        // Receding contact angle

                        if (NV_DOT(veltan, vof_n) > 0.0) // if dot product of tangential veloctity and ...
                        {
                            Ca = mu*NV_MAG(veltan)/sigma; // capillary number
                            theta_n =  DEGREES(cbrt(CUB(RADIANS(theta_r)) - 72.0*Ca));

                            C_CONTACT_ANGLE(c0, t0) = theta_n;
                        }
                        else
                        {
                            // Advancing

                            Ca = mu*NV_MAG(veltan)/sigma; // capillary number

                            // Kistler's law

                            f_Hoff_inverse = CUB(RADIANS(theta_ed))/72.0;
                            x_hoff = Ca + f_Hoff_inverse;
                            theta_n = DEGREES(acos(1.0 - 2.0*tanh(5.16*pow((x_hoff/(1.0 + 1.31*pow(x_hoff, 0.99))), 0.706))));

                            C_CONTACT_ANGLE(c0, t0) = theta_n;
                        }
                        C_VTFS(c0, t0) = NV_DOT(veltan, vof_n);
                    }
                    end_f_loop(f, t_m)
                }
            }
        }
    }
}

DEFINE_PROFILE(contact_angle, t, i)
{
    face_t f;
    cell_t c0;
    Thread *t0;

    begin_f_loop(f, t)
    {
        c0 = F_C0(f, t);
        t0 = F_C0_THREAD(f, t);

        F_PROFILE(f, t, i) = RADIANS(C_CONTACT_ANGLE(c0, t0));
    }
    end_f_loop(f, t)
}

DEFINE_EXECUTE_AFTER_CASE(set_name, libname)
{
    Set_User_Memory_Name(VOF_G_X, "VOF Gradient - x");
    Set_User_Memory_Name(VOF_G_Y, "VOF Gradient - y");
    Set_User_Memory_Name(VOF_G_Z, "VOF Gradient - z");
    Set_User_Memory_Name(VTFS, "Vtfs");
    Set_User_Memory_Name(CONTACT_ANGLE, "VOF Gradient - x");
    Set_User_Memory_Name(REL_VEL_U, "VOF Gradient - x");
    Set_User_Memory_Name(REL_VEL_V, "VOF Gradient - x");
    Set_User_Memory_Name(REL_VEL_W, "VOF Gradient - x");
    Set_User_Memory_Name(C_ucell, "Magnitude of cell velocity");
    Set_User_Memory_Name(C_lambda, "lambda");
    Set_User_Memory_Name(C_nw, "Magnitude of normal vector");
    Set_User_Memory_Name(C_tw, "Magnitude of tangential vector");
    Message0("Done. \n");
}

DEFINE_ON_DEMAND(check_settings)
{
    Message0("\n");
    Message0("Model parameter values: \n");
    Message0("\n");
    Message0("Hydrophilic: \n");
    Message0("theta_ed = %e \n", theta_ed_ar[HYDROPHILIC]);
    Message0("mu = %e \n", mu_ar[HYDROPHILIC]);
    Message0("sigma = %e \n", sigma_ar[HYDROPHILIC]);
    Message0("theta_a = %e \n", theta_a_ar[HYDROPHILIC]);
    Message0("theta_r = %e \n", theta_r_ar[HYDROPHILIC]);
    Message0("\n");
    Message0("Hydrophobic: \n");
    Message0("theta_ed = %e \n", theta_ed_ar[HYDROPHOBIC]);
    Message0("mu = %e \n", mu_ar[HYDROPHOBIC]);
    Message0("sigma = %e \n", sigma_ar[HYDROPHOBIC]);
    Message0("theta_a = %e \n", theta_a_ar[HYDROPHOBIC]);
    Message0("theta_r = %e \n", theta_r_ar[HYDROPHOBIC]);
    Message0("\n");
}