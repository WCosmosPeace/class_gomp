#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "loadtable.h"
#include <errno.h>

#define POINTS 220  //DarkHistory redshift points.

static double DH_z[POINTS], DH_data[POINTS][3]; //DH_data[redshift][0/1/2](0: dxHII/dz, 1: dxHeII/dz, 2: dT/dz).
static int DH_loaded = 0;   //The default value is 0; (1: decay).
static int interp_method = 1; //(0: linear and 1: log)
static struct DarkHistory_loaded_info last_loaded_params = {0};

double Calc_dxedz_dTdz(int T_or_xe, double z, struct DarkHistory_decay *pDH){

    if (pDH->DM_decay_flag == 1){

        double exp_DM_lifetime = log10(pDH->DM_lifetime);

        char filename[1024];
        FILE *infile;

        const char* data_path = getenv("DH_DATA_DIR");

        int i;
        double z0, z1;

        double y0, y1, result;

        if (data_path == NULL) {
            fprintf(stderr, "Error: Environment variable `DH_DATA_DIR` is not set.\n");
            fprintf(stderr, "Please set this variable, for example: `export DH_DATA_DIR=/DarkHistory/data\n");
            exit(1); 
        }

        int needs_reload = 0;

        if (!DH_loaded) {
            needs_reload = 1;
        } else {
            if (pDH->DM_mass != last_loaded_params.DH_loaded_mass ||
                pDH->DM_lifetime != last_loaded_params.DH_loaded_lifetime ||
                pDH->Cosmo_hubble_constant != last_loaded_params.DH_loaded_H0 ||
                pDH->Cosmo_Omega_DM != last_loaded_params.DH_loaded_Omega_cdm ||
                pDH->Cosmo_Omega_baryon != last_loaded_params.DH_loaded_Omega_b ||
                pDH->Cosmo_mnu != last_loaded_params.DH_loaded_mnu ||
                pDH->Gomp_tilt != last_loaded_params.DH_loaded_tilt ||
                pDH->Gomp_lna_pivot != last_loaded_params.DH_loaded_lna_pivot) 
            {
                needs_reload = 1;
            }
        }

        if (needs_reload) {
            snprintf(filename, sizeof(filename),
                    "%s/../decay_e/%f_%f_%f_%f_%f_%f_%f_%f.txt",
                    data_path,
                    pDH->DM_mass,
                    exp_DM_lifetime,
                    pDH->Cosmo_hubble_constant,
                    pDH->Cosmo_Omega_DM,
                    pDH->Cosmo_Omega_baryon,
                    pDH->Cosmo_mnu,
                    pDH->Gomp_tilt,
                    pDH->Gomp_lna_pivot
            );
            printf("Attempting to open file: %s\n", filename);

            infile = fopen(filename, "r");
            if (infile == NULL) {
                fprintf(stderr, "Error opening file '%s': %s\n", filename, strerror(errno));
                return -1;
            }

            char line[1024];
            while (fgets(line, sizeof(line), infile)) {
                if (line[0] != '#') {
                    break;
                }
            }

            if (sscanf(line, "%lf %lf %lf %lf", &DH_data[0][0], &DH_data[0][1], &DH_data[0][2], &DH_z[0]) != 4) {
                fprintf(stderr, "Error reading first data line from file.\n");
                fclose(infile);
                return -1;
            }

            for (i = 1; i < POINTS; i++) {
                if (fscanf(infile, "%lf %lf %lf %lf", &DH_data[i][0], &DH_data[i][1], &DH_data[i][2], &DH_z[i]) != 4) {
                    fprintf(stderr, "Error reading data from file.\n");
                    fclose(infile);
                    return -1;
                }
            }
            
            fclose(infile);
            DH_loaded = 1; 

            last_loaded_params.DH_loaded_mass = pDH->DM_mass;
            last_loaded_params.DH_loaded_lifetime = pDH->DM_lifetime;
            last_loaded_params.DH_loaded_H0 = pDH->Cosmo_hubble_constant;
            last_loaded_params.DH_loaded_Omega_cdm = pDH->Cosmo_Omega_DM;
            last_loaded_params.DH_loaded_Omega_b = pDH->Cosmo_Omega_baryon;
            last_loaded_params.DH_loaded_mnu = pDH->Cosmo_mnu;
            last_loaded_params.DH_loaded_tilt = pDH->Gomp_tilt;
            last_loaded_params.DH_loaded_lna_pivot = pDH->Gomp_lna_pivot;            
        }

        // Find the correct index
        int i0 = -1, i1 = -1;
        for (i=0; i < POINTS - 1; i++) {
            if (DH_z[i] >= z && DH_z[i+1] <= z) { // DH_z is descending
                i0 = i;
                i1 = i + 1;
                break;
            }
        }

        // Check if z is within the valid range
        if (i0 == -1 || i1 == -1) {
            return 0;
        }

        z0 = DH_z[i0];
        z1 = DH_z[i1];

        if (T_or_xe == 0) {           // dxHII/dz
            y0 = DH_data[i0][0];
            y1 = DH_data[i1][0];
            result = two_pts_interpolation(interp_method, z0, z1, y0, y1, z);
        } else if (T_or_xe == 1) {    // dxHeII/dz
            y0 = DH_data[i0][1];
            y1 = DH_data[i1][1];
            result = two_pts_interpolation(interp_method, z0, z1, y0, y1, z);
        } else if (T_or_xe == 2) {    // dT/dz
            y0 = DH_data[i0][2];
            y1 = DH_data[i1][2];
            result = two_pts_interpolation(interp_method, z0, z1, y0, y1, z);
        }

        return result;
    }
    else {
        fprintf(stderr, "Error: The `decay_flag` value can only be 0 or 1. Please check your input.\n");
        return -1; 
    }
}


double two_pts_interpolation(int use_log, double x0_, double x1_, double y0_, double y1_, double xin_) {
    double x0 = x0_, x1 = x1_, xin = xin_;
    double y0 = y0_, y1 = y1_;

    if (x0 == x1) {
        return y0;
    }

    if (use_log == 1) {
        x0 = log10(x0);
        x1 = log10(x1);
        xin = log10(xin);

        if (y0 >= 0 || y1 >= 0) {
             fprintf(stderr, "Error: Cannot perform log interpolation on non-negative values. Because `dxHII/dz` and `dT/dz` are negative values.\n");
             return 0;
        }
        
        y0 = log10(-y0);
        y1 = log10(-y1);
    }
    double k = (y1 - y0) / (x1 - x0);
    double ans = y0 + k * (xin - x0);

    if (use_log == 1) {
        ans = -pow(10.0, ans);
    }

    return ans;
}
