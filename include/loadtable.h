#ifndef LOADTABLE_H
#define LOADTABLE_H

struct DarkHistory_decay{
    int DM_decay_flag;

    /** DM decay parameters */
    double DM_mass;                  // [GeV]
    double DM_lifetime;              // [s]

    /** Cosmological parameters */
    double Cosmo_hubble_constant;    // [km/s/Mpc]
    double Cosmo_Omega_DM;          
    double Cosmo_Omega_baryon;
    double Cosmo_mnu;                // [eV]

    /** Gomp reionization parameters */
    double Gomp_tilt;
    double Gomp_lna_pivot;
};


struct DarkHistory_loaded_info {
    double DH_loaded_mass;
    double DH_loaded_lifetime;
    double DH_loaded_H0;
    double DH_loaded_Omega_cdm;
    double DH_loaded_Omega_b;
    double DH_loaded_mnu;
    double DH_loaded_tilt;
    double DH_loaded_lna_pivot;
};

#ifdef __cplusplus
extern "C" {
#endif

    double Calc_dxedz_dTdz(int T_or_xe, double z, struct DarkHistory_decay *pDH);
    double two_pts_interpolation(int use_log, double x0_, double x1_, double y0_, double y1_, double xin_);

#ifdef __cplusplus 
}
#endif


#endif /* LOADTABLE_H */