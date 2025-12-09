# Purpose: Integrates the DarkHistory and CLASS codes to calculate CMB anisotropies and thermal history
# influenced by Dark Matter (DM) decay, specifically designed for use as a Cobaya Theory class.

import sys
import os
import argparse
import yaml
from cobaya.theory import Theory
from copy import deepcopy
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
import numpy as np

from cobaya.conventions import Const

H_units_conv_factor = {
    "1/Mpc": 1,
    "km/s/Mpc": Const.c_km_s
}

# =========================
# DarkHistory Path Setup (MUST be done before importing darkhistory)
# Setup paths to locate DarkHistory data and source directory based on DH_DATA_DIR environment variable.
# =========================
base_data_dir = os.environ.get("DH_DATA_DIR")
if base_data_dir is None:
    raise EnvironmentError("Please set DH_DATA_DIR environment variable before running.")

# DarkHistory root directory is the parent of DH_DATA_DIR
base_dir = os.path.dirname(base_data_dir)
sys.path.append(base_dir)


# DarkHistory
import darkhistory.physics as phys
import darkhistory.main as main
from darkhistory.config import load_data

# CLASS
from classy import Class

def run_dh_class(cfg):
    # Parameter extraction and validation.
    required_keys = [
        "n_s",
        "A_s",
        "H0",
        "omega_cdm",
        "omega_b",
        "log10_mnu",
        "mDM",
        "lifetime_exponent",
        "alpha_gomp",
        "beta_gomp",
    ]

    for k in required_keys:
        if k not in cfg:
            raise ValueError(f"Missing required parameter: {k}")

    # -------------------------
    # Fix: compute h and omegas from cfg (previous code used cosmo_params before definition)
    # -------------------------
    h = float(cfg["H0"]) / 100.0
    omega_b = float(cfg["omega_b"])
    omega_cdm = float(cfg["omega_cdm"])
    Omega_b = omega_b / h**2
    Omega_cdm = omega_cdm / h**2

    cosmo_params = {
        "H0": float(cfg["H0"]),
        "Omega_cdm": Omega_cdm,
        "Omega_b": Omega_b,
        "mnu": 10 ** float(cfg["log10_mnu"]),
    }

    mDM = float(cfg["mDM"])
    lifetime_exponent = float(cfg["lifetime_exponent"])
    lifetime = 10 ** lifetime_exponent
    lna_pivot = float(cfg["alpha_gomp"])
    tilt = float(cfg["beta_gomp"])
    inj_particle = cfg.get("inj_particle", "decay_e")

    # Set cosmology params globally for DarkHistory.
    phys.set_cosmology_params(cosmo_params)

    # Construct unique filename for caching results.
    M = mDM / 1e9
    tau_str = "{:.6f}".format(lifetime_exponent)

    decay_dir = os.path.join(base_dir, inj_particle)
    os.makedirs(decay_dir, exist_ok=True)

    filename = os.path.join(
        decay_dir,
        f"{M:.6f}_{tau_str}_{cosmo_params['H0']:.6f}_{cosmo_params['Omega_cdm']:.6f}_{cosmo_params['Omega_b']:.6f}_{cosmo_params['mnu']:.6f}_{tilt:.6f}_{lna_pivot:.6f}.txt"
    )

    # =========================
    # DarkHistory Execution
    # Runs the DarkHistory evolution if the result file does not exist, calculating energy deposition.
    # =========================
    binning = load_data('binning')
    dep_tf_data = load_data('dep_tf')
    ics_tf_data = load_data('ics_tf')
    if not os.path.exists(filename):
        try:
            _ = main.evolve(
                mDM=mDM, DM_process='decay',
                lifetime=lifetime,
                primary=inj_particle,
                start_rs=3000, end_rs=4,
                coarsen_factor=30,
                backreaction=True,
                reion_switch=True,
                reion_rs=800,
                lna_pivot=lna_pivot, tilt=tilt,
                binning=binning,
                dep_tf_data=dep_tf_data,
                ics_tf_data=ics_tf_data
            )
        except Exception as e:
            # record error to file for debugging
            with open("/tmp/dh_errors.log", "a") as f:
                f.write(f"DH ERROR for params {cfg}: {e}\n")
            raise RuntimeError(f"DarkHistory failed: {e}")

    # =========================
    # CLASS Execution
    # Runs the CLASS code using parameters consistent with DarkHistory to compute CMB spectra.
    # =========================

    # Configuration for CLASS, including Gomp-like reionization and DM decay source.
    Gomp_DM = {
        'output': 'tCl,pCl,lCl,mPk',
        'lensing': 'yes',
        
        'n_s': float(cfg["n_s"]),
        'A_s': float(cfg["A_s"]),

        'h': h,
        'omega_b': omega_b,
        'omega_cdm': omega_cdm,

        # massive neutrino
        'N_ur': 2.0308,
        'T_ncdm': 0.71611,
        'N_ncdm': 1,
        'm_ncdm': cosmo_params["mnu"],

        # gomp reionization
        'reio_parametrization': 'reio_gomp_noSR',
        'z_reio': 800,
        'beta_gomp': tilt,
        'alpha_gomp': lna_pivot,

        # DM decay
        'DH_DM_decay_flag': 1,
        'DH_DM_decay_mass': M,
        'DH_DM_decay_lifetime_exponent': lifetime_exponent,

        # high accuracy settings described in Figure 48 of 2503.14454
        'T_cmb': 2.7255,
        'YHe': 'BBN',
        'non_linear': 'hmcode',
        'hmcode_version': '2020',
        'recombination': 'HyRec',
        'modes': 's',
        'l_max_scalars': 9500,
        'delta_l_max': 1800,
        'P_k_max_h/Mpc': 100.,
        'l_logstep': 1.025,
        'l_linstep': 20,
        'perturbations_sampling_stepsize': 0.05,
        'l_switch_limber': 30.,
        'hyper_sampling_flat': 32.,
        'l_max_g': 40,
        'l_max_ur': 35,
        'l_max_pol_g': 60,
        'ur_fluid_approximation': 2,
        'ur_fluid_trigger_tau_over_tau_k': 130.,
        'radiation_streaming_approximation': 2,
        'radiation_streaming_trigger_tau_over_tau_k': 240.,
        'hyper_flat_approximation_nu': 7000.,
        'transfer_neglect_delta_k_S_t0': 0.17,
        'transfer_neglect_delta_k_S_t1': 0.05,
        'transfer_neglect_delta_k_S_t2': 0.17,
        'transfer_neglect_delta_k_S_e': 0.17,
        'accurate_lensing': 1,
        'start_small_k_at_tau_c_over_tau_h': 0.0004,
        'start_large_k_at_tau_h_over_tau_k': 0.05,
        'tight_coupling_trigger_tau_c_over_tau_h': 0.005,
        'tight_coupling_trigger_tau_c_over_tau_k': 0.008,
        'start_sources_at_tau_c_over_tau_h': 0.006,
        'l_max_ncdm': 30,
        'tol_ncdm_synchronous': 1.e-6
    }

    try:
        G = Class()
        G.set(Gomp_DM)
        G.compute()
    except Exception as e:
        with open("/tmp/class_errors.log", "a") as f:
            f.write(f"CLASS ERROR for params {cfg}: {e}\n")
        raise RuntimeError(f"CLASS failed: {e}")

    derived = G.get_current_derived_parameters(['z_reio', 'tau_reio', 'conformal_age', 'sigma8',
                                               'theta_s_100', 'Omega_m', 'YHe', 'age'])

    T_cmb = G.T_cmb()
    rs_drag = G.rs_drag()

    return {
        'T_cmb': T_cmb,
        'z_reio': derived['z_reio'],
        'tau_reio': derived['tau_reio'],
        'conformal_age': derived['conformal_age'],
        'sigma8': derived['sigma8'],
        'theta_s_100': derived['theta_s_100'],
        'Omega_m': derived['Omega_m'],
        'YHe': derived['YHe'],
        'age': derived['age'],
        'rs_drag': rs_drag,
        'G': G
    }


# =======================================
# Cobaya Theory Class
# Provides an interface for Cobaya to run the DH+CLASS calculation and retrieve results (Cls and derived parameters).
# =======================================
class DHClassTheory(Theory):

    # Define parameters that Cobaya will sample.
    params = {
        "n_s": None,
        "A_s": None,
        "H0": None,
        "omega_cdm": None,
        "omega_b": None,

        "log10_mnu": None,

        "mDM": None,
        "lifetime_exponent": None,

        "alpha_gomp": None,
        "beta_gomp": None
    }
    
    def get_can_provide_params(self):

        return [
            "z_reio",
            "tau_reio",
            "conformal_age",
            "sigma8",
            "theta_s_100",
            "Omega_m",
            "YHe",
            "age",
            "rs_drag",
            "rdrag"
        ]



    def calculate(self, state, want_derived=True, **params_values_dict):
        derived = run_dh_class(params_values_dict)
        cl_dict = derived['G'].lensed_cl()

        state['Cl'] = {
            'tt': cl_dict['tt'],
            'ee': cl_dict['ee'],
            'te': cl_dict['te'],
            'bb': cl_dict['bb'],
            'pp': cl_dict['pp'],
            'tp': cl_dict['tp'],
            'ell': cl_dict['ell']
        }

        # Save T_cmb to current state for get_Cl usage
        state["T_cmb"] = derived.get("T_cmb", 2.7255)

        state["derived"] = {
            "z_reio": derived.get("z_reio"),
            "tau_reio": derived.get("tau_reio"),
            "conformal_age": derived.get("conformal_age"),
            "sigma8": derived.get("sigma8"),
            "theta_s_100": derived.get("theta_s_100"),
            "Omega_m": derived.get("Omega_m"),
            "YHe": derived.get("YHe"),
            "age": derived.get("age"),
            "rs_drag": derived.get("rs_drag"),
            "rdrag": derived.get("rs_drag")
        }

        state["_class_instance"] = derived.get("G")


        return state
    

    def get_result(self, result_name, **kwargs):

        if result_name == "rdrag":
            return np.array([self.current_state.get("rs_drag")])

        class_instance = self.current_state.get("_class_instance")
        if class_instance is None:
            raise RuntimeError("CLASS instance not found")

        z = kwargs.get("z")
        if z is None:
            raise ValueError(f"{result_name} requires z")

        # Always treat z as list
        z_arr = np.atleast_1d(z)

        if result_name == "angular_diameter_distance":
            return np.array([class_instance.angular_distance(zi) for zi in z_arr])

        if result_name == "Hubble":
            return np.array([class_instance.Hubble(zi) for zi in z_arr])

        raise ValueError(f"Result '{result_name}' not provided by DHClassTheory.")


    def get_angular_diameter_distance(self, z):
        class_instance = self.current_state.get("_class_instance")
        if class_instance is None:
            raise RuntimeError("CLASS instance not found.")

        if np.isscalar(z):
            return np.array([class_instance.angular_distance(z)])

        return np.array([class_instance.angular_distance(zi) for zi in z])



    def get_Hubble(self, z, units="km/s/Mpc"):
        class_instance = self.current_state.get("_class_instance")
        if class_instance is None:
            raise RuntimeError("CLASS instance not found.")

        if np.isscalar(z):
            hz = class_instance.Hubble(z)
            hz *= H_units_conv_factor[units]
            return np.array([hz]) 

        hz = np.array([class_instance.Hubble(zi) for zi in z])
        hz *= H_units_conv_factor[units]
        return hz



    def get_can_provide(self):
        return[
            "angular_diameter_distance", 
            "Hubble",
            "rdrag"
        ]

    def get_Cl(self, ell_factor=False, units="FIRASmuK2"):
        try:
            cls = deepcopy(self.current_state["Cl"])
        except Exception:
            raise RuntimeError(f"No {'Cl'} found in current_state. Have you computed them?")

        # ℓ array
        ell = cls["ell"]
        ell_factor_array = ((ell * (ell + 1)) / (2 * np.pi))[2:] if ell_factor else 1.0

        T_cmb = self.current_state.get("T_cmb", 2.7255)  # Default to 2.7255 K if not set

        if units in ["1", 1]:
            units_factor = 1.0
        elif units in ["muK2"]:
            units_factor = T_cmb * 1.0e6
        elif units in ["K2"]:
            units_factor = T_cmb
        elif units in ["FIRASmuK2"]:
            units_factor = 2.7255e6
        elif units in ["FIRASK2"]:
            units_factor = 2.7255
        else:
            raise RuntimeError(f"Units '{units}' not recognized.")

        for cl_name in cls:
            if cl_name == "ell":
                continue
            # Determine power: count 't','e','b' in the name
            units_power = float(sum(cl_name.count(p) for p in ["t", "e", "b"]))
            # Apply units_factor
            cls[cl_name][2:] *= units_factor ** units_power

            # Apply ell_factor
            if ell_factor:
                if "p" not in cl_name:
                    cls[cl_name][2:] *= ell_factor_array
                elif cl_name == "pp":
                    cls[cl_name][2:] *= ell_factor_array ** 2 * (2 * np.pi)
                elif "p" in cl_name:
                    cls[cl_name][2:] *= ell_factor_array ** (3 / 2) * np.sqrt(2 * np.pi)

        return cls
