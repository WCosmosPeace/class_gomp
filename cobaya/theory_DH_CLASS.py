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


# =========================
# Step 1: Input parameters
# Defines command-line arguments and loads configuration from YAML, allowing CLI arguments to override YAML settings.
# =========================
def parse_args():
    parser = argparse.ArgumentParser(description="Run DarkHistory + CLASS with custom parameters")

    # YAML
    parser.add_argument("--yaml", type=str, help="Path to YAML configuration file")

    # Cosmology
    parser.add_argument("--H0", type=float)
    parser.add_argument("--Omega_cdm", type=float)
    parser.add_argument("--Omega_b", type=float)
    parser.add_argument("--mnu", type=float)

    # Dark matter decay
    parser.add_argument("--mDM", type=float)
    parser.add_argument("--lifetime_exponent", type=float, help="log10(DM lifetime [s])")
    parser.add_argument("--inj_particle", type=str, default="decay_e")

    # Reionization parameters
    parser.add_argument("--lna_pivot", type=float)
    parser.add_argument("--tilt", type=float)

    return parser.parse_args()


def load_config(args):
    """Merges YAML configuration and command-line arguments; CLI overwrites YAML."""
    config = {}

    if args.yaml:
        with open(args.yaml, "r") as f:
            config = yaml.safe_load(f)

    for key, value in vars(args).items():
        if value is not None and key != "yaml":
            config[key] = value

    return config


# =========================
# Step 2: Main physics computation function
# Executes DarkHistory and CLASS calculations for a given set of cosmological and DM parameters.
# =========================
def run_dh_class(cfg):
    # Parameter extraction and validation.
    required_keys = ["H0", "Omega_cdm", "Omega_b", "log10_mnu", "mDM", "lifetime_exponent", "lna_pivot", "tilt"]
    for k in required_keys:
        if k not in cfg:
            raise ValueError(f"Missing required parameter: {k}")

    # Prepare cosmology parameters for DarkHistory.
    cosmo_params = {
        "H0": float(cfg["H0"]),
        "Omega_cdm": float(cfg["Omega_cdm"]),
        "Omega_b": float(cfg["Omega_b"]),
        "mnu": 10**float(cfg["log10_mnu"]),
    }

    mDM       = float(cfg["mDM"])
    lifetime_exponent = float(cfg["lifetime_exponent"])
    lifetime  = 10**lifetime_exponent    
    lna_pivot = float(cfg["lna_pivot"])
    tilt      = float(cfg["tilt"])
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
                lna_pivot=lna_pivot, tilt=tilt
            )
        except Exception as e:
            # record error to file for debugging
            with open("/tmp/dh_errors.log","a") as f:
                f.write(f"DH ERROR for params {cfg}: {e}\n")
            raise RuntimeError(f"DarkHistory failed: {e}")


    # =========================
    # CLASS Execution
    # Runs the CLASS code using parameters consistent with DarkHistory to compute CMB spectra.
    # =========================
    h = cosmo_params["H0"]/100
    omega_cdm = cosmo_params["Omega_cdm"] * h**2
    omega_b   = cosmo_params["Omega_b"] * h**2

    # Configuration for CLASS, including Gomp-like reionization and DM decay source.
    Gomp_DM = {
        'output' : 'tCl,pCl,lCl,mPk',
        'lensing': 'yes',
        'l_max_scalars': 3000, 

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
        'DH_DM_decay_lifetime_exponent': lifetime_exponent
    }

    try:
        G = Class()
        G.set(Gomp_DM)
        G.compute()
    except Exception as e:
        with open("/tmp/class_errors.log","a") as f:
            f.write(f"CLASS ERROR for params {cfg}: {e}\n")
        raise RuntimeError(f"CLASS failed: {e}")

    derived = G.get_current_derived_parameters(['z_reio','tau_reio','conformal_age'])
    T_cmb = G.T_cmb()
    return {
        'T_cmb': T_cmb,
        'z_reio': derived['z_reio'],
        'tau_reio': derived['tau_reio'],
        'conformal_age': derived['conformal_age'],
        'G': G
    }


# =======================================
# Cobaya Theory Class
# Provides an interface for Cobaya to run the DH+CLASS calculation and retrieve results (Cls and derived parameters).
# =======================================
class DHClassTheory(Theory):

    # Define parameters that Cobaya will sample.
    params = {
        "H0": None,
        "Omega_cdm": None,
        "Omega_b": None,
        "log10_mnu": None,
        "mDM": None,
        "lifetime_exponent": None,
        "lna_pivot": None,
        "tilt": None
    }

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
            'ell':cl_dict['ell']
        }       

        state["derived"] = {
            "z_reio": derived["z_reio"],
            "tau_reio": derived["tau_reio"],
            "conformal_age": derived["conformal_age"]
        }
        
        return state

    def get_can_provide_params(self):
        return [
            "z_reio",
            "tau_reio",
            "conformal_age"
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
            units_power = float(sum(cl_name.count(p) for p in ["t","e","b"]))
            # Apply units_factor
            cls[cl_name][2:] *= units_factor**units_power

            # Apply ell_factor
            if ell_factor:
                if "p" not in cl_name:
                    cls[cl_name][2:] *= ell_factor_array
                elif cl_name == "pp":
                    cls[cl_name][2:] *= ell_factor_array**2 * (2 * np.pi)
                elif "p" in cl_name:
                    cls[cl_name][2:] *= ell_factor_array**(3/2) * np.sqrt(2 * np.pi)

        return cls


if __name__ == "__main__":
    args = parse_args()
    cfg = load_config(args)

    run_dh_class(cfg)
