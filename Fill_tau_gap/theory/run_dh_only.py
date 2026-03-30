# run_dh_only.py
import os
import sys
import json
import traceback
import numpy as np

base_data_dir = os.environ.get("DH_DATA_DIR")
if base_data_dir is None:
    raise EnvironmentError("Please set DH_DATA_DIR")

base_dir = os.path.dirname(base_data_dir)
sys.path.append(base_dir)

logpath = os.path.join(os.getcwd(), "dh_errors.log")

import darkhistory.physics as phys
import darkhistory.main as main
from darkhistory.config import load_data


def main_entry(cfg_json):
    cfg = json.loads(cfg_json)

    h = cfg["H0"] / 100.0
    omega_b = cfg["omega_b"]
    omega_cdm = cfg["omega_cdm"]

    mnu = max(float(cfg["mnu"]), 1e-5)
    if not np.isfinite(mnu):
        raise ValueError(f"Non-finite mnu in run_dh_only: {cfg['mnu']}")

    cosmo_params = {
        "H0": cfg["H0"],
        "Omega_b": omega_b / h**2,
        "Omega_cdm": omega_cdm / h**2,
        "mnu": mnu,
    }

    phys.set_cosmology_params(cosmo_params)

    mDM = cfg["mDM"]
    lifetime = 10 ** cfg["lifetime_exponent"]
    lna_pivot = cfg["alpha_gomp"]
    tilt = cfg["beta_gomp"]
    inj_particle = cfg.get("inj_particle", "decay_e")

    M = mDM / 1e9
    tau_str = f"{cfg['lifetime_exponent']:.6f}"

    decay_dir = os.path.join(base_dir, inj_particle)
    os.makedirs(decay_dir, exist_ok=True)

    filename = os.path.join(
        decay_dir,
        f"{M:.6f}_{tau_str}_{cosmo_params['H0']:.6f}_"
        f"{cosmo_params['Omega_cdm']:.6f}_{cosmo_params['Omega_b']:.6f}_"
        f"{cosmo_params['mnu']:.6f}_{tilt:.6f}_{lna_pivot:.6f}.txt"
    )

    if os.path.exists(filename):
        with open(filename, "r") as f:
            lines = f.readlines()

        if len(lines) == 224 and not any("nan" in line.lower() for line in lines):
            return 0
        else:
            print(f"{filename} corrupted (lines={len(lines)}). Recomputing.")
            os.remove(filename)

    main.evolve(
        mDM=mDM,
        DM_process="decay",
        lifetime=lifetime,
        primary=inj_particle,
        start_rs=3000,
        end_rs=4,
        coarsen_factor=30,
        backreaction=True,
        reion_switch=True,
        reion_rs=800,
        lna_pivot=lna_pivot,
        tilt=tilt
    )

    return 0


if __name__ == "__main__":
    try:
        cfg_json = sys.argv[1]
        main_entry(cfg_json)
    except Exception:
        with open(logpath, "a") as f:
            f.write(traceback.format_exc())
        sys.exit(1)

