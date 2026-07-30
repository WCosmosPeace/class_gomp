# `DarkHistory` Integration for Decaying Dark Matter

This branch, `DH_CLASS_integration`, extends `class_gomp` to include the ionization and thermal histories produced by decaying dark matter using `DarkHistory`.

The repository provides:

* a Jupyter notebook demonstrating the complete
  `DarkHistory → class_gomp → CMB spectra` workflow;
* example `Cobaya` configuration files for parameter inference with the extended model.

The required history table is generated automatically by the notebook and cached for subsequent runs.

## Repository contents

### Jupyter notebook

```text
notebooks/Gomp_DM_decay_CMB_example.ipynb
```

The notebook demonstrates how to:

1. define a decaying-dark-matter model and a Gomp reionization history;
2. run `DarkHistory` to calculate the ionization and thermal history;
3. save or reuse the resulting `DarkHistory` cache;
4. pass the corresponding model parameters to `class_gomp`;
5. calculate the CMB angular power spectra;
6. compare the Gomp and Gomp + DM-decay models.

The included example uses the decay channel $\chi\rightarrow e^+e^-$, with $m_\chi=1~\mathrm{GeV}$.

### Cobaya configurations

The Cobaya configuration files are stored in:

```text
cobaya/
```

These files provide example likelihood, parameter, prior, and sampler settings for parameter inference with `class_gomp`. External likelihood data and machine-specific paths must be configured before use.

The Cobaya configurations are not required to run the Jupyter notebook.

## Directory structure

The integration expects the `class_gomp` and `DarkHistory` directories to be located at the same level:

```text
project-directory/
├── class_gomp/
│   ├── notebooks/
│   │   └── Gomp_DM_decay_CMB_example.ipynb
│   └── cobaya/
└── DarkHistory/
    ├── darkhistory/
    ├── data/
    │   ├── decay/
    │   │   └── decay_electron_gamma.txt
    │   └── README.md
    └── decay_e/
```

The `data/` directory contains the DarkHistory transfer-function data.
The `data/decay/decay_electron_gamma.txt` file contains the secondary photon injection spectra produced by dark matter decay into electron–positron pairs.

The `decay_e/` directory stores the ionization and thermal histories generated for the decay channel $\chi \rightarrow e^+e^-$. 
A valid cached history is automatically reused when the same model parameters are requested again.

## Environment setup

Set `DH_DATA_DIR` to the DarkHistory data directory before starting Jupyter:

```bash
export DH_DATA_DIR=/absolute/path/to/DarkHistory/data
```

For example, if the two repositories are located under `/home/user/project/`, use:

```bash
export DH_DATA_DIR=/home/user/project/DarkHistory/data
```

The environment variable must be set before the Jupyter kernel starts. After changing it, restart the kernel.

## Compile the Python wrapper

Compile the `classy` Python wrapper from this `class_gomp` branch and make sure that it is importable from the Python environment used by Jupyter and Cobaya.

The following check should import the wrapper compiled from this repository rather than a standard CLASS installation:

```bash
python -c "from classy import Class; import classy; print(classy.__file__)"
```

Also verify the compatible DarkHistory installation:

```bash
python -c "import darkhistory; print(darkhistory.__file__)"
```

## Notes

* The notebook compares the two models at their respective posterior-mean parameter values. Their spectral difference therefore contains both the direct DM-decay contribution and shifts in the fitted cosmological and Gomp parameters.
* Public Planck PR3 spectra are used only for visualization in the notebook.
* DarkHistory cache files are generated locally and reused when available; they are not distributed with the repository.
* Cobaya configuration files and run scripts are provided, while generated chains, checkpoints, logs, and scheduler outputs are not included.