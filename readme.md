# EC_NDJ_Dynamics

Computational framework for Endothelial Cell (EC) Notch-Delta-Jagged (NDJ) dynamics, combining spatial CompuCell3D simulations with intracellular signaling models and parameter search workflows.

## What this repository contains

- **CC3D simulation model** for EC behavior and interactions.
- **NDJ intracellular model** (Antimony/SBML based).
- **Search/optimization scripts** for parameter exploration.
- **Documentation** for model rationale and parameters.

## Project layout

- Core orchestration: `EC_Connect_core.py`
- Parameter search workflows: `EC_Connect_search.py`
- NDJ signaling model definition: `NED_HES_Antimony.py`
- Environment specifications: `config_cc3d.yml`, `config_tellurium.yml`
- CC3D project and steppables: `CC3D_NED/`
- Documentation and model notes: `docs/`
- Generated figures and simulation outputs: `nb_Figures/`

## Quick start

### 1) Create environment

For full CC3D workflows:

```bash
conda env create -f config_cc3d.yml
conda activate NDJ_CC3D
```

For SBML/Tellurium-only workflows:

```bash
conda env create -f config_tellurium.yml
conda activate NDJ_Dynamics_tellurium
```

### 2) Run a simulation/search from Python

```python
from EC_Connect_search import simple_sim_run, default_params

store_path = "path/to/results.zarr"
simple_sim_run(store_path, default_params, "experiment_name")
```

### 3) Open the CC3D project

Open `CC3D_NED/EC_Connect_v11.cc3d` in CompuCell3D Player/Twedit++.

## Documentation

Detailed model notes and parameter references are in `docs/`:

- [NED buildup notebook explanations](docs/NED_buildup_nb_explanations_V2.0.md)
- [NED parameters table](docs/NED_buildup_parmeters_table.md)
- [HES1 autoregulation sensitivity analysis (PDF)](docs/HES1_Autoregulation_SA.pdf)

## Data

All original data generated for this project is available on Zenodo: [https://doi.org/10.5281/zenodo.18865760](https://doi.org/10.5281/zenodo.18865760)

## Citation

If you use this model, please cite:

**Chesnais et al.**
*A spatialised agent-based model of NOTCH signalling pathway in Endothelial Cells predicts emergent heterogeneity due to continual dynamic phenotypic adjustments.*
bioRxiv (2022). DOI: [10.1101/2022.08.06.503043](https://www.biorxiv.org/content/10.1101/2022.08.06.503043v1)

## License

Open source under the MIT License. See `LICENSE` for details.