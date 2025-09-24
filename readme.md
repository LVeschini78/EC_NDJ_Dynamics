# EC_NDJ_Dynamics

A computational framework for modeling Endothelial Cell (EC) Notch-Delta-Jagged (NDJ) dynamics using CompuCell3D and systems biology approaches.

## Overview

This repository contains a multi-scale computational model that simulates the spatiotemporal dynamics of Notch signaling in endothelial cells. The model integrates:

- **Cellular-level dynamics** using CompuCell3D (CC3D) for spatial cell-based modeling
- **Molecular-level dynamics** using SBML/Antimony for intracellular signaling networks
- **Parameter optimization** using evolutionary algorithms (Nevergrad)
- **Statistical validation** against experimental HUVEC data

## Key Features

### 1. Multi-Scale Modeling
- **Spatial cellular dynamics**: Cell morphology, division, adhesion, and interactions
- **Intracellular signaling**: Notch-Delta-Jagged pathway with HES1 autoregulation
- **Cis/trans interactions**: Both lateral inhibition and self-inhibition mechanisms

### 2. Notch Signaling Components
- **Receptors**: NOTCH1, NOTCH4
- **Ligands**: DLL4 (Delta-like 4), JAG1 (Jagged1)
- **Transcription factors**: HES1 with autoregulatory feedback
- **NICD**: Notch intracellular domain dynamics

### 3. Parameter Optimization
- **Nevergrad integration**: Differential Evolution optimizer
- **Parallel execution**: Multi-core parameter search capabilities
- **Statistical validation**: Kolmogorov-Smirnov test against HUVEC experimental data

## Repository Structure

```
EC_NDJ_Dynamics/
├── config_cc3d.yml           # Conda environment for CC3D simulations
├── config_tellurium.yml      # Conda environment for SBML modeling
├── EC_Connect_core.py         # Core simulation engine
├── EC_Connect_search.py       # Parameter optimization framework
├── NED_HES_Antimony.py       # SBML model builder for Notch pathway
├── NED_buildup_nb.ipynb      # Jupyter notebook for model development
├── huvec_data.json           # Experimental HUVEC reference data
├── CC3D_NED/                 # CompuCell3D simulation files
│   ├── EC_Connect_v11.cc3d   # CC3D project file
│   └── Simulation/           # CC3D steppables and configurations
│       ├── EC_Connect_V10.py # Main CC3D simulation script
│       ├── EC_Connect_V10.xml # CC3D configuration
│       ├── NDJ_V9.py        # Notch-Delta-Jagged model parameters
│       ├── NDJ_SBML_Step.py # SBML integration steppable
│       ├── NDJ_Interactions.py # Cell-cell interaction logic
│       ├── NDJ_CellMorph.py # Cell morphology dynamics
│       ├── NDJ_Vis_DataTrack.py # Data visualization and tracking
│       └── *.piff           # Cell shape and pattern files
├── nb_Figures/              # Generated analysis figures
├── Output_CC3D/            # Simulation output data
└── LICENSE                 # MIT License
```

## Installation & Requirements

### Option 1: Full CC3D Environment (Recommended)
For complete functionality including spatial cellular modeling:

```bash
# Clone the repository
git clone https://github.com/LVeschini78/EC_NDJ_Dynamics.git
cd EC_NDJ_Dynamics

# Create and activate CC3D environment
conda env create -f config_cc3d.yml
conda activate NDJ_CC3D
```

This environment includes:
- **CompuCell3D 4.6**: Spatial cell-based modeling framework
- **Nevergrad**: Parameter optimization with evolutionary algorithms  
- **Zarr**: High-performance array storage
- **Scientific Python stack**: NumPy, pandas, matplotlib, scikit-image

### Option 2: Tellurium-Only Environment
For SBML modeling and analysis without spatial simulations:

```bash
# Create and activate Tellurium environment
conda env create -f config_tellurium.yml
conda activate NDJ_Dynamics_tellurium
```

This lighter environment includes:
- **Tellurium**: SBML/Antimony modeling and simulation
- **Scientific Python stack**: NumPy, pandas, matplotlib, scikit-image

### Environment Notes
- **Python 3.10**: Required for CompuCell3D compatibility
- **Zarr storage**: Optimized for large simulation datasets
- **Multiprocessing**: Enabled for parallel parameter searches

## Usage

### 1. Running Single Simulations
```python
from EC_Connect_search import simple_sim_run, default_params

# Run with default parameters
store_path = "path/to/your/output.zarr"
simple_sim_run(store_path, default_params, "experiment_name")
```

### 2. Parameter Optimization
```python
from EC_Connect_search import run_parallel_optimization

# Run evolutionary parameter search
store_path = "path/to/your/output.zarr"
run_parallel_optimization(store_path)
```

### 3. Two-Parameter Grid Search
```python
from EC_Connect_search import run_two_params_search_parallel

# Search over DLL4 production (ps_bd4) and Jagged inducibility (ps_ind)
run_two_params_search_parallel(
    store_path,
    param_1='ps_bd4', range_1=(50, 150),
    param_2='ps_ind', range_2=(0, 100),
    n_points=10,
    max_workers=4
)
```

### 4. Interactive Model Development
Open `NED_buildup_nb.ipynb` in Jupyter for:
- Model parameter exploration
- SBML model visualization  
- Results analysis and plotting

## Model Parameters

### Key Biological Parameters
- `cr`: Cis interaction ratio (lateral inhibition strength)
- `ps_bd4`: Basal DLL4 production rate
- `ps_bn1`: Basal NOTCH1 production rate  
- `ps_bj1`: Basal JAG1 production rate
- `ps_ind`: JAG1 inducibility by HES1
- `ps_ja`: JAG1 induction threshold
- `ps_ni`: NICD production threshold
- `ps_Kdni`: NICD degradation rate
- `ps_kconv`: Cis interaction conversion efficiency

## Experimental Validation

The model is validated against HUVEC (Human Umbilical Vein Endothelial Cells) experimental data:
- **Reference data**: `huvec_data.json` contains normalized expression measurements
- **Statistical test**: Kolmogorov-Smirnov test for distribution comparison
- **Optimization objective**: Minimize KS statistic between simulation and experimental data

## Key Publications & Methods

This work builds upon established computational approaches for:
- **Notch signaling modeling**: Collier et al., Sprinzak et al.
- **CompuCell3D framework**: Swat et al.
- **Evolutionary optimization**: Nevergrad (Facebook Research)

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions about the model implementation or collaboration opportunities, please open an issue or contact the repository maintainer.

## Acknowledgments

- CompuCell3D development team
- Tellurium/Antimony SBML community
- Facebook Research Nevergrad team

## Citation

If you use this code or model in your research, please cite:

**Chesnais et al.** "A spatialised agent-based model of NOTCH signalling pathway in Endothelial Cells predicts emergent heterogeneity due to continual dynamic phenotypic adjustments." *bioRxiv* (2022). doi: [10.1101/2022.08.06.503043](https://www.biorxiv.org/content/10.1101/2022.08.06.503043v1)

```bibtex
@article{chesnais2022spatialised,
  title={A spatialised agent-based model of NOTCH signalling pathway in Endothelial Cells predicts emergent heterogeneity due to continual dynamic phenotypic adjustments},
  author={Chesnais, L and others},
  journal={bioRxiv},
  pages={2022--08},
  year={2022},
  publisher={Cold Spring Harbor Laboratory},
  doi={10.1101/2022.08.06.503043},
  url={https://www.biorxiv.org/content/10.1101/2022.08.06.503043v1}
}
```