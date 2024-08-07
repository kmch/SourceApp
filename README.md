# :pill: SourceApp :non-potable_water:
This repository contains material accompanying the paper by Chrapkiewicz et al. (2024) [Apportioning sources of chemicals of emerging concern along an urban river with inverse modelling](https://doi.org/10.1016/j.scitotenv.2024.172827), published in Science of the Total Environment. 

![graphical_abstract](https://github.com/kmch/SourceApp/assets/19067370/81e1b40c-6bb9-4b74-aa37-ce7c884980f3)

At the centre of our workflow is the `faster-unmixer` software by [Barnes \& Lipp (2024)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023WR036159). It can be installed following the guidelines below.

## Structure of the repository
- [data](data): contains publicly available data sufficient to reproduce the results,
- [examples](examples): contains a minimal working example, MWE (see 'Installation' below for how to run it).

The example uses two input files (provided under `data`):
- flow directions across Thames basin (in ESRI D8 format) 
-  an Excel file with concentrations of contaminants of emerging concern measured across Thames basin in London by Egli et al.  (2023), see: https://doi.org/10.1016/j.envint.2023.108210.

## Installation
To install the software required to run the MWE, follow the steps below. The commands are supposed to be run in a Unix-like terminal (tested in `bash` and `zsh` shells), unless otherwise specified. Make sure you have `anaconda` installed before you proceed. 

1. Clone the repository and navigate to its directory:
```bash
git clone git@github.com:kmch/SourceApp.git
cd SourceApp
```

2. Create a conda environment using the provided yaml file:
```bash
my_env=mwe; conda env create --name $my_env --file=mwe.yaml && conda activate $my_env
```
From now on, all the commands need to be run in this conda environment. Make sure that the `python` command points to the executable of this environment, not the system-wide Python installation. Check this by running:
```bash
which python
```
It should output a path looking more or less like this:
```bash
$HOME/anaconda3/envs/mwe/bin
```
where $HOME will be your home path.

3. Install the `autocatchments` package, following the guidelines on https://github.com/AlexLipp/autocatchments/tree/name listed below:
```bash
git clone git@github.com:AlexLipp/autocatchments.git
cd autocatchments
pip install -e .
```

4. Install the `faster-unmixer` package:
```bash
git clone --recurse-submodules https://github.com/r-barnes/faster-unmixer/ 
cd faster-unmixer
git submodule update --init --recursive
pip install -e .
```

5. Run the MWE:
```bash
cd SourceApp/examples
python mwe.py
```
Open the generated PNG files. The top-right subplot of each 6-subplot figure is supposed to be blank.

### Troubleshooting
If you experience troubles with Python failing to import some modules, check if they are installed correctly by running:
```bash
conda activate $my_env && conda list | grep riverpy
```
where `riverpy` should be replaced by the module you want to check. If the above command returns an empty line, this means the module is not installed in your environment. To solve this problem, make sure that the `python` command points to the executable of this environment (see the step 2. of this section).

If you encounter any run-time errors from the calls to `funmixer`, it may mean the version you cloned has some modifications that made it incompatible with the rest of the software. In this case, change to the tested version of the package by running:
```bash
cd faster-unmixer
git checkout 939c9f89f27ccefc4a172eb065c0d6f1f0de0438
```
 
## Citing
If you use any part of this repository in your research, please cite the following paper:

-  Chrapkiewicz, Kajetan and Lipp, Alex and Barron, Leon Patrick and Barnes, Richard and Roberts, Gareth, *Apportioning sources of chemicals of emerging concern along an urban river with inverse modelling*. Science of the Total Environment. [DOI:10.1016/j.scitotenv.2024.172827](https://doi.org/10.1016/j.scitotenv.2024.172827).
