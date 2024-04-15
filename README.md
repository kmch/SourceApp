# :pill: :non-potable_water: SourceApp
This repository contains material accompanying the manuscript  *Apportioning sources of chemicals of emerging concern along an
urban river with inverse modelling* by Chrapkiewicz et al. (submitted to Science of the Total Environment), preprint available: [link](https://doi.org/10.31223/X52T22)


## Structure of the repository
- [data](data): contains publicly available data sufficient to reproduce the results,
- [examples](examples): contains a minimal working example, MWE (see 'Installation' below for how to run it).

## Installation
To install the software required to run the MWE, follow the steps below. The commands are to supposed to be run in a Unix-like terminal, unless otherwise specified. Make sure you have anaconda installed before you proceed. 

1. Clone the repository and navigate to its directory:
```bash
git clone git@github.com:kmch/SourceApp.git
cd SourceApp
```

2. Create a conda environment using the provided yaml file:
```bash
my_env=mwe; conda env create --name $my_env --file=mwe.yaml && conda activate $my_env
```
From now on, all the commands need to be run in this conda environment.

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
Open the generated png files.

## Citing
If you use any part of this repository in your research, please cite the following paper:

-  Chrapkiewicz, Kajetan and Lipp, Alex and Barron, Leon Patrick and Barnes, Richard and Roberts, Gareth, *Apportioning sources of chemicals of emerging concern along an urban river with inverse modelling*. Available at SSRN: https://ssrn.com/abstract=4697195 or http://dx.doi.org/10.2139/ssrn.4697195.
