# :pill: :non-potable_water: Wandle-App
This repository contains material accompanying the manuscript  *Apportioning sources of chemicals of emerging concern along an
urban river with inverse modelling* by Chrapkiewicz et al. (submitted to Science of the Total Environment), preprint available: [link](https://doi.org/10.31223/X52T22)


## Structure of the repository
- [data](data): contains publicly available data necessary to reproduce most of the results,
- [riverpy](riverpy): (**switch to branch riverpy**) contains the open-source framework of the data analysis,
- [notebooks](notebooks): contains Jupyter notebooks with the workflow of data analysis.

## Installing `riverpy`
Make sure to first install `jupconfig`, `plotea`, `iogeo` repositories of `kmch`'s and `autocatchments`, `compysitional`, `faster-unmixer` and `thames-sewage` of Alex Lipp's.


```python
conda create -n riverpy python=3.10
```

## Citing
If you use any part of this repository in your research, please cite the following paper:

-  Chrapkiewicz K, Lipp AG, Barron LP, Barnes R, and Roberts GG, *Apportioning Sources of Chemicals of Emerging Concern Along an Urban River with Inverse Modelling*. Submitted to Science of the Total Environment. Preprint available at https://doi.org/10.31223/X52T22 and http://dx.doi.org/10.2139/ssrn.4697195.
