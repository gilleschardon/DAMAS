# DAMAS and CMF-NNLS

Implementation of CMF-NNLS and DAMAS-NNLS, python and Matlab

This code reproduces the results of
Theoretical analysis of the DAMAS algorithm and efficient implementation for large-scale problems


## Demos

Simple examples : setup of the problem, inversion, plots

* demo_cmf_nnls.m
* demo_cmf_nnls.py

## Figures

Reproduction of the figures of the paper

* FIGS_2D_EXP.m
* FIGS_3D_EXP.m
* FIGS_2D_SIM.m (this code is non deterministic, don't be surprised to get slightly different figures than in the paper)

## Matlab code

* cmf_nnls.m: basic CMF-NNLS (can be used for low-noise level and noise level estimation)
* cmf_nnls_dr.m: CMF-NNLS with diagonal removal (same results as with noise level estimation, slightly faster)
* dictionary.m: builds the matrix D
* proddamas.m, proddamasDR.m, proddamastranspose.m: efficient implementations of matrix products
* damas.c, damas_rand.c: mex implementation of the DAMAS algorithm (don't use it! CMF-NNLS is faster)
* damas_nnls.m: DAMAS-NNLS

## Python code

* damas.py: implementation of CMF-NNLS and DAMAS-NNLS (standard and with diagonal removal)