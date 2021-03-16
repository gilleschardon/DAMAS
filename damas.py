#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 14:36:35 2020

@author: gilleschardon
"""

import numpy as np
import scipy.linalg as la
import time
import matplotlib.pyplot as plt

# some fast matrix products

# beamforming
def proddamastranspose(D, Data):    
    x = np.real( np.sum( (D.conj().T @ Data) * D.T, 1))    
    return x

# product by the DAMAS matrix
def proddamas(D, x, support=None):    
    if support is None:
        z = D @ (x * D.conj()).T
    else:
        z = D[:, support] @ (x[support] * D[:, support].conj()).T    
    return proddamastranspose(D, z)

def proddamasdr(D, x, support=None):    
    if support is None:
        z = D @ (x * D.conj()).T
    else:
        z = D[:, support] @ (x[support] * D[:, support].conj()).T
    z = z - np.diag(np.diag(z))    
    return proddamastranspose(D, z)
    
# local unconstrained least-squares problem
def solve_ls(D, bf, support):
    Gram = np.abs(D[:, support].conj().T @ D[:, support]) ** 2
    return la.solve(Gram, bf[support], assume_a='pos'), Gram

# local unconstrained least-squares problem for diagonal removal
def solve_ls_dr(D, bf, support):
    aD2 = np.abs(D[:, support] ** 2)
    Gram = np.abs(D[:, support].conj().T @ D[:, support]) ** 2 - aD2.T @ aD2
    return la.solve(Gram, bf[support], assume_a='pos'), Gram

# same, for DAMAS
def solve_damas_ls(D, b, support):    
    AR = np.abs(D.conj().T @ D[:, support]) ** 2
    Gram = AR.T @ AR
    return la.solve(Gram, b[support], assume_a='pos'), Gram

def solve_damas_ls_dr(D, b, support):
    DR = D[:, support]        
    aDR2 = np.abs(D)**2
    aDR2r = np.abs(DR)**2
    AR = np.abs(D.conj().T @ DR)**2 - aDR2.T @ aDR2r
    Gram = AR.T @ AR
    return la.solve(Gram, b[support], assume_a='pos'), Gram

## CMF-NNLS
def cmf_nnls_lh(D, Data):
    
    prodA = lambda D, x, support : proddamas(D, x, support)
    solve = lambda D, b, support : solve_ls(D, b, support)
    
    bf = proddamastranspose(D, Data)

    return lawson_hanson(prodA, solve, D, bf)
        
## CMF-NNLS-DR
def cmf_nnls_dr_lh(D, Data):
    
    Data = Data - np.diag(np.diag(Data))
  
    prodA = lambda D, x, support : proddamasdr(D, x, support)
    solve = lambda D, b, support : solve_ls_dr(D, b, support)
    
    bf = proddamastranspose(D, Data)

    return lawson_hanson(prodA, solve, D, bf)

## DAMAS-NNLS
def damas_nnls_lh(D, Data):
    
    
    bf = proddamastranspose(D, Data)

    b = proddamas(D, bf)
    
    prodA = lambda D, x, support : proddamas(D, proddamas(D, x, support))
    solve = lambda D, b, support : solve_damas_ls(D, b, support)
    

    return lawson_hanson(prodA, solve, D, b)
  
## DAMAS-NNLS-DR      
def damas_nnls_dr_lh(D, Data):
    
    DataDR = Data - np.diag(np.diag(Data))

    bf = proddamastranspose(D, DataDR)

    b = proddamasdr(D, bf)
    
    prodA = lambda D, x, support : proddamasdr(D, proddamasdr(D, x, support))
    solve = lambda D, b, support : solve_damas_ls_dr(D, b, support)
    

    return lawson_hanson(prodA, solve, D, b)
        


## Lawson-Hanson solver, custom implementation

# ||Mx - y||_2^2, with M described by D
# prodA(D, x, support) = M*x
# solve(D, b, support) solve the local LS problem
# b = M'*y
def lawson_hanson(prodA, solve, D, b, verbose = True):
    
    T0 = time.perf_counter()
    
    n = D.shape[1]
    R = np.ones([n], dtype=bool) # active set
    N = np.arange(n);
    x = np.zeros([n])
    
    # residual
    
    w = b
    
    it = 0
    
    while np.any(R) and (np.max(w[R]) > 0):
        if verbose:
            print(f"iter {it} tol {np.max(w[R]):.2}")
        it = it + 1

        # update of the active set        
        idx = np.argmax(w[R])
        Ridx = N[R]
        idx = Ridx[idx]        
        R[idx] = 0
        
        # least-square in the passive set        
        s = np.zeros(x.shape)      
        s[~R], Gram = solve(D, b, ~R)
        
        # removal of negative coefficients
        while np.min(s[~R]) <= 0:
            
            Q = (s <= 0) & (~R)
            alpha = np.min(x[Q] / (x[Q] - s[Q]))
            x = x + alpha * (s - x)
            R = (( x <= np.finfo(float).eps) & ~R) | R
            
            s = np.zeros(x.shape)
            s[~R], Gram = solve(D, b,  ~R)

            
        # update of the solution
        x = s
        # update of the residual
        w = b - prodA(D, x, ~R)

        
    try:
        la.cholesky(Gram)
        if verbose:
            unique = True
            print(f'Solution is unique, T = {time.perf_counter() - T0:.2}')
    except e:
        if verbose:
            unique = False
            print(f'Solution is not unique, T = {time.perf_counter() - T0:.2}')
        
            
    return x, unique



def beamforming_noweight(D, Data):

    return np.real(  np.sum((D.conj().T @ Data) * D.T, 1))

def beamforming_weight(D, Data):
    Dbf = D / np.sum(np.abs(D)**2, 0);
    return np.real(  np.sum((Dbf.conj().T @ Data) * Dbf.T, 1))

def orthogonal_beamforming(D, Data, K):
    Dbf = D / np.sum(np.abs(D)**2, 0);
    w, v = la.eigh(Data)    
    obf = np.zeros([D.shape[1]])
    for u in range(K):
        l = w[-u-1]
        vec = v[:, -u-1]
        bfvec = l * np.abs(Dbf.conj().T @ vec)
        idx = np.argmax(bfvec)
        obf[idx] = obf[idx] + bfvec[idx]
    return obf

def CLEAN(D, Data, K, phi):
    Dbf = D / np.sum(np.abs(D)**2, 0);
    
    Data = np.copy(Data)

    cbf = np.zeros([D.shape[1]])
    
    for u in range(K):
        BF = beamforming_noweight(Dbf, Data)
        idx = np.argmax(BF)
        a = phi * BF[idx]
        
        cbf[idx] = cbf[idx] + a
        Data = Data - a * D[:, idx:idx+1] @ D[:, idx:idx+1].conj().T
    
    return cbf

def CLEANSC(D, Data, K, phi):
    T0 = time.perf_counter()

    Dbf = D / np.sum(np.abs(D)**2, 0);
    
    Data = np.copy(Data)
    
    cbf = np.zeros([D.shape[1]])
    
    for u in range(K):
        BF = beamforming_noweight(Dbf, Data)
        idx = np.argmax(BF)
        p = BF[idx]
        a = phi * p
        cbf[idx] = cbf[idx] + a
        
        h = Data @ Dbf[:, idx:idx+1] / p
        Data = Data - a * h @ h.conj().T
    print(f'T = {time.perf_counter() - T0:.2}')

    return cbf
