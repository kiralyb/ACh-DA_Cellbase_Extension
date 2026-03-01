# ACh-DA_Cellbase_Extension
Matlab analyses for cholinergic and dopaminergic single unit and fiber photometry recordings

## Description
This repository contains MATLAB code for data analysis on the interaction of the cholinergic and dopaminergic systems in broadcasting prediction error signals during associative learning. The figures and supplementary figures of 'Cholinergic–dopaminergic interplay underlies prediction error broadcasting' by Király et al. can be produced with the MAIN_ACH_DA function relying on three datasets recorded from mice performing an auditory go/no-go task:

  - Fiber photometry recordings of acetylcholine (ACh) and dopamine (DA) release from the basolateral amygdala (BLA) and the ventral stiatum (VS).

  - Single unit recordings from the horizontal limb of the diagonal band of Broca (HDB) and the ventral tegmental area (VTA) combined with the optogenetical tagging of basal forebrain cholinergic neurons (BFCNs) and midbrain dopaminergic neurons (DANs).

  - Fiber photometry recordings of ACh and DA release from the BLA and the VS during the chemogenetic supression of BFCNs.

## Content

- .m files for data analysis
- license file

## Installation

Move the .m files on your MATLAB path (should take around a few seconds). No further installation is needed, but Cellbase needs to be initilized (initcb.m). 

## Dependencies

This code is an extension to CELLBASE (https://github.com/hangyabalazs/CellBase).

## System requirements
The code was tested on the following configuration:
- Windows 10 64bits
- Intel i7
- 32 GB RAM
