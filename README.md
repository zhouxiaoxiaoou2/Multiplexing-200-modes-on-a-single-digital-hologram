# Reproduction of "Multiplexing 200 modes on a single digital hologram" 

This repository contains my reproduction of the experiments and results from：
 (Rosales-Guzmán C, Bhebhe N, Mahonisi N, et al. Multiplexing 200 spatial modes with a single hologram[J]. Journal of Optics, 2017, 19(11): 113501.)

## Contents
- Implementation of model described in the paper
- Training scripts

## Notes
This is a personal reproduction project for learning and research purposes.  
All credits for the original idea and experiments go to the authors of the paper.  

## What This Paper is About
The core idea of the paper is:
- Objective: Encode and superpose many different optical field modes simultaneously within a single digital hologram.

- Significance: Greatly increases the capacity for optical information transmission and processing (applications in optical communications and optical computing).

- Mode Types: Focuses on Hermite–Gaussian (HG) and Laguerre–Gaussian (LG) modes. Both are exact solutions to the laser field equations and can act as independent information channels.

- Method: By using a Spatial Light Modulator (SLM) or computer-generated hologram (CGH), different modes can be encoded, generated, and detected simultaneously.

## Example: Reproducing Figure 8 MATLAB codes
- Defines a spatial grid (X,Y) and polar coordinates (ρ,φ).
- Implements Hermite–Gaussian (HG) and Laguerre–Gaussian (LG) mode functions.
- Randomly generates 42 modes (half HG, half LG) with varied parameters.
- Arranges them into a 6×7 multiplexed grid.
- Displays the multiplexed image, illustrating how multiple modes can coexist

## Example output
[Click here to view Figure 8 result](results/figure8.png)