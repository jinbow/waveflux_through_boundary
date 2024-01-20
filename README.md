# Diagnose wave energy flux from MITgcm simulations

Some codes originally used for the project led to Mazloff et al. 2020. This code is not yet meant to be reproducible but only as reference for algorithm. Generalizing it is on the to-do list. 

Mazloff, M. R., Cornuelle, B.,Gille, S. T., & Wang, J. (2020). Theimportance of remote forcing forregional modeling of internal waves.Journal of Geophysical Research:Oceans,125, e2019JC015623. https://doi.org/10.1029/2019JC015623

### Calculation Procedures
The algorithm follows Nash, J. D., M. H. Alford, and E. Kunze, 2005: Estimating Internal Wave Energy Fluxes in the Ocean. J. Atmos. Oceanic Technol., 22, 1551–1570, https://doi.org/10.1175/JTECH1784.1. 

#### 1. Density and Pressure Calculation:

- **Density Calculation**: 
  - The potential density is calculated from salinity (S) and temperature (T) data. Using the JMD95 equation of state referenced to p=0 $\rho=\text{densjmd95}(S, T, 0)$
  - densjmd95 comes from MITgcm utils: http://mitgcm.org/download/daily_snapshot/MITgcm/utils/python/MITgcmutils/MITgcmutils/jmd95.py
- **Pressure Calculation**: 
  - The pressure at each depth level is often calculated by integrating the weight of the water column above it. This is represented by the formula: $$P_n = \sum_{i=0}^{n} \rho_i \times \Delta z_i$$ MITgcm can also save hydrostatic pressure as a diagnostic (PHIHYD). Using PHIHYD*rhoconst gives you pressure anomaly (deviation from rhoRef(z)) that can be directly used in the flux calculation. 
#### 2. Perturbation Calculation:
- Perturbations in variables like pressure and velocity are computed by removing the temporal linear trend and a low-frequency component. One can use band-pass filter to signal out a certain tidal band. The procedure typically involves:
  - Removing the temporal linear trend.
  - Applying a high-pass (band-pass) filter to remove mesoscale variations.
#### 3. Computing Flux:
  - The flux is computed as the product of the velocity and pressure perturbations $u'$ and $p'$, respectively.
    $$Flux = u'p'$$ on level n.
  - The flux is scaled by the vertical grid spacing (DRF) and the horizontal grid spacing (DXG for meridional or DYG for zonal), multiplied by the cell fraction (hFacW, hFacS, or hFacC): $$Flux=Flux \times DRF \times (\text{DXG or DYG}) \times \text{hFacW or hFacS}$$
