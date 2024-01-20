# Diagnose wave energy flux from MITgcm simulatons

Some codes used for the project in Mazloff et al. 2020. This code is not meant to reproducible but only as reference for algorithm.  

Mazloff, M. R., Cornuelle, B.,Gille, S. T., & Wang, J. (2020). Theimportance of remote forcing forregional modeling of internal waves.Journal of Geophysical Research:Oceans,125, e2019JC015623. https://doi.org/10.1029/2019JC015623

### Calculation Procedures
Follows Nash et al. (2005)
#### 1. Density and Pressure Calculation:
- **Density Calculation**: 
  - The potential density is typically calculated from salinity (S) and temperature (T) data. Using the JMD95 equation of state, it is given by:
  
    \[ \text{Density} = \text{function}(S, T, \text{pressure}) \]

- **Pressure Calculation**: 
  - The pressure at each depth level is often calculated by integrating the weight of the water column above it. This is represented by the formula:

    \[ \text{Pressure at level } n = \text{constant} \times \sum_{i=0}^{n} (\text{Density at level } i \times \text{cell thickness at level } i) \]

#### 2. Perturbation Calculation:
- Perturbations in variables like pressure and velocity are computed by removing the mean or background state from the actual state. The procedure typically involves:
  - Removing the temporal linear trend.
  - Applying a high-pass filter to remove mesoscale variations.

#### 3. Energy Flux Calculation:
- **Velocity and Pressure Perturbations**: 
  - Let \( U_p \) and \( P_p \) represent the perturbations in the velocity and pressure fields, respectively.

- **Computing Flux**:
  - The flux is computed as the product of the velocity perturbation and the average of pressure perturbations at adjacent levels:

    \[ \text{Flux} = U_p \times \frac{P_{p, \text{level } n} + P_{p, \text{level } n+1}}{2} \]

- **Scaling with Grid Spacing and Cell Fraction**:
  - The flux is scaled by the vertical grid spacing (DRF) and the horizontal grid spacing (DXG or DYG), multiplied by the cell fraction (hFacW, hFacS, or hFacC):

    \[ \text{Scaled Flux} = \text{Flux} \times \text{DRF} \times (\text{DXG or DYG}) \times h \]

- **Aggregation Over the Grid**:
  - The final flux is obtained by aggregating these scaled fluxes, either by summing or averaging over certain dimensions:

    \[ \text{Total Flux} = \sum (\text{Scaled Flux}) \]
    \[ \text{Mean Flux} = \text{mean}(\text{Scaled Flux}) \]

