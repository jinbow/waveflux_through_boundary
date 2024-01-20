"""
Calculate the internal wave flux through a line in MITgcm.
The algorithm follows Nash et al. (2005)

jmd95.densjmd95 was used to calcualte density from the model field.
Download jmd95 following this link http://mitgcm.org/download/daily_snapshot/MITgcm/utils/python/MITgcmutils/MITgcmutils/jmd95.py 

"""

import popy
import numpy as np
from scipy import signal
import xarray as xr
xrod=xr.open_dataset

def saveh5(filename, var_name, data, attrs=None):
    """
    Append or update a variable in a NetCDF file.

    Parameters:
    - filename: str
        The path to the NetCDF file.
    - var_name: str
        The name of the variable to append or update.
    - data: array-like
        The data to be stored in the variable.
    - attrs: dict, optional
        A dictionary containing attributes to be added to the variable.

    Notes:
    - If the variable already exists in the file, it will be updated with the new data.
    - If the variable does not exist, it will be created and added to the file.
    
    Example:
    saveh5('example.nc', 'new_var', np.random.rand(10, 10), {'units': 'meters'})
    """
    # Open the existing NetCDF file
    with xr.open_dataset(filename, mode='a') as ds:
        # Assign the data to the variable
        ds[var_name] = data

        # Add attributes to the variable if provided
        if attrs is not None and isinstance(attrs, dict):
            ds[var_name].attrs.update(attrs)

        # Save the changes back to the file
        ds.to_netcdf(filename, mode='a')


def loadh5(filename, variable=None):
    """
    Load data from an HDF5 file using xarray.

    Parameters:
    - filename: str
        The path to the HDF5 file.
    - variable: str, optional
        The specific variable to load from the HDF5 file. If None, the entire file is loaded.

    Returns:
    - data: xarray.Dataset or xarray.DataArray
        The loaded data as an xarray Dataset (if variable is None) or DataArray (if variable is specified).

    Raises:
    - FileNotFoundError if the file does not exist.
    - KeyError if the specified variable is not found in the file.
    """
    # Open the HDF5 file using xarray
    ds = xr.open_dataset(filename, engine='h5netcdf')

    # If a specific variable is requested, return that variable
    if variable is not None:
        if variable in ds:
            return ds[variable]
        else:
            raise KeyError(f"Variable '{variable}' not found in file '{filename}'.")

    # Otherwise, return the entire dataset
    return ds


def c_pden():
    """
    Calculates and saves potential density data for specific locations.

    This function iterates over two locations: 'sth' and 'wst'. For each location,
    it loads salinity (S) and temperature (T) data from a specified .mat file,
    calculates the potential density using the JMD95 equation of state, and then
    saves the computed potential densities in a new .h5 file.

    The function operates on a file named 'BoundaryPropertiesForJinbo.mat' located
    in the '/nobackupp2/jwang23/projects/regional.mitgcm/' directory. It replaces
    potential density values with zero where salinity is zero to handle missing data.

    The resulting potential densities are saved in a file named with the location
    identifier and '_pden.h5' suffix, replacing the original '.mat' suffix.

    Note:
        - This function does not return any value.
        - It relies on the 'popy.io' and 'popy.jmd95' modules for data loading,
          potential density calculation, and data saving.
    """
    import jmd95
    dout={}
    for loc in ['sth', 'wst']:
        fn = '/nobackupp2/jwang23/projects/regional.mitgcm/BoundaryPropertiesForJinbo.mat'
        ss = xrod(fn)[f'S_{loc}']
        tt = xrod(fn)[f'T_{loc}']

        pden = jmd95.densjmd95(ss, tt, 0)
        pden[ss == 0] = 0
        fnout=fn.replace('.mat', '_pden.h5')
        
        dout[f'pden_{loc}']=pden.copy()
        
    xr.Dataset(dout).to_netcdf(fnout)

    return

def butter_highpass(data, cutoff, fs, order=5, axis=0):
    """
    Apply a high-pass Butterworth filter to the given data.

    This function designs a high-pass filter using the Butterworth design and
    applies it to the given data using a forward-backward filter method. The
    Butterworth filter is characterized by its maximally flat frequency response 
    in the passband.

    Parameters:
    - data: array_like
        The input data to be filtered.
    - cutoff: float
        The cutoff frequency of the filter in Hz. Frequencies below this value will be attenuated.
    - fs: float
        The sampling frequency of the input data in Hz.
    - order: int, optional (default is 5)
        The order of the filter. A higher order results in a sharper transition at the cutoff frequency.
    - axis: int, optional (default is 0)
        The axis along which the filter is applied in the input data array.

    Returns:
    - y: ndarray
        The filtered output with the same shape as the input data.

    Notes:
    - This function requires the `signal` module from the SciPy library.
    - The filter design is stable and avoids phase distortion by applying the filter forwards and backwards.
    """
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='highpass', analog=False)
    y = signal.filtfilt(b, a, data, axis=axis)

    return y

def c_perturbation(dd, hfac):
    """
    Calculate the perturbation field of a 3D variable.

    This function computes the perturbation field by:
    1. Removing the temporal linear trend from the data.
    2. Applying a high-pass filter to remove mesoscale features.

    The function also takes into account the vertical structure of the numerical grid by 
    weighting the perturbation calculation with the vertical distribution of 
    the variable (as provided by hfac).

    Parameters:
    - dd: ndarray
        A 3D numpy array representing the variable data over time. The expected
        dimensions are time, vertical level, and horizontal location.
    - hfac: ndarray
        A 3D numpy array representing the vertical fraction of each grid cell
        that is open water. The dimensions should match the vertical and 
        horizontal dimensions of dd.

    Returns:
    - dpp: ndarray
        The 3D array of the calculated perturbation field, with the same shape as dd.

    Note:
    - The function assumes the time dimension is the first axis in the input array.
    - Zero values in dd are treated as masked values and are excluded from the calculation.
    - The function uses the 'signal.detrend' method from the SciPy library for detrending
      and a custom 'butter_highpass' function for high-pass filtering.
    """
    msk = dd == 0
    dp = signal.detrend(dd, axis=0)
    totalH = hfac.sum(axis=0)[np.newaxis, :]  # Add time axis [nt, nx]
    totalH[totalH == 0] = np.inf

    dpm = (dp * hfac[np.newaxis, ...]).sum(axis=1) / totalH

    dp = dp - dpm[:, np.newaxis, :]
    dpp = np.apply_along_axis(butter_highpass, 0, dp, 1., 24.)
    dpp[msk] = 0
    del dp, dpm

    return dpp

class perturbation:
    def __init__(self, data, cutoff, fs, hfac, drf):
        """
        Initialize the perturbation class.

        This class is designed to compute the perturbation field of a 3D variable by
        applying a high-pass filter and removing the mean in the vertical dimension.

        Parameters:
        - data: ndarray
            A 3D numpy array (time, vertical level, horizontal location) representing the variable data.
        - cutoff: float
            The cutoff frequency for the high-pass filter in Hz.
        - fs: float
            The sampling frequency of the data in Hz.
        - hfac: ndarray
            A 3D numpy array (vertical level, horizontal y, horizontal x) representing the vertical fraction of each grid cell that is open water.
        - drf: ndarray
            A 1D numpy array representing the vertical thickness of each grid cell.
        """
        self.data = data
        self.cutoff = cutoff
        self.fs = fs
        self.hfac = hfac
        self.drf = drf
        return

    def highpass(self):
        """
        Apply a high-pass filter to the data.

        This method applies a high-pass filter to the data attribute of the class.
        The filter specifications are based on the cutoff frequency and the sampling
        frequency set during the class initialization.
        """
        self.d_highpass = butter_highpass(self.data, self.cutoff, self.fs)
        return

    def remove_zmean(self):
        """
        Remove the mean in the vertical dimension and apply a high-pass filter.

        This method first computes the weighted mean in the vertical dimension,
        using the vertical fraction of each grid cell (hfac) and their thickness (drf).
        It then subtracts this mean from the data to obtain the perturbation field.
        Finally, it applies a high-pass filter to this perturbation field.
        """
        msk = self.data == 0
        hfac = self.hfac * self.drf[:, np.newaxis, np.newaxis]  # final array has dimension (nz, ny, nx)
        totalH = hfac.sum(axis=0)[np.newaxis, np.newaxis, ...]  # add time axis to make totalH [nt, nz, ny, nx]
        totalH[totalH == 0] = np.inf

        d_zmean = (self.data * hfac[np.newaxis, ...]).sum(axis=1)[:, np.newaxis, :, :] / totalH  # dpm dimension (nt, 1, ny, nx)
        perturb = self.data - d_zmean
        dhighpass = butter_highpass(perturb, 1., 24.)
        self.perturb = dhighpass
        return


def calc_pressure():
    """
    Calculate and save the pressure field for different locations.

    This function calculates the pressure field based on potential density data and other
    related variables. It reads potential density data and other variables such as sea surface height
    and cell fraction data from specific HDF5 files and computes the pressure field using these inputs.

    The function operates on files named with '_pden.h5' and '.h5' suffixes derived from a base filename.
    It calculates the pressure for each location specified in the potential density data and saves 
    the pressure data back to the base HDF5 file with a specific variable name format.

    Steps:
    1. Load the 'DRF' variable from the base file and reshape it.
    2. Load potential density and other relevant data from the HDF5 files.
    3. Calculate the masked pressure field by integrating the product of potential density, 
       sea surface height, and cell fraction data.
    4. Save the calculated pressure field in the base HDF5 file.

    Notes:
    - The function assumes the existence of specific variables in the HDF5 files, such as 'H_', 'hFacC_', etc.
    - The pressure field is calculated for each location found in the potential density data.
    - Zero values in potential density data are treated as masked values and are excluded from calculations.
    """
    fn1 = fn.replace('.mat', '_pden.h5')
    fn0 = fn.replace('.mat', '.h5')

    drf = loadh5(fn, 'DRF').ravel()[np.newaxis, :, np.newaxis, np.newaxis]
    d1 = loadh5(fn1)
    
    for var in d1.keys():
        dd = loadh5(fn) 

        loc = var.split('_')[-1]
        ssh = dd['H_%s' % loc][:][:, np.newaxis, :, :]
        hfac = dd['hFacC_%s' % loc][:]

        d = np.ma.masked_equal(d1[var][:], 0)

        pden = d * (drf * hfac[np.newaxis, ...])
        pden[:, 0, :, :] = pden[:, 0, ...] + d[:, 0, ...] * ssh[:, 0, ...]
        
        pressure = 9.8 * pden.cumsum(axis=1)
        pressure[d.mask] = 0
        saveh5(fn,f'P_{loc}',pressure.astype('>f4'), {'name':'pressure','unit':'Pa'})
        del dd, d1
    return


def calc_allperturb():
    """
    Calculate and save the perturbation fields for various variables.

    This function processes variables from a specified .mat file and calculates their perturbation fields.
    It applies a high-pass filter to each variable to isolate fluctuations with periods shorter than 24 hours,
    and then removes the mean in the vertical dimension to obtain the perturbation field.

    The function identifies the variables based on their names and shapes, applies the appropriate cell fraction
    data (hFacC, hFacW, hFacS) depending on where they are on the c-grid, and then computes the perturbation
    using the 'perturbation' class.

    The resulting perturbation fields are saved in a new file with a '_perturb.h5' suffix.

    - The function expects the .mat file to be located at a specific path. In our case for CCS: west, north and south
    - The variable names should follow a specific naming convention (e.g., 'S_', 'T_', 'P_', 'U_', 'V_').
    - The shape of the variables in the file is expected to be 4-dimensional.
    - The function assumes that a zero value in the data indicates a masked or invalid value as in MITgcm convension
    - It relies on the 'perturbation' class to compute the perturbation fields.
    """
    fn = '/nobackupp2/jwang23/projects/regional.mitgcm/BoundaryPropertiesForJinbo.mat'
    d = loadh5(fn)
    varns = d.keys()

    drf = loadh5(fn, 'DRF').ravel()

    for var in varns:
        locs = var.split('_')
        if len(locs) == 2 and len(d[var].shape) == 4:
            loc = locs[-1]
            dd = d[var][:]
            if var[:2] in ['S_', 'T_', 'P_']:
                hfac = d['hFacC_%s' % loc][:]
            elif 'U_' in var:
                hfac = d['hFacW_%s' % loc][:]
            elif 'V_' in var:
                hfac = d['hFacS_%s' % loc][:]

            # High-pass filter, period < 24 hours
            aa = perturbation(dd, 1., 24., hfac.squeeze(), drf.squeeze())
            aa.remove_zmean()
            perturb = aa.perturb
            perturb[dd == 0] = 0
            del d, aa

            saveh5(fn.replace('.mat', '_perturb.h5'), var + '_p', perturb.astype('>f4'))
            d = loadh5(fn)

    return

def c_flux():
    """
    Calculate and save flux data for specified boundary.

    This function processes perturbation data for velocity (U or V) and pressure (P) components
    from a specified HDF5 file ('_perturb.h5'). It calculates the flux for each component by
    multiplying the velocity perturbation with the average pressure perturbation at adjacent
    levels. The calculation is performed for specified locations ('sth', 'nth', 'wst').

    The resulting flux data is saved back into the same HDF5 file with a specific naming convention.

    - The function expects the source HDF5 files to follow a specific naming convention and contain
      the necessary variables (perturbation of U, V, P and hFacS, hFacC, hFacW for each location).
    - The function assumes that the cell fraction data (hFac*) is used to mask or adjust the flux calculation.
    """
    fnn = fn.replace('.mat', '_perturb.h5')

    for loc in ['sth', 'nth', 'wst']:
        # Load velocity and pressure perturbation data
        if loc in ['sth', 'nth']:
            uvp = loadh5(fnn, 'V_%s_p' % loc)
            hfac1 = loadh5(fn, 'hFacS_%s' % loc)
        else:  # 'wst'
            uvp = loadh5(fnn, 'U_%s_p' % loc)
            hfac1 = loadh5(fn, 'hFacW_%s' % loc)

        uvp = uvp * hfac1[np.newaxis, ...]
        pp = loadh5(fnn, 'P_%s_p' % loc)
        hfac = loadh5(fn, 'hFacC_%s' % loc)
        pp = pp * hfac[np.newaxis, ...]

        # Calculate flux
        if loc in ['sth', 'nth']:
            aa = uvp[:, :, 0, :] * (pp[:, :, 0, :] + pp[:, :, 1, :]) / 2.0
            aa = aa * (hfac1 == 1)[np.newaxis, :, :, 0]
        else:  # 'wst'
            aa = uvp[:, :, :, 0] * (pp[:, :, :, 0] + pp[:, :, :, 1]) / 2.0
            aa = aa * (hfac1 == 1)[np.newaxis, :, :, 0]

        # Save the calculated flux
        saveh5(fnn, 'flux_%s' % loc, aa.astype('>f4'))

    return

def save_flux4Matt():
    """
    Save averaged flux data to a MATLAB file.

    This function reads the flux data for different locations ('wst', 'nth', 'sth') from
    a specified HDF5 file ('_perturb.h5') and calculates the mean of these flux datasets
    across the first axis. It then saves these averaged flux data into a MATLAB file with
    a '_flux.mat' suffix.

    The resulting MATLAB file contains three variables: 'flux_wst', 'flux_nth', and 'flux_sth',
    corresponding to the averaged flux data for the respective locations.

    - The function expects the HDF5 file to follow a specific naming convention and contain
      the required flux data variables.
    - The MATLAB file is saved in the same directory as the HDF5 file.

    - This function requires the 'scipy.io' module for saving data in MATLAB format.
    - It also relies on the 'popy.io' module for loading data from HDF5 files.
    """
    from scipy.io import savemat

    fnn = fn.replace('.mat', '_perturb.h5')
    dd = popy.io.loadh5(fnn)
    
    savemat(fn.replace('.mat', '_flux.mat'), {
        'flux_wst': dd['flux_wst'][:].mean(axis=0),
        'flux_nth': dd['flux_nth'][:].mean(axis=0),
        'flux_sth': dd['flux_sth'][:].mean(axis=0)
    })
    
    return 

def p_flux():
    """
    Plot and save figures
    """
    from matplotlib import gridspec
    fnn=fn.replace('.mat','_perturb.h5')
    dd=popy.io.loadh5(fn)
    rc=dd['RC'][0,:]
    drf=dd['DRF'][0,:]
    gs=gridspec.GridSpec(1,3)
    
#################################################
    flux_wst=popy.io.loadh5(fnn,'flux_wst').mean(axis=0)
    y=dd['YC_wst'][:,0]
    dy=dd['DYG_wst'][:,0]
    hfac=dd['hFacW_wst'][:,:,0]
    flux_wstfac=flux_wst*drf[:,np.newaxis]*hfac
    print(flux_wstfac.sum(axis=0).std()/1e3)
    print((flux_wstfac*dy[np.newaxis,:]).sum())
    print((flux_wstfac*dy[np.newaxis,:]).sum(axis=0))
    ax=plt.subplot(gs[0])
    plt.contourf(y,rc, flux_wst*10,levels=np.linspace(-1,1,11),extend='both',cmap=plt.cm.bwr) 
    plt.colorbar()
    plt.xlabel('latitude')
    plt.ylim(-5000,0)
#################################################
    flux_nth=popy.io.loadh5(fnn,'flux_nth').mean(axis=0)
    x=dd['XG_nth'][0,:]
    dx=dd['DXG_nth'][0,:]
    hfac=dd['hFacS_nth'][:,0,:]
    flux_nthfac=flux_nth*drf[:,np.newaxis]*dx[np.newaxis,:]*hfac
    print(flux_nthfac.sum())
    ax=plt.subplot(gs[1])
    plt.contourf(x,rc, -flux_nth*10,levels=np.linspace(-1,1,11),extend='both',cmap=plt.cm.bwr) 
    plt.colorbar()
    plt.xlabel('longitude')
    plt.xlim(x[x>0].min(),x[x>0].max())
    plt.ylim(-5000,0)

#################################################
    flux_sth=popy.io.loadh5(fnn,'flux_sth').mean(axis=0)
    x=dd['XC_sth'][0,:]
    dx=dd['DXG_sth'][0,:]
    hfac=dd['hFacS_sth'][:,0,:]
    flux_sthfac=flux_sth*drf[:,np.newaxis]*dx[np.newaxis,:]*hfac
    print(flux_sthfac.sum())

    ax=plt.subplot(gs[2])
    plt.contourf(x,rc, flux_sth*10,levels=np.linspace(-1,1,11),extend='both',cmap=plt.cm.bwr) 
    plt.colorbar()
    plt.xlabel('longitude')
    plt.ylim(-5000,0)


#################################################
    plt.tight_layout()
    plt.savefig('figures/boundary_fluxes.png',dpi=300)
    plt.show()

    return 



if __name__=="__main__":
    import pylab as plt

    fn='/nobackupp2/jwang23/projects/regional.mitgcm/BoundaryPropertiesForJinbo.mat'

    #c_pden()
    #calc_pressure()
    #calc_allperturb()
    #c_flux()
    p_flux()
    #save_flux4Matt()
