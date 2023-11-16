
import xarray as xr
import xroms
from glob import glob


# Call this function to load a file 
def loadRun(filename, chunks, old=False, avg=False):
    if old:
        ds, grid = loadOldRun(filename, chunks)
    else:
        ds, grid = loadNewRun(filename, chunks, avg)
    return ds, grid


def loadOldRun(filename, chunks):
    ds = xr.open_dataset(filename, chunks=chunks)
    ds['zeta'] = 0*ds.hbl # didn't save zeta...oops.
    ds = ds.rename(time='ocean_time')
    ds, grid = xroms.roms_dataset(ds)
    ds['temp'] = (-(1000+ds.rho)/(1000+ds.rho.isel(ocean_time=0, s_rho=-1).mean()) + 1)/2e-4 + 20
    return ds, grid

def loadNewRun(filename, chunks, avg=False):
    # note here that filename should be the directory of history files
    if avg:
        files = sorted(glob(filename+'jet_avg.*.nc'))
    else:
        files = sorted(glob(filename+'jet_his.*.nc'))
#     dsWC = xr.open_mfdataset(files, combine='nested', concat_dim='time', data_vars='minimal', chunks=chunks, coords='minimal', compat='override', parallel=True)
    dsWC = xr.open_mfdataset(files, combine='nested', concat_dim='time', data_vars='minimal', chunks=chunks, coords='minimal', compat='override', parallel=True, decode_cf = False, engine='h5netcdf')
    
    ds = xr.open_dataset('../data/model/jet_his.nc')

    # For some reason I have to correct the grid metrics that aren't being output with this CROCO run.
    dsWC['s_rho'] = ('s_rho', ds.s_rho)
    dsWC['s_w'] = ('s_w', ds.s_w)
    dsWC['Cs_r'] = ds.Cs_r
    dsWC['Cs_w'] = ds.Cs_w
    dsWC['h'] = (('eta_rho', 'xi_rho'), 4000  + 0*dsWC.zeta.isel(time=0))
    dsWC['pm'] = (('eta_rho', 'xi_rho'), 1/500 + 0*dsWC.zeta.isel(time=0))
    dsWC['pn'] = (('eta_rho', 'xi_rho'), 1/500 + 0*dsWC.zeta.isel(time=0))
    dsWC['pm_u'] = (('eta_rho', 'xi_u'), 1/500 + 0*dsWC.u.isel(time=0, s_rho=0))
    dsWC['pn_v'] = (('eta_v', 'xi_rho'), 1/500 + 0*dsWC.v.isel(time=0, s_rho=0))

    #dsWC['f'] = dsNC.f ## XXX THIS IS A HACK THAT WON"T ALWAYS WORK
    #dsWC = dsWC.rename(time='ocean_time')
    print('here')
    dsWC, grid = xroms.roms_dataset(dsWC, Vtransform=1)
    print('and here')
    dsWC['rho'] = dsWC.temp
    dsWC = dsWC.swap_dims({'time':'ocean_time'})
    dsWC['temp'] = (-(1000+dsWC.rho)/(1000+dsWC.rho.isel(ocean_time=1, s_rho=-1).mean()) + 1)/2e-4 + 20
    return dsWC, grid


def loadZarrRun(filename, chunks='auto'):
    dsWC = xr.open_zarr(filename, chunks=chunks)
    ds = xr.open_dataset('../data/model/jet_his.nc')

    # For some reason I have to correct the grid metrics that aren't being output with this CROCO run.
    dsWC['s_rho'] = ('s_rho', ds.s_rho)
    dsWC['s_w'] = ('s_w', ds.s_w)
    dsWC['Cs_r'] = ds.Cs_r
    dsWC['Cs_w'] = ds.Cs_w
    dsWC['h'] = (('eta_rho', 'xi_rho'), 4000  + 0*dsWC.hbl.isel(time=0))
    dsWC['pm'] = (('eta_rho', 'xi_rho'), 1/500 + 0*dsWC.hbl.isel(time=0))
    dsWC['pn'] = (('eta_rho', 'xi_rho'), 1/500 + 0*dsWC.hbl.isel(time=0))
    dsWC['pm_u'] = (('eta_rho', 'xi_u'), 1/500 + 0*dsWC.u.isel(time=0, s_rho=0))
    dsWC['pn_v'] = (('eta_v', 'xi_rho'), 1/500 + 0*dsWC.v.isel(time=0, s_rho=0))

    #dsWC = dsWC.rename(time='ocean_time')
    dsWC, grid = xroms.roms_dataset(dsWC, Vtransform=1)

    dsWC['rho'] = dsWC.temp
    dsWC = dsWC.swap_dims({'time':'ocean_time'})
    dsWC['temp'] = (-(1000+dsWC.rho)/(1000+dsWC.rho.isel(ocean_time=1, s_rho=-1).mean()) + 1)/2e-4 + 20
    return dsWC, grid


def newLoad(filename, chunks='auto', avg=False):
    # note here that filename should be the directory of history files
    from glob import glob

    
    if avg:
        files = sorted(glob(filename+'jet_avg.*.nc'))
    else:
        files = sorted(glob(filename+'jet_his.*.nc'))
    dsWC = xr.open_mfdataset(files, combine='nested', concat_dim='time', data_vars='minimal', chunks=chunks, coords='minimal', compat='override', parallel=True)

    ds = xr.open_dataset('../data/model/jet_his.nc')

    dsWC, grid = xroms.roms_dataset(dsWC, Vtransform=1) ### IS THIS CORRECT???? NEED TO VERIFY!! Correct 8/26
    dsWC['rho'] = dsWC.temp
    dsWC = dsWC.swap_dims({'time':'ocean_time'})
    dsWC['temp'] = (-(1000+dsWC.rho)/(1000+dsWC.rho.isel(ocean_time=1, s_rho=-1).mean()) + 1)/2e-4 + 20
    return dsWC, grid