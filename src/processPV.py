# To RUN on command line:
# Deactivate all envs
# source activate /homes/metofac/wenegrat/.conda/envs/CFB_EBF_3
# check > which python
# > python processPV.py (may need to do twice)
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from dask.distributed import Client
import dask
import xarray as xr
import xgcm
from dask.diagnostics import ProgressBar
import sys
import matplotlib.patches as patches
sys.path.append('/homes/metofac/wenegrat/xroms/')
#sys.path.append('/homes/metofac/wenegrat/pyspec/')

#from xroms import open_roms_netcdf_dataset
import cmocean.cm as cmo
import xroms
#from pyspec import spectrum as spec
from scipy import integrate as integrate
from timeit import default_timer as timer

#%%
plt.rcParams['text.usetex'] = True
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['font.size'] = 14
plt.rcParams['contour.negative_linestyle'] = 'solid'

import sys
sys.path.append("../src/")
import load_model_runs as lruns

from importlib import reload  
import time


client = Client(processes=False)
client


## LOAD THE NEW RUNS
reload(lruns)
#dsWC, gridWC = lruns.loadRun('../data/model/JET_ML_CFB/', chunks={'time':1}, old=False, avg=False)
dsWC, gridWC = lruns.newLoad('../data/model/JET_NML_CFB_BFLUX/', chunks={'time':1})
runname='JET_NML_CFB_BFLUX_NEW'
dsWC

dsNC, gridNC = lruns.newLoad('../data/model/JET_NML_NOCFB_BFLUX/', chunks={'time':1})


def calcIntegratedPV(dst, ts, yl, zl=slice(43,100)):
    dst = dst.swap_dims({'ocean_time':'time'})
    dst['buoyancy'] = xroms.buoyancy(dst.rho, rho0=dst.rho0)
    t = time.time()
    PV = dst.isel(time=ts).xroms.ertel
    #PV_areaweigh = PV*dst.dz0*dst.dA 
    PV_areaweigh = PV*(dst.dz.isel(time=ts))*dst.dA
    PV_sum = PV_areaweigh.isel(eta_rho=yl, s_rho=zl).sum().compute()
    return PV_sum
    
yl = slice(1000, 4002-1000)#FZ
#yl = slice(900, 1100) #WZ
#yl = slice(1900, 2100) #SZ
nt = 31 
print('Processing WC')
PV_WC = np.zeros((nt,))
for i in range(0,nt):
    tic = time.perf_counter()
    PV_WC[i] = calcIntegratedPV(dsWC, i, yl)
    toc = time.perf_counter()
    print(f'Processing step {i}/{nt} took: {(toc-tic)/60} min', end="\r") # Should take approximately 1.2 min/step

print('Processing NC')
PV_NC = np.zeros((nt,))
for i in range(0,nt):
    tic = time.perf_counter()
    PV_NC[i] = calcIntegratedPV(dsNC, i, yl)
    toc = time.perf_counter()
    print(f'Processing step {i}/{nt} took: {(toc-tic)/60} min', end="\r") # Should take approximately 1.2 min/step
    
#save to interim     
dspv = xr.Dataset(
        data_vars = dict(
        pv_wc=(['time'], PV_WC),
        pv_nc=(['time'], PV_NC),    
        ),
        coords = dict(
        time=dsWC.time.values[0:nt]
        )
)
dspv.attrs['type'] = 'noflux'
dspv.attrs['yl_l'] = yl.start
dspv.attrs['yl_r'] = yl.stop
#dspv.attrs['zl_l'] = zl.start
dspv.to_netcdf('../data/interim/pvintegrated_500depth_FZ_bflux_Zt.nc')    
