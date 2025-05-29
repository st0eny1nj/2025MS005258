import numpy as np
from netCDF4 import Dataset as ds
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as colors
from freq_funcs import *
#from sklearn.metrics import mean_squared_error

# comparing different fields of dyamond2 model runs
# also comparing tropical domain mean and root mean square error

idir     =  "/home/rangke/work/thesis/input/vari/od"
odir     =  "/home/rangke/work/thesis/figs"
odir_tex =  "/home/rangke/doc_works/thesis"

model_list  =  ["arpege","gem","geos","grist","gsam","icon","icon","ifs","mpas","scream","shield","um","e5"]
isres_list  =  ["2km"   ,"5km","3km" ,"5km"  ,"4km" ,"2km" ,"5km" ,"4km","3km" ,"3km"   ,"3km","5km","qd"]
color_list  =  ["C0","C1","C2","C3","C4","C5","C6","C7","C8","C9","lightsteelblue","olive","teal"]


#model_list  =  ["grist","e5"]
#isres_list  =  ["5km","qd"]
#color_list  =  ["m","cyan"]


obs_list    =  ["ce","imrg","e5"]
obs_isrlist =  ["od","01d","qd"]

var_list    =  ["lwnt","pr","pw"]
fvnm_list   =  ["Outgoing longwave radiation","Precipitation rate","Precipitable water"]
uts_list    =  ["W m\u207b\u00b2","mm day\u207b\u00b9","mm"]


rad_vars    =  ["lwnt"]
pr_vars     =  ["pr"]
pw_vars     =  ["pw"]
noa_vars    =  ["lhf","shf","skt"]

dsize       =  len(model_list)
vsize       =  len(var_list)


msk_ocn     =  "n"
msk_land    =  "y"

reg         =  "trp"
osres       =  "od"
tres        =  "havg"

fext        =  "pdf"

obsn        =  "obs"
obs_isr     =  "od"

lwid        =  2.0

if    (msk_ocn == "y" or msk_land == "y"):
   if    (msk_ocn == "y"):
      lsm_ifn  =  "lsm.od.trp.ocn1.havg.nc"
      lsm_ifp  =  "{0:s}/{1:s}".format(idir,lsm_ifn)
      tlbl     =  "Land only"
      opt      =  "land.only"
      print("  land only")
   elif  (msk_land== "y"):
      lsm_ifn  =  "lsm.od.trp.land1.havg.nc"
      lsm_ifp  =  "{0:s}/{1:s}".format(idir,lsm_ifn)
      tlbl     =  "Ocean only"
      opt      =  "ocn.only"
      print("  ocean only")

   lsm_id      =  ds("{0:s}".format(lsm_ifp),"r",format="NETCDF4")
   ilsm        =  lsm_id["lsm"][288:]
   ilsm_mskd   =  np.ma.masked_equal(ilsm,1)

else:

   tlbl        =  "Ocean and land"
   opt         =  "ocn.and.land"
   print("  including ocean and land...")



#for v in range(vsize):
for v in range(0,2):
   
   vnm         =  var_list   [v]
   fvnm        =  fvnm_list  [v]
   uts         =  uts_list   [v]

   obsn        =  obs_list   [v]
   obs_isr     =  obs_isrlist[v]


   ofn         =  "zonal.avg.{0:s}.{1:s}.{2:s}".format(vnm,opt,fext)
   ofp         =  "{0:s}/{1:s}".format(odir,ofn)
   ofp_tex     =  "{0:s}/{1:s}".format(odir_tex,ofn)

   print("output : {0:s}\n".format(ofp))
#   print("output : {0:s}".format(ofp_tex))

   if    (vnm == "lwnt"):
      obsn     =  "ce"
      isres    =  "od"
      src      =  "CERES"
      sres     =  "1.0\u00b0"
   elif  (vnm == "pw"  ):
      obsn     =  "e5"
      isres    =  "qd"
      src      =  "ERA5"
      sres     =  "0.25\u00b0"
   elif  (vnm == "pr"  ):
      obsn     =  "imrg"
      isres    =  "01d"
      src      =  "IMERG V06"
      sres     =  "0.1\u00b0"


   fig,axs     =  plt.subplots(figsize=(13,6))


   if (vnm == "lhf"):
      print("  skippig LHF OBS...")
   else:
      obs_ifn  =  "{0:s}.{1:s}.{2:s}.{3:s}.{4:s}.{5:s}.nc".format(obsn,obs_isr,vnm,reg,osres,tres)
      obs_ifp  =  "{0:s}/{1:s}".format(idir,obs_ifn)
   
      obs_lbl  =  "{0:s} {1:s}".format(src,sres)
   
      obs_id   =  ds("{0:s}".format(obs_ifp),"r",format="NETCDF4")
      obs_ivar =  obs_id[vnm][288:]
      obs_ilat =  obs_id["lat"][:]

      if    (msk_ocn == "y" or msk_land == "y"):
         obs_ivar =  np.where(ilsm_mskd.mask == False, obs_ivar, np.nan)

      # obs monthly average
      obs_mavg =  np.nanmean(obs_ivar,axis=0)

      # obs zonal average
      obs_zoavg=  np.nanmean(obs_mavg,axis=1)
   
      # obs zonal average plot
      axs.plot(obs_ilat,obs_zoavg,label=obs_lbl,c="black",lw=lwid+1,zorder=vsize)

      obs_id.close()

   for m in range(dsize):

      model    =  model_list[m]
      isres    =  isres_list[m]
      lcolor   =  color_list[m]



      if (model == "e5"):
         src      =  "era5"
      else:
         src      =  model

      if    (isres == "qd"):
         sres  =  "0.25\u00b0"
      elif  (isres == "od" and vnm in rad_vars):
         sres  =  "1.0\u00b0"
      elif  (isres == "od" and vnm in pr_vars):
         sres  =  "0.1\u00b0"
      elif  (isres == "od" and vnm in pw_vars):
         sres  =  "0.25\u00b0"
      elif  (isres == "od" and vnm in noa_vars):
         sres  =  "0.25\u00b0"
      elif  (isres == "01d" and vnm in pr_vars):
         sres  =  "0.1\u00b0"
      else:
         sres  =  isres

      print("  doing {0:8s} {1:5s}".format(src,sres))
 
      ifn      =  "{0:s}.{1:s}.{2:s}.{3:s}.{4:s}.{5:s}.nc".format(model,isres,vnm,reg,osres,tres)
      ifp      =  "{0:s}/{1:s}".format(idir,ifn)

      m_lbl    =  "{0:s} {1:s}".format(src.upper(),sres)

      id       =  ds("{0:s}".format(ifp    ),"r",format="NETCDF4")

      ivar     =  id[vnm  ][288:]
      #ivar     =  id[vnm  ][:]
      ilat     =  id["lat"][:]
      ilon     =  id["lon"][:]

      if    (msk_ocn == "y" or msk_land == "y"):
         ivar  =  np.where(ilsm_mskd.mask == False, ivar, np.nan)
         ivar[ivar == -9999.0] = np.nan

         # NOTE!!!
         # Note that once you apply the np.where() method, NaN in the ivar will become fill value (-9999.0) again !
         # So, you need to change the -9999.0 value back into NaN in order for the data to be properly masked !! 

      # model monthly average
      mavg     =  np.nanmean(ivar,axis=0)

      # model zonal average
      zoavg    =  np.nanmean(mavg,axis=1)

      # model zonal average plot
      axs.plot(ilat,    zoavg,label=m_lbl,c=lcolor,lw=lwid)

      id.close()

   axs.set_title("{0:s}".format(tlbl))
   axs.set_ylabel("{0:s} ({1:s})".format(fvnm,uts))
   axs.set_xlabel("Latitude (degrees)")

#   axs.set_ylim(0,60)
   axs.yaxis.get_ticklocs(minor=True)
   axs.minorticks_on()

   plt.legend(bbox_to_anchor=(1.16, 1.02),fancybox=True, shadow=True)
   plt.subplots_adjust(left=0.06,right=0.87,top=0.95,bottom=0.10)
   plt.savefig("{0:s}".format(ofp))
#   plt.show() 

   plt.close()
