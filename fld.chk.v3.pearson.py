import numpy as np
from netCDF4 import Dataset as ds
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
import matplotlib.colors as colors
from scipy import stats 
from freq_funcs import *

# comparing different fields of dyamond2 model runs
# also comparing tropical domain mean and root mean square error

idir     =  "/home/rangke/work/thesis/input/vari/od"
odir     =  "/home/rangke/work/thesis/figs"
odir_tex =  "/home/rangke/doc_works/thesis"

model_list  =  ["arpege","gem","geos","grist","gsam","icon","icon","ifs","mpas","scream","shield","um","e5","obs"]
isres_list  =  ["2km"   ,"5km","3km" ,"5km"  ,"4km" ,"2km" ,"5km" ,"4km","3km" ,"3km"   ,"3km","5km","qd","od"]

obs_list    =  ["ce","mtpw","imrg"]
obs_isrlist =  ["od","qd","01d"]

var_list    =  ["lwnt","pw","pr"]
fvnm_list   =  ["Outgoing longwave radiation","Precipitable water","Precipitation rate"]
#uts_list    =  ["W/m\u00b2","mm","mm/day"]
uts_list    =  ["W m\u207b\u00b2","mm","mm day\u207b\u00b9"]
#idx_list    =  [(0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(2,0),(2,1),(2,2),(3,0),(3,1),(3,2),(4,1)] # 5 x 3 grid
idx_list    =  [(1,0),(1,1),(2,0),(2,1),(3,0),(3,1),(4,0),(4,1),(5,0),(5,1),(6,0),(6,1),(0,1),(0,0)] # 7 x 2 grid

maxv_list   =  [330,77 ,40]
cntv_list   =  [270,30 ,10 ]
minv_list   =  [120 ,0  ,-0.1]
tix_loc_arr =  [
               [80,130,180,220,270,320],
               [0,15,30,45,60,75],
               [0,10,20,30,40],
               ]
#               [0,20,40,60,80,100,120,140,160],
# box lat,lon coordinates for regional analysis
# lat : [low,high]
# lon : [low,high]
# % pacific ocean region needs two separate boxes
# region order : africa, amazon, indian ocean, pacific ocean 1, pacific ocean 2
# coordinate order : left top, right top, left low, right low

box_coords  =  [
               [(-10,5),(15,30)],
               [(-10,5),(-75,-50)],
               [(-10,10),(60,90)],
               [(-10,10),(-180,-120)],
               [(-10,10),(160,180)],
               ]

boxes       =  len(box_coords)



tfull_list  =  ["Feb. 2020"]
tidx_list   =  [0,1,2,3,4]

rad_vars    =  ["lwnt"]
pr_vars     =  ["pr"]
pw_vars     =  ["pw"]
noa_vars    =  ["lhf","shf","skt"]

dsize       =  len(model_list)
vsize       =  len(var_list)

exceptions  =  ["skt","lhf","shf"]


reg         =  "trp"
osres       =  "od"
tres        =  "havg"

ilat        =  np.arange(-30,30+1,1)
ilon        =  np.arange(-180,180,1)

lon,lat     =  np.meshgrid(ilon,ilat)


title_fs    =  8
tick_fs     =  8
unit_fs     =  8
color_map   =  "Spectral_r"
fext        =  "png"


clvl        =  40

margs       =  {
               "llcrnrlat" :np.min(ilat),
               "urcrnrlat" :np.max(ilat),
               "llcrnrlon" :np.min(ilon),
               "urcrnrlon" :np.max(ilon),
               "projection":"cyl",
               "fix_aspect":"False",
               }

ncol        =  2
nrow        =  7

obsn        =  "obs"
obs_isr     =  "od"


for v in range(vsize):
#for v in [0,2]:

   vnm         =  var_list   [v]
   fvnm        =  fvnm_list  [v]
   isres       =  isres_list [v]
   uts         =  uts_list   [v]
   maxv        =  maxv_list  [v]
   cntv        =  cntv_list  [v]
   minv        =  minv_list  [v]
   tix_loc     =  tix_loc_arr[v]
   tix_lbl     =  numarr2strarr(tix_loc)

   obsn        =  obs_list   [v]
   obs_isr     =  obs_isrlist[v]


   #print("minv : {0:5.2f} | vcnt : {1:5.2f} | vmax : {2:5.2f}".format(minv,cntv,maxv))

   clvls       =  np.linspace(minv,maxv,clvl)
   cmp_lvls    =  colors.TwoSlopeNorm(vmin=minv,vcenter=cntv,vmax=maxv)

   ofn         =  "fld.{0:s}.{1:s}.{2:s}".format(vnm,"feb.mavg",fext)
   ofp         =  "{0:s}/{1:s}".format(odir,ofn)
   ofp_tex     =  "{0:s}/{1:s}".format(odir_tex,ofn)

   print("output : {0:s}".format(ofp))
   print("output : {0:s}".format(ofp_tex))

   fig,axs     =  plt.subplots(nrow,ncol,figsize=(8,9))

   for m in range(dsize):

      model    =  model_list[m]
      isres    =  isres_list[m]

      nr       =  idx_list[m][0]
      nc       =  idx_list[m][1]


      if    (model == "obs" and vnm == "lwnt"):
         model    =  "ce"
         isres    =  "od"
         src      =  "CERES"
      elif  (model == "obs" and vnm == "pw"  ):
         model    =  "mtpw"
         isres    =  "qd"
         src      =  "MIMIC TPW2"
      elif  (model == "obs" and vnm == "pr"  ):
         model    =  "imrg"
         isres    =  "01d"
         src      =  "IMERG"
      elif  (model == "e5"):
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


      obs_ifn  =  "{0:s}.{1:s}.{2:s}.{3:s}.{4:s}.{5:s}.nc".format(obsn,obs_isr,vnm,reg,osres,tres)
      obs_ifp  =  "{0:s}/{1:s}".format(idir,obs_ifn)


      ifn      =  "{0:s}.{1:s}.{2:s}.{3:s}.{4:s}.{5:s}.nc".format(model,isres,vnm,reg,osres,tres)
      ifp      =  "{0:s}/{1:s}".format(idir,ifn)

      obs_id   =  ds("{0:s}".format(obs_ifp),"r",format="NETCDF4")
      id       =  ds("{0:s}".format(ifp    ),"r",format="NETCDF4")

      #obs_ivar =  obs_id[vnm][:]
      obs_ivar =  obs_id[vnm][288:]
      obs_avg  =  np.nanmean(obs_ivar,axis=0)

      ivar     =  id[vnm  ][288:]
      #ivar     =  id[vnm  ][:]
      ilat     =  id["lat"][:]
      ilon     =  id["lon"][:]

      ivar_avg =  np.nanmean(ivar,axis=0)
      ivar_mn  =  np.nanmean(ivar)


      # only considering non-masked values
      val_obs  =  np.where(ivar_avg.mask == False, obs_avg,np.nan)
      val_mod  =  np.where(ivar_avg.mask == False,ivar_avg,np.nan)

      val_obs_nn  =  val_obs.flatten()
      val_mod_nn  =  val_mod.flatten()

      pear     =  stats.pearsonr(val_obs_nn[~np.isnan(val_obs_nn)],val_mod_nn[~np.isnan(val_mod_nn)])


      # for scikit-learn module
#      mod_rmse =  mean_squared_error(val_obs[~np.isnan(val_obs)], val_mod[~np.isnan(val_mod)], squared=False)

      print("model : {0:6s} | mean : {1:5.2f} {2:s} | pear : stat = {3:5.2f} | pval = {4:5.2e}".format(model,ivar_mn,uts,pear[0],pear[1]))

      fig_title=  "{0:s} {1:s} (mean = {2:5.2f} {3:s}, corr. = {4:5.2f})".format(src.upper(),sres,ivar_mn,uts,pear[0])


      if (vnm == "pr"):
         ivar[ivar < 0.0] = 0.0


      axs[nr,nc].set_title(fig_title,fontsize=title_fs)
      bmap     =  Basemap(**margs,ax=axs[nr,nc])
      bmap.drawcoastlines()



      cntf     =  bmap.contourf(lon,lat,ivar_avg,clvl,levels=clvls,ax=axs[nr,nc],norm=cmp_lvls,latlon=True,cmap=color_map)

      # drawing boxes for regional averages
      for b in range(boxes):

         coords   =  box_coords[b]

         lat_hi   =  coords[0][1]
         lat_lo   =  coords[0][0]
         lon_hi   =  coords[1][1]
         lon_lo   =  coords[1][0]

         lat_coords  =  [lat_lo,lat_hi,lat_hi,lat_lo]
         lon_coords  =  [lon_lo,lon_lo,lon_hi,lon_hi]

         #x,y   =  bmap(lon_coords,lat_coords)
         #xy    =  zip(x,y)
         #poly  =  Polygon(list(xy), facecolor="black", alpha=1.0, fill=False, lw=1)
         #axs[nr,nc].add_patch(poly)

         draw_poly_sub(lat_coords,lon_coords,bmap,axs[nr,nc])


      id.close()
      obs_id.close()

   
   #print("rmse_{0:s} = {1:}".format(vnm,rmse_arr))
   #print("mean_{0:s} = {1:}".format(vnm,mean_arr))

   # colorbar specs
   #fig.add_axes([xloc, yloc, length, height])

   #cb = plt.colorbar(cntf,location = "bottom", ax = axs[0,1], ticks=tix_loc,cax = fig.add_axes([0.52, 0.93, 0.46, 0.020]))
   cb = plt.colorbar(cntf,location = "bottom", ax = axs[6,0], ticks=tix_loc,cax = fig.add_axes([0.01, 0.05, 0.98, 0.020]))
   cb.set_label("{0:s} ({1:s})".format(fvnm,uts),fontsize=8)
   cb.ax.tick_params(labelsize=8)
   cb.ax.set_xticklabels(tix_lbl)

   # 5 x 3 grid
   #axs[4,0].axis("off")
   #axs[4,2].axis("off")

   # 7 x 2 grid
   #axs[0,1].axis("off")

   #plt.subplots_adjust(left=0.005,right=0.995,bottom=0.06,top=0.999,hspace=0.005,wspace=0.01) # for 5 x 3 grid
   plt.subplots_adjust(left=0.005,right=0.995,bottom=0.07,top=0.999,hspace=0.1,wspace=0.01)



   
   plt.savefig("{0:s}".format(ofp))
#   plt.savefig("{0:s}".format(ofp_tex))
#   plt.show()

   plt.figure().clear()

