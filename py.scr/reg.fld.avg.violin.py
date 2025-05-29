import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as ds
from freq_funcs import tparse
from datetime import datetime,timedelta
import seaborn as sb

#idir  =  "/home/rangke/work/thesis/input/final/3havg/10sn"
idir  =  "/home/rangke/work/thesis/input/vari/10sn"
odir  =  "/home/rangke/work/thesis/figs"

# 3havg
osres          =  "qd"
models         =  ["arpege","gem","geos","grist","gsam","icon","icon","ifs","mpas","scream","shield","um" ,"e5","obs" ]
isress         =  ["2km"   ,"5km","3km" ,"5km"  ,"4km" ,"2km" ,"5km" ,"4km","3km" ,"3km"   ,"3km"   ,"5km","qd","qd"  ]
osress         =  ["qd"    ,"qd" ,"qd"  ,"qd"   ,"qd"  ,"qd"  ,"qd"  ,"qd" ,"qd"  ,"qd"    ,"qd"    ,"qd" ,"qd","qd"  ]

modfnm         =  ["ARPEGE 2km","GEM 5km","GEOS 3km","GRIST 5km","GSAM 4km","ICON 2km","ICON 5km","IFS 4km","MPAS 3km","SCREAM 3km","SHIELD 3km","UM 5km","ERA5 0.25\u00b0","Observation"]
clist          =  ["C0","C1","C2","C3","C4","C5","C6","C7","C8","C9","lightsteelblue","olive","teal","black"]
lwlist         =  [4,4,4,4,4,4,4,4,4,4,4,4,4,6]
rclist         =  [(0,0),(0,1),(1,0),(1,1),(2,0),(2,1)]


reglist        =  ["amz","afr","pco","ino","lnd","ocn"]
regs           =  len(reglist)
tres           =  "havg"


evnm           =  ["lhf","shf","swnt"]

vnm            =  "pw"     # pr, lhf, shf

fext           =  "png"

if (vnm == "pw"):
   pos            =  [1,2,3,4,5,6,7,8,9,10,11,12,13]
   entries        =  len(models) - 1
   zolist         =  np.arange(0,entries,1)

   print(entries,len(zolist),zolist)
else:
   pos            =  [1,2,3,4,5,6,7,8,9,10,11,12,13,14]
   entries        =  len(models)
   zolist         =  np.arange(0,entries,1)

   print(entries,len(zolist),zolist)

fig,axs        =  plt.subplots(3,2,figsize=(12,18))


if    (vnm  == "pr"):
   fvnm     =  "Precipitation rate"
   uts      =  "mm day\u207b\u00b9"
   vmax     =  36
   vmin     =  -1
elif  (vnm  == "lhf"):
   fvnm     =  "Latent heat flux"
   uts      =  "W m\u207b\u00b2"
   vmax     =  400
   vmin     =  50
elif  (vnm  == "shf"):
   fvnm     =  "Sensible heat flux"
   uts      =  "W m\u207b\u00b2"
   vmax     =  400
   vmin     =  50
elif  (vnm  == "swnt"):
   fvnm     =  "SWNT"
   uts      =  "W m\u207b\u00b2"
   vmax     =  400
   vmin     =  50
elif  (vnm  == "lwnt"):
   fvnm     =  "Outgoing longwave radiation"
   uts      =  "W m\u207b\u00b2"
   vmax     =  400
   vmin     =  50
elif  (vnm  == "pw"):
   fvnm     =  "Precipitable water"
   uts      =  "mm"
   vmax     =  90
   vmin     =  -1


print("  {0:s}".format(fvnm))

ofn         =  "reg.avg.{0:s}.{1:s}".format(vnm,fext)
ofp         =  "{0:s}/{1:s}".format(odir,ofn)




# min
# max
# mean
# std
# median
var_stats   =  np.empty((regs,entries,5))


for r in range(regs):
#for r in range(0,1):

   reg         =  reglist[r]

   col         =  rclist [r][0]
   row         =  rclist [r][1]

   mt_arr      =  []

   for m in range(entries):
   
      model    =  models[m]
      isres    =  isress[m]
      osres    =  osress[m]
      color    =  clist [m]
      zodr     =  zolist[m]
      lw       =  lwlist[m]
      mfnm     =  modfnm[m]

   
      if    (model == "obs" and vnm == "pr"):
         mfnm        =  "IMERG 0.1\u00b0"
      elif  (model == "obs" and vnm == "pw"):
         mfnm        =  "MIMIC TPW2"
      elif  (model == "obs" and vnm in ["lwnt","swnt"]):
         mfnm        =  "CERES 1.0\u00b0"
         isres       =  "od"
         osres       =  "od"
   
   
      if    (reg == "amz"):
         lsm_ifn  =  "lsm.amz0.{0:s}.nc".format(osres)
         tlbl     =  "Amazon biome"
         mskv     =  1
      elif  (reg == "pco"):
         lsm_ifn  =  "lsm.pco0.{0:s}.nc".format(osres)
         tlbl     =  "Pacific ocean"
         mskv     =  1
      elif  (reg == "afr"):
         lsm_ifn  =  "lsm.afr0.{0:s}.nc".format(osres)
         tlbl     =  "African continent"
         mskv     =  1
      elif  (reg == "ino"):
         lsm_ifn  =  "lsm.ino0.{0:s}.nc".format(osres)
         tlbl     =  "Indian ocean"
         mskv     =  1
      elif  (reg == "lnd"):
         lsm_ifn  =  "lsm.ocn1.lnd0.10sn.{0:s}.nc".format(osres)
         tlbl     =  "Land only (global)"
         mskv     =  0
      elif  (reg == "ocn"):
         lsm_ifn  =  "lsm.ocn1.lnd0.10sn.{0:s}.nc".format(osres)
         tlbl     =  "Ocean only (global)"
         mskv     =  1
   
   
      lsm_ifp     =  "{0:s}/{1:s}".format(idir,lsm_ifn)
      lsm_id      =  ds("{0:s}".format(lsm_ifp)) 
   
      lsm_id.set_auto_mask(False)
   
      ilsm        =  lsm_id["lsm"][:]
   
      print("  region : {0:20s} | model : {1:12s}".format(tlbl,mfnm))
   

      if (vnm in ["lhf","shf"] and model == "obs"):
         print("        skipping obs for heat flux")
         continue


      ifn            =  "{0:s}.{1:s}.{2:s}.10sn.{3:s}.{4:s}.nc".format(model,isres,vnm,osres,tres)
      ifp            =  "{0:s}/{1:s}".format(idir,ifn)

      fld_id         =  ds("{0:s}".format(ifp))

      ilat           =  fld_id["lat"][:]
      ilon           =  fld_id["lon"][:]

      if (model == "grist"):
         ivar           =  fld_id[vnm]  [24:]
         idx            =  np.where((ilsm[24:] == mskv) & (ivar >= 0))
         #idx            =  np.where(ilsm[24:] == mskv)
      else:
         ivar           =  fld_id[vnm]  [:]
         idx            =  np.where((ilsm == mskv) & (ivar >= 0))
         #idx            =  np.where(ilsm == mskv)


      sel_dat        =  ivar[idx].compressed()
      mt_arr.append(sel_dat)


      var_stats[r,m,0]   =  np.min   (sel_dat)
      var_stats[r,m,1]   =  np.max   (sel_dat)
      var_stats[r,m,2]   =  np.std   (sel_dat)
      var_stats[r,m,3]   =  np.mean  (sel_dat)
      var_stats[r,m,4]   =  np.median(sel_dat)

   vio_parts      =  axs[col,row].violinplot(mt_arr,positions=pos,showmeans=True,showmedians=True,widths=0.8)

   for comp, coll in vio_parts.items():
      if comp == "bodies":
         for pc, co in zip(coll, clist):
            pc.set_facecolor(co)
            pc.set_edgecolor(co)

      elif  (comp == "cmedians"):
         coll.set_color(clist)
         coll.set_linestyle(":")


      else:
         coll.set_color(clist)
      

   axs[col,row].set_xticks([],[])
   axs[col,row].grid()
   axs[col,row].set_title("{0:s}".format(tlbl))
   axs[col,row].set_ylabel("{0:s} ({1:s})".format(fvnm,uts))
   #axs[col,row].set_yscale("log")
   axs[col,row].set_ylim(vmin,vmax)

   lsm_id.close()
   fld_id.close()


if (vnm == "lwnt"):
  axs[2,0].set_ylim(-30,450)
elif  (vnm == "pw"):
  axs[2,0].set_ylim(-5,85)


for r in range(regs):
#for r in range(0,1):
   model_var   =  var_stats[r]
   reg         =  reglist[r]

   for m in range(entries):
      mfnm     =  modfnm   [m]
      min_val  =  model_var[m,0]
      max_val  =  model_var[m,1]
      std_val  =  model_var[m,2]
      avg_val  =  model_var[m,3]
      med_val  =  model_var[m,4]

      print(
            "{0:s} | {1:20s} | {2:15s} | ".format(reg,fvnm,mfnm),
            "max = {0:7.2f} {1:s} | ".format(max_val,uts),
            "min = {0:7.2f} {1:s} | ".format(min_val,uts),
            "std = {0:7.2f} {1:s} | ".format(std_val,uts),
            "avg = {0:7.2f} {1:s} | ".format(avg_val,uts),
            "med = {0:7.2f} {1:s} | ".format(med_val,uts),
            )

   print("")


for m in range(entries):
   mfnm     =  modfnm[m]
   co       =  clist [m]

   if    (mfnm == "Observation" and vnm == "pw"):
      mfnm  =  "MIMIC-TPW2"
   elif  (mfnm == "Observation" and vnm == "pr"):
      mfnm  =  "IMERG V06"
   elif  (mfnm == "Observation" and vnm == "lwnt"):
      mfnm  =  "CERES 1.0\u00b0"
      
   axs[2,1].scatter(zolist,np.ones(entries) * -100,label=mfnm,alpha=0.8,c=co)

axs[2,1].legend(fontsize=10,bbox_to_anchor=(0.9, -0.10),ncol=7)
plt.subplots_adjust(left=0.05,right=0.99,bottom=0.07,top=0.98,hspace=0.15,wspace=0.25)
plt.savefig("{0:s}".format(ofp))
#plt.show()
