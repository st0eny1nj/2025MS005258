import numpy as np
from netCDF4 import Dataset as ds


def olr2tb(in_olr):
   sig   =  5.67e-8
   a     =  1.228
   b     =  -1.106e-3
   tf =  (in_olr / sig)**(0.25)

   out_tb1   =  (-a + np.sqrt(a**2 + 4 * b * tf)) / (2 * b)
#   out_tb2   =  (-a - np.sqrt(a**2 + 4 * b * tf)) / (2 * b)

   return out_tb1
#   return out_tb2


ireg        =  "20sn"

#idir     =  "/home/rangke/work/thesis/input/3havg/tsrt"
#odir     =  "/home/rangke/work/thesis/input/3havg/tsrt"

idir     =  "/home/rangke/work/thesis/input/vari/{0:s}".format(ireg)
odir     =  "/home/rangke/work/thesis/input/vari/{0:s}".format(ireg)


osres       =  "qd"
tres        =  "havg"
ivnm        =  "lwnt"
ovnm        =  "tb"


models      =  ["arpege","geos","gem","grist","gsam","icon","icon","ifs","mpas","scream","shield","um","e5"]
isress      =  ["2km","3km","5km","5km","4km","2km","5km","4km","3km","3km","3km","5km","qd"]

#models      =  ["e5"]
#isress      =  ["qd"]

model_siz   =  len(models)

for m in range(model_siz):

   model    =  models[m]
   isres    =  isress[m]

   ifp      =  "{0:s}/{1:s}.{2:s}.{3:s}.{4:s}.{5:s}.{6:s}.nc".format(idir,model,isres,ivnm,ireg,osres,tres)
   ofp      =  "{0:s}/{1:s}.{2:s}.{3:s}.{4:s}.{5:s}.{6:s}.nc".format(odir,model,isres,ovnm,ireg,osres,tres)
   
   print("")
   print("lwnt in : {0:s}".format(ifp))
   print("tb  out : {0:s}".format(ofp))
   print("")

   id       =  ds("{0:s}".format(ifp),"r",format="NETCDF4")
   od       =  ds("{0:s}".format(ofp),"w",format="NETCDF4")
   
   
   iolr     =  id["lwnt"]  [:]
   ilat     =  id["lat"]   [:]
   ilon     =  id["lon"]   [:]
   itime    =  id["time"]
    
   
   otb      =  olr2tb(iolr)
   
   dtime    =  od.createDimension("time", None     )
   dlat     =  od.createDimension("lat" , len(ilat))
   dlon     =  od.createDimension("lon" , len(ilon))
   
   vtime    =  od.createVariable("time","f4",("time")            )
   vlat     =  od.createVariable("lat" ,"f4",("lat")             )
   vlon     =  od.createVariable("lon" ,"f4",("lon")             )
   vtb      =  od.createVariable("tb"  ,"f4",("time","lat","lon"),fill_value=-9999.0,compression="zlib")
   
   
   vtime.standard_name  =  "time"
   vtime.units          =  itime.units
   
   vlat.standard_name   =  "lat"
   vlat.units           =  "degree_north"
   
   vlon.standard_name   =  "lon"
   vlon.units           =  "degree_east"
   
   vtb.standard_name    =  "brightness temperature"
   vtb.units            =  "K"
   
   vtime[:]             =  itime[:]
   vlat[:]              =  ilat
   vlon[:]              =  ilon
   vtb[:]               =  otb
   
   id.close()
   od.close()
