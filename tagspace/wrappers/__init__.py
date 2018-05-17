import inspect
import warnings
import h5py

def getwrapperattrs(group,callable,kwargdict={'num':10}):
   """
   Adds attributes to HDF5 group according to parameters 
   """
   properties = inspect.getargspec(callable)
   allargs = properties[0]
   defaults = properties[3]
   arglen = len(allargs) - len(defaults)
   if arglen > 0:
      warnings.warn('This function has non-keyword arguments and not all reproduction information will be stored')
   for i in range(len(defaults)):
      print('adding attribute {0}'.format(allargs[i+arglen]))
      if allargs[i+arglen] in kwargdict.keys():
         try:
            group.attrs.create('{0}'.format(allargs[i+arglen]), kwargdict[allargs[i+arglen]],dtype=h5py.special_dtype(vlen=str))
         except RuntimeError:
            pass
      elif allargs[i+arglen] not in kwargdict.keys():
         try:
            group.attrs.create('{0}'.format(allargs[i+arglen]), defaults[i],dtype=h5py.special_dtype(vlen=str))
         except RuntimeError:
            pass
   return None
