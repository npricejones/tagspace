import inspect
import warnings

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
		if allargs[i+arglen] in kwargdict.keys():
			group.attrs['{0}'.format(allargs[i+arglen])] = kwargdict[allargs[i+arglen]]
		elif allargs[i+arglen] not in kwargdict.keys():
			group.attrs['{0}'.format(allargs[i+arglen])] = defaults[i]
	return None