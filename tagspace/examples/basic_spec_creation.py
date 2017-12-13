import h5py
import numpy as np
from tagspace.clusters.makeclusters import makeclusters
from tagspace.wrappers.genfns import normalgeneration
from tagspace.data.spectra import psmspectra

clusters = makeclusters()

clusters.create_abundances()

mem = np.sum(clusters.nummembers)

teffdict = {'num':mem,'numprop':1,'centers':4500*np.ones(mem),'stds':100*np.ones(mem)}
loggdict = {'num':mem,'numprop':1,'centers':2.6*np.ones(mem),'stds':0.1*np.ones(mem)}
vturbdict = {'num':mem,'numprop':1,'centers':6*np.ones(mem),'stds':np.ones(mem)}
clusters.get_photosphere(nummembers=mem,params=[['TEFF'],['LOGG'],['VTURB']],
						 genfns=[normalgeneration,normalgeneration,
						 		 normalgeneration],
						 kwargdicts=[teffdict,loggdict,vturbdict])
datafile = h5py.File(clusters.synfilename,'r+')
datakey = datafile['normalgeneration'].keys()[-1]
abundances = datafile['normalgeneration/'+datakey]
specinfo = psmspectra(mem,clusters.photosphere)