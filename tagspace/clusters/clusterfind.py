import os
import numpy as np
from galpy.util import multi as ml
from astropy.io import fits
from tagspace import tagdir
from tagspace.data import gettimestr
from tagspace.data.spectra import spectra

class tag(object):
	def __init__(self,data=None,starinfo=None,maxcores=1,save=True):
		if isinstance(data,str):
			if data[0] == '/' or data[0] == '~':
				hdulist = fits.open(data)
			else:
				hdulist = fits.open(tagdir+'/'+data)
			self.datainfo = hdulist[0].header
			self.starparams = hdulist[self.datainfo['param_index']].data
			if 'abundance_index' in self.datainfo.keys():
				self.abundances = hdulist[self.datainfo['abundance_index']].data
			if 'spectra_index' in self.datainfo.keys():
				self.spectra = spectra()
				self.spectra.specs = hdulist[self.datainfo['spectra_index']].data
		self.maxcores = maxcores
		self.save = save
		return None

	def cluster_wrapper(self,i):
		repeatlist = []
		for r in range(self.repeats):
			predict = self.clusterfn(**self.kwargs).fit_predict(self.clusterdata[i])
			repeatlist.append(predict)
		return repeatlist

	def cluster(self,datatype='abundances',clusterfn,repeats=1,**kwargs):
		if isinstance(datatype,str):
			if datatype == 'abundances'
				self.clusterdata = self.abundances
			elif datatype == 'spectra':
				self.clusterdata = self.spectra.specs
		self.clusterfn = clusterfn
		self.repeats = repeats
		self.kwargs = kwargs
		self.labels_pred = ml.parallel_map(self.cluster_wrapper,range(self.instances),
										   numcores=self.maxcores)
		if self.save:
			self.savecluster()

	def savecluster(self):
		directory = '{0}/{1}/'.format(self.directory,self.clusterfn.__name__)
		if not os.path.isdir(directory):
			os.system('mkdir -p {0}'.format(directory))
		name = '{0}_{1}repeats_{2}'.format(self.datainfo['name'],self.repeats,gettimestr())
		np.savez(directory+name, labels_true=np.array(self.labels_true), labels_pred=np.array(self.labels_pred))
