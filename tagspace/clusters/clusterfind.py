import os
import numpy as np
from galpy.util import multi as ml
from astropy.io import fits
from tagspace import tagdir
from tagspace.data import gettimestr
from tagspace.data.spectra import spectra
from tagspace.clusters.makeclusters import makeclusters

class tag(makeclusters):
	"""
	Using existing data, find clusters using a user-specified cluster finding algorithm
	structured after the template of scikit-learn (e.g. a class with a fit_predict 
	function)

	"""
	def __init__(self,centergenfn=normalgeneration,
				 starinfo=None,maxcores=1,datatype='abundances',
				 instances=1,save=True,readdata=False,fname=None,
				 **kwargs):
		"""
		Accepts a set of stellar information, including temperatures and surface gravities
		as well as spectra or chemical abundances, in addition to a list of assignment indices
		that label each star with the number of the cluster to which it belongs.
		"""
		makeclusters.__init__(genfn=normalgeneration,instances=instances,
							  save=save,maxcores=maxcores,readdata=readdata,
							  fname=fname,**kwargs)
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

	def savetag(self):
		directory = '{0}/{1}/'.format(self.directory,self.clusterfn.__name__)
		if not os.path.isdir(directory):
			os.system('mkdir -p {0}'.format(directory))
		name = '{0}_{1}repeats_{2}'.format(self.datainfo['name'],self.repeats,gettimestr())
		np.savez(directory+name, labels_true=np.array(self.labels_true), labels_pred=np.array(self.labels_pred))
