import numpy as np
import pandas as pd
import periodictable as pt
from tagspace import tagdir
from tagspace.data.spectra import spectra

class makeclusters(object):
	def __init__(self,genfn=normalgeneration,instances=1,readdata=False,
				 fname=None,numcluster=20,numelem=10,
				 elems=np.array([6,7,8,11,12,13,14,16,19,20]),
				 **kwargs):
		if readdata:
			if filename[0] == '/' or filename[0] == '~':
				hdulist = pd.read_hdf(filename)
			else:
				hdulist = pd.read_hdf(tagdir+'/'+filename)
		elif not readdata:
			self.instances = instances
			self.numcluster = numcluster
			if isinstance(self.numcluster,(int,float)):
				self.numcluster = np.array([self.numcluster]*self.instances)
			self.elems = elems
			self.names = np.zeros(self.elems.shape,dtype='S2')
			for e in range(len(elems)):
				try:
					atmnum = int(elems[e])
					self.names =
				except ValueError:
					self.name[e] = elems[e]
					self.elems[e] = pt.elements[elems[e]].number
			self.numelem = len(elems)
			for i in range(self.instances):
				clustercenters = self.centergenfn(numcluster=self.numcluster[i],
														  numelem=self.numelem,**kwargs)
				clusterdata = pd.DataFrame(self.clustercenters[i],columns = )
		return None

	def _element_list(self):


	def create_abundances(self,genfn=normalgeneration, **kwargs):
		self.datatype = 'abundances'
		self.membergenfn = genfn
		return None

	def savecluster(self,name=None):
		directory = '{0}/{1}/{2}/{3}'.format(self.directory,
											 self.centergenfn.__name__,
											 self.membergenfn.__name__,
											 self.datatype)
		if not os.path.isdir(directory):
			os.system('mkdir -p {0}'.format(directory))
		name = '{0}_{1}'.format(self.name,gettimestr())
		self.instancedata.to_hdf(directory+'/'+name)
