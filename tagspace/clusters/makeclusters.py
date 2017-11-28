import numpy as np
import pandas as pd
from tagspace import tagdir,ptnumdict,ptsymdict
from tagspace.data.spectra import spectra

class makeclusters(object):
	def __init__(self,genfn=normalgeneration,instances=1,readdata=False,
				 fname=None,numcluster=20,numelem=10,
				 elems=np.array([6,7,8,11,12,13,14,16,19,20]),
				 **kwargs):
		if readdata:
			if filename[0] == '/' or filename[0] == '~':
				self.clusterdata = pd.read_hdf(filename)
			else:
				self.clusterdata = pd.read_hdf(tagdir+'/'+filename)
		elif not readdata:
			self.instances = instances
			self.numcluster = numcluster
			if isinstance(self.numcluster,(int,float)):
				self.numcluster = np.array([self.numcluster]*self.instances)
			self.elems = elems
			self.elemnames = np.zeros(self.elems.shape,dtype='S2')
			for e in range(len(elems)):
				try:
					self.elems[e] = int(elems[e])
					self.elemnames[e] = ptsymdict[elems[e]]
				except ValueError:
					self.elemnames[e] = elems[e]
					self.elems[e] = ptnumdict[elems[e]]
			self.numelem = len(elems)
			self.clusterdata = pd.DataFrame()
			for i in range(self.instances):
				indx = [i]*self.numbercluster[i]
				dfdict = {'labels_true':pd.Series(np.arange(self.numcluster[i],dtype=int),
												  index=indx)}
				clustercenters = self.centergenfn(numcluster=self.numcluster[i],
												  numelem=self.numelem,**kwargs)
				for abun in range(clustercenters.shape[1]):
					dfdict[self.elemnames[abun]] = pd.Series(clustercenters[:,abun],index=indx)
				tempdf = pd.DataFrame(dfdict)
				self.clusterdata.append(tempdf)
		return None


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
