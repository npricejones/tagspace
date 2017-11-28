import numpy as np
import pandas as pd
from tagspace.wrappers.genfns import normalgeneration
from tagspace import tagdir,ptnumdict,ptsymdict
#from tagspace.data.spectra import spectra

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
			self.name=''
			self.centergenfn = genfn
			self.instances = instances
			self.numcluster = numcluster
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
			self.centerdata = pd.DataFrame()
			self.centernames = ['center-{0}'.format(i) for i in self.elemnames]
			entry = 0
			for i in range(self.instances):
				indx = [i]*self.numcluster
				clustercenters = self.centergenfn(num=self.numcluster,
												  numelem=self.numelem,**kwargs)
				abundf = pd.DataFrame(clustercenters,columns=self.centernames,index=np.arange(entry,entry+self.numcluster,dtype=int))
				abundf['instances'] = pd.Series(indx,index=abundf.index)
				abundf['labels_true'] = pd.Series(np.arange(self.numcluster,dtype=int),index=abundf.index)
				self.centerdata = self.centerdata.append(abundf)
				entry += self.numcluster
		return None


	def create_abundances(self,genfn=normalgeneration,nummembers=20,**kwargs):
		self.datatype = 'abundances'
		self.membergenfn = genfn
		self.nummembers = nummembers
		self.clusterdata = pd.DataFrame()
		if isinstance(self.nummembers,(int,float)):
			self.nummembers = np.array([self.nummembers]*self.numcluster)
		entry = 0
		for i in range(self.instances):
			centers = self.centerdata[self.centerdata['instances']==i]
			for c in range(len(centers)):
				# parallel spot
				indx = [i]*self.nummembers[c]
				label = [c]*self.nummembers[c]
				# passing means implies normal
				center = np.array(centers[self.centernames][centers['labels_true']==c])[0]
				clustermembers = self.membergenfn(num=self.nummembers[c],
												  numelem=self.numelem,
												  centers=center, **kwargs)
				abundf = pd.DataFrame(clustermembers,columns=self.elemnames,index=np.arange(entry,entry+self.nummembers[c],dtype=int))
				abundf['instances'] = pd.Series(indx,index=abundf.index)
				abundf['labels_true'] = pd.Series(label,index=abundf.index)
				self.clusterdata = self.clusterdata.append(abundf)
				entry += self.nummembers[c]
		return None

	def savecluster(self,name=None):
		directory = '{0}/{1}/{2}/{3}'.format(self.directory,
											 self.centergenfn.__name__,
											 self.membergenfn.__name__,
											 self.datatype)
		if not os.path.isdir(directory):
			os.system('mkdir -p {0}'.format(directory))
		name = '{0}{1}'.format(gettimestr(),self.name)
		self.instancedata.to_hdf(directory+'/'+name)
