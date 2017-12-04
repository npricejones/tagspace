import h5py
import copy
import numpy as np
import pandas as pd
from tagspace import tagdir
from tagspace.data import gettimestr
from tagspace.wrappers import getwrapperattrs
from tagspace.wrappers.genfns import normalgeneration
from tagspace import tagdir,ptnumdict,ptsymdict
#from tagspace.data.spectra import spectra

class makeclusters(object):
	"""
	"""
	def __init__(self,genfn=normalgeneration,instances=1,readdata=False,
				 filename=None,numcluster=20,numelem=10,maxcores=1,
				 elems=np.array([6,7,8,11,12,13,14,16,19,20]),
				 **kwargs):
		"""
		"""
		if readdata:
			if filename[0] == '/' or filename[0] == '~':
				self.datafile = h5py.File(filename,'r')
			else:
				self.datafile = h5py.File(tagdir+'/'+filename,'r')
		elif not readdata:
			self.datafile = h5py.File(tagdir+'/synthetic_clusters.hdf5','w')
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
			self.centerpaths = []
			self.centerdata= np.zeros((self.instances,self.numcluster,self.numelem))
			kwargdict = copy.deepcopy(kwargs)
			for i in range(self.instances):
				currenttime = gettimestr()
				instance = self.datafile.create_group(currenttime+'/'+self.centergenfn.__name__)
				self.centerpaths.append(instance)
				clustercenters = self.centergenfn(num=self.numcluster,
												  numelem=self.numelem,**kwargs)
				self.centerdata[i] = clustercenters
				kwargdict.update({'num':self.numcluster,'numelem':self.numelem})
				getwrapperattrs(instance,self.centergenfn,kwargdict=kwargdict)
				instance.attrs['centers'] = clustercenters
				instance.attrs['elemnames'] = self.elemnames
				instance.attrs['elemnums'] = self.elems
		return None


	def create_abundances(self,genfn=normalgeneration,nummembers=20,readdata=False,path=None,**kwargs):
		"""
		"""
		self.datatype = 'abundances'
		if readdata:
			pass
			# call function from __init__ to get info
			# path could specify some date info
		elif not readdata:
			self.membergenfn = genfn
			self.nummembers = nummembers
			self.members = np.zeros((self.instances,np.sum(self.nummembers),self.numelem))
			self.labels_true = -np.ones(self.instances,np.sum(self.nummembers))
			kwargdict = copy.deepcopy(kwargs)
			if isinstance(self.nummembers,(int,float)):
				self.nummembers = np.array([self.nummembers]*self.numcluster)
			for i in range(self.instances):
				centers = self.centerdata[i]
				instance = self.centerpaths[i].create_group(self.membergenfn.__name__)
				members = np.zeros((np.sum(self.nummembers),self.numelem))
				labels_true = -np.ones(np.sum(self.nummembers))
				starpos = 0
				for c in range(len(centers)):
					# parallel spot
					label = [c]*self.nummembers[c]
					clustermembers = self.membergenfn(num=self.nummembers[c],
													  numelem=self.numelem,
													  centers=centers, **kwargs)
					members[starpos:starpos+self.nummembers[c]] = clustermembers
					labels_true[starpos:starpos+self.nummembers[c]] = label
					starpos+= self.nummembers[c]
				self.members[i] = members
				self.labels_true[i] = labels_true
				instance['clustermembers'] = members
				memberinfo = instance['clustermembers']
				memberinfo.attrs['labels_true'] = labels_true
				memberinfo.attrs['datatype'] = self.datatype
				memberinfo.attrs['elemnames'] = self.elemnames
				memberinfo.attrs['elemnums'] = self.elems
				kwargdict.update({'num':self.nummembers,'numelem':self.numelem})
				getwrapperattrs(instance,self.membergenfn,kwargdict=kwargdict)
		return None
