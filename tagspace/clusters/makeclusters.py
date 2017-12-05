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
		self.maxcores = maxcores
		if readdata:
			if filename[0] == '/' or filename[0] == '~':
				self.filename = filename
			else:
				self.filename = tagdir+'/'+filename
			datafile = h5py.File(self.filename,'r')
		elif not readdata:
			self.filename = tagdir+'/synthetic_clusters.hdf5'
			self.datafile = h5py.File(self.filename,'w')
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
			self.timestamps = []
			self.centerdata= np.zeros((self.instances,self.numcluster,self.numelem))
			instance = self.datafile.create_group(self.centergenfn.__name__)
			kwargdict = copy.deepcopy(kwargs)
			for i in range(self.instances):
				currenttime = gettimestr()
				self.timestamps.append(currenttime)
				clustercenters = self.centergenfn(num=self.numcluster,
												  numelem=self.numelem,**kwargs)
				self.centerdata[i] = clustercenters
				instance['center_abundances_{0}'.format(currenttime)] = clustercenters
				centerinfo = instance['center_abundances_{0}'.format(currenttime)]
				kwargdict.update({'num':self.numcluster,'numelem':self.numelem})
				getwrapperattrs(centerinfo,self.centergenfn,kwargdict=kwargdict)
				centerinfo.attrs['elemnames'] = self.elemnames
				centerinfo.attrs['elemnums'] = self.elems
			self.datafile.close()
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
			self.datafile = h5py.File(self.filename,'w')
			self.membergenfn = genfn
			self.nummembers = nummembers
			if isinstance(self.nummembers,(int,float)):
				self.nummembers = np.array([self.nummembers]*self.numcluster)
			self.abundances = np.zeros((self.instances,np.sum(self.nummembers),self.numelem))
			self.labels_true = -np.ones((self.instances,np.sum(self.nummembers)))
			kwargdict = copy.deepcopy(kwargs)
			for i in range(self.instances):
				centers = self.centerdata[i]
				instance = self.datafile.create_group(self.centergenfn.__name__+'/'+self.membergenfn.__name__)
				self.members = np.zeros((np.sum(self.nummembers),self.numelem))
				labels_true = -np.ones(np.sum(self.nummembers))
				starpos = 0
				for c in range(len(centers)):
					# parallel spot
					label = [c]*self.nummembers[c]
					clustermembers = self.membergenfn(num=self.nummembers[c],
													  numelem=self.numelem,
													  centers=centers, **kwargs)
					self.members[starpos:starpos+self.nummembers[c]] = clustermembers
					labels_true[starpos:starpos+self.nummembers[c]] = label
					starpos+= self.nummembers[c]
				self.abundances[i] = self.members
				self.labels_true[i] = labels_true
				instance['member_abundances_{0}'.format(self.timestamps[i])] = self.members
				memberinfo = instance['member_abundances_{0}'.format(self.timestamps[i])]
				memberinfo.attrs['labels_true'] = labels_true
				memberinfo.attrs['elemnames'] = self.elemnames
				memberinfo.attrs['elemnums'] = self.elems
				kwargdict.update({'num':self.nummembers,'numelem':self.numelem})
				getwrapperattrs(memberinfo,self.membergenfn,kwargdict=kwargdict)
			self.datafile.close()
		return None
