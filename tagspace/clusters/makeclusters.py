import os
import h5py
import copy
import inspect
import numpy as np
import pandas as pd
from tagspace import tagdir
from tagspace.data import gettimestr
from tagspace.data.spectra import psmspectra
from tagspace.wrappers import getwrapperattrs
from tagspace import tagdir,ptnumdict,ptsymdict
from numpy.lib.recfunctions import append_fields
from tagspace.wrappers.genfns import normalgeneration,choice2dgeneration,lineargeneration
#from tagspace.data.spectra import spectra

class makeclusters(object):
	"""
	Used to generate synthetic cluster abundances and spectra.
	"""
	def __init__(self,genfn=normalgeneration,instances=1,readdata=False,
				 filename=None,numcluster=20,maxcores=1,
				 elems=np.array([6,7,8,11,12,13,14,16,19,20]),
				 **kwargs):
		"""
		Reads cluster center abundances from file or generates cluster 
		center abundances

		genfn:		Function used to find cluster center abundances 
					(defaults to choosing from a normal distribution)
		instances:	The number of clusters to create with the same 
					parameters (defaults to one)
		readdata:	Boolean that looks for existing data if set True 
					(redundant with filename******)
		filename:	If readdata, read cluster center information from this 
					path. Unless the path starts are root '/' or home '~', 
					assumes data lives in the TAGSPACE_DATA environment 
					variable
		numcluster:	number of clusters to generate (defaults to 20)
		maxcores:	Maxmimum number of cores to use for parallel processes
		elems:		List or array of atomic numbers or symbols for elements
					to be generated (defaults to: carbon, nitrogen, oxygen,
					sodum, magnesium, aluminum, silicon, sulfur, potassium,
					and calcium)
		**kwargs:	Passed to genfn

		Returns None
		"""
		self.maxcores = maxcores
		# Read data and assign it to relevant class attributes (the latter
		# not yet implemented)
		if readdata:
			
			if filename[0] == '/' or filename[0] == '~':
				self.filename = filename
			
			else:
				self.filename = tagdir+'/'+filename
			
			datafile = h5py.File(self.filename,'r')

		# Generate data
		elif not readdata:
			# Assign basic class attributes
			self.centergenfn = genfn
			self.synfilename = tagdir+'synthetic_clusters/centergen_'+self.centergenfn.__name__
			if not os.path.isdir(self.synfilename):
				os.system('mkdir -p {0}'.format(self.synfilename))
			self.synfilename += '/clustering_data.hdf5'
			if os.path.isfile(self.synfilename):
				self.datafile = h5py.File(self.synfilename,'r+')
			elif not os.path.isfile(self.synfilename):
				self.datafile = h5py.File(self.synfilename,'w')
			self.instances = instances
			self.numcluster = numcluster
			self.elems = elems
			self.elemnames = np.zeros(self.elems.shape,dtype='S2')

			# Confirm that form in which the list of elements was received,
			# either atomic numbers or symbols
			for e in range(len(elems)):
				try:
					self.elems[e] = int(elems[e])
					self.elemnames[e] = ptsymdict[elems[e]]
				except ValueError:
					self.elemnames[e] = elems[e]
					self.elems[e] = ptnumdict[elems[e]]
			self.numelem = len(elems)

			# Data will be uniquely stored by formation timestamp, store as an 
			# array for convenience
			self.timestamps = np.zeros(self.instances,dtype='S100')

			# Cluster center data to be stored in one large array
			self.centerdata= np.zeros((self.instances,self.numcluster,self.numelem))
			# Update kwarg dictionary with parameters that are class attributes 
			kwargs.update({'num':self.numcluster,'numprop':self.numelem})

			# Create cluster center abundances for each requested instance
			for i in range(self.instances):
				# Get unique timestamp for this instance
				currenttime = gettimestr()
				self.timestamps[i] = currenttime
				# Generate cluster center abundances
				clustercenters = self.centergenfn(**kwargs)
				self.centerdata[i] = clustercenters

				# Assign attributes to HDF5 dataset
				dsetname = 'center_abundances_{0}'.format(currenttime)
				self.datafile[dsetname] = clustercenters
				centerinfo = self.datafile[dsetname]
				getwrapperattrs(centerinfo,self.centergenfn,kwargdict=kwargs)
				centerinfo.attrs['elemnames'] = self.elemnames
				centerinfo.attrs['atmnums'] = self.elems
			self.datafile.close()
		return None


	def create_abundances(self,genfn=normalgeneration,nummembers=20,
						  readdata=False,path=None,**kwargs):
		"""
		Given a generation function and its kwargs, create cluster members.

		genfn:		Function used to find cluster member abundances 
					(defaults to choosing from a normal distribution)
		nummembers:	Number of stars per cluster; integer or iterable of 
					integers with length of makeclusters.numclusters
		readdata:	Boolean that looks for existing data if set True 
					(redundant with filename******)
		path:		Some information about where to get files ********
		**kwargs:	Passed to genfn

		Returns None
		"""
		if readdata:
			pass
			# call function from __init__ to get info
			# path could specify some date info
		elif not readdata:
			# Assign more class attributes
			self.datafile = h5py.File(self.synfilename,'r+')
			self.membergenfn = genfn
			self.nummembers = nummembers

			# Make nummembers into an interable
			if isinstance(self.nummembers,(int)):
				self.nummembers = np.array([self.nummembers]*self.numcluster)

			# Create holder arrays
			self.abundances = np.zeros((self.instances,
										np.sum(self.nummembers),
										self.numelem))
			self.labels_true = -np.ones((self.instances,
										 np.sum(self.nummembers)))
			
			# Deep copy kwargs so it can be updated with user specified kwargs
			kwargdict = copy.deepcopy(kwargs)
			for i in range(self.instances):

				# Get cluster centers for this instance
				centers = self.centerdata[i]

				# Create group to store data 
				dsetpath = self.membergenfn.__name__
				if dsetpath not in self.datafile:
					instance = self.datafile.create_group(dsetpath)
				elif dsetpath in self.datafile:
					instance = self.datafile[dsetpath]

				# Create array to store members and true labels
				self.members = np.zeros((np.sum(self.nummembers),self.numelem))
				labels_true = -np.ones(np.sum(self.nummembers))
				
				# Create members for each cluster
				starpos = 0
				for c in range(len(centers)):
					# parallel spot
					label = [c]*self.nummembers[c]
					clustermembers = self.membergenfn(num=self.nummembers[c],
													  numprop=self.numelem,
													  centers=centers[c], **kwargs)
					self.members[starpos:starpos+self.nummembers[c]] = clustermembers
					labels_true[starpos:starpos+self.nummembers[c]] = label
					starpos+= self.nummembers[c]
				
				# Update class properties
				self.abundances[i] = self.members
				self.labels_true[i] = labels_true

				# Create data set in file
				instance['member_abundances_{0}'.format(self.timestamps[i])] = self.members
				memberinfo = instance['member_abundances_{0}'.format(self.timestamps[i])]
				# Assign attributes
				memberinfo.attrs['datatype'] = 'abundances'
				memberinfo.attrs['labels_true'] = labels_true
				memberinfo.attrs['elemnames'] = self.elemnames
				memberinfo.attrs['atmnums'] = self.elems
				kwargdict.update({'num':self.nummembers,'numprop':self.numelem})
				getwrapperattrs(memberinfo,self.membergenfn,kwargdict=kwargdict)
			self.datafile.close()
		return None

	def get_photosphere(self,readdata=False,path=None,nummembers=20,
					    params=[['TEFF','LOGG'],['VTURB']],
					    genfns=[choice2dgeneration,lineargeneration],
						kwargdicts=[{'num':10,'numprop':2,
								  'readdata':'{0}/APOGEE/rc_teff_logg.npy'},
								 {'num':10,'numprop':1,'indeps':'LOGG',
								  'slope':0.32,
								  'intercept':2.48+0.1*np.random.randn(20)}]):
		if readdata:
			pass
		# Make nummembers into an interable
		if isinstance(self.nummembers,(int)):
			self.nummembers = np.array([self.nummembers]*self.numcluster)
		paramlabels =  [x.upper() for sublist in params for x in sublist]
		# If self.photosphere doesn't exist, initlialize it
		self.photosphere = getattr(self,'photosphere',
								   np.zeros(np.sum(self.nummembers),
								   			dtype=[(paramlabels[0],float)]))
		
		for param in range(len(params)):
			genfn = genfns[param]
			genresults = genfn(**kwargdicts[param])
			for p in range(len(params[param])):
				# Check for possiblity that this attribute is being updated
				try:
					#setattr(self,param[p],genresults[:,p])
					self.photosphere[params[param][p]] = genresults[:,p]
				except ValueError:
					self.photosphere = append_fields(self.photosphere,
													 params[param][p],
													 genresults[:,p])
		return None


	def create_spectra(self,readdata=False,path=None,nummembers = 20,
					   params=[['TEFF','LOGG'],['VTURB']],
					   genfns=[choice2dgeneration,lineargeneration],
					   kwargdicts=[{'num':10,'numprop':2,
								'readdata':'{0}/APOGEE/rc_teff_logg.npy'},
							   {'num':10,'numprop':1,'indeps':'LOGG',
								'slope':0.32,
								'intercept':2.48+0.1*np.random.randn(20)}],
					   specclass=psmspectra,specfn='member',**kwargs):
		"""
		Given a generation function and its kwargs, create cluster members.

		genfn:		Function used to find cluster member spectra
					(defaults to choosing from a normal distribution)
		nummembers:	Number of stars per cluster; integer or iterable of 
					integers with length of makeclusters.numclusters
		readdata:	Boolean that looks for existing data if set True 
					(redundant with filename******)
		path:		Some information about where to get files ********
		**kwargs:	Passed to specclass.specfn

		Returns None
		"""
		if readdata:
			pass
		elif not readdata:
			# split up kwargs for functions
			genfnposkwargs = inspect.getargspec(genfn)[0]
			specgenposkwargs = inspect.getargspec(specgen)[0]
			genfnkwargs = {}
			specgenkwargs = {}
			for kwarg in kwargs.keys():
				if kwarg in genfnposkwargs:
					genfnkwargs[kwarg] = kwargs[kwarg]
				elif kwarg in specgenposkwargs:
					specgenkwargs[kwarg] = kwargs[kwarg]

			self.datafile = h5py.File(self.synfilename,'w')
			self.membergenfn = genfn
			self.nummembers = nummembers	

			# Make nummembers into an interable
			if isinstance(self.nummembers,(int)):
				self.nummembers = np.array([self.nummembers]*self.numcluster)

			self.get_photosphere(nummembers=self.nummembers,params=params,
								 genfns=genfns,kwargdicts=kwargdicts)

			for i in range(self.instances):
				self.sinfo = specclass(self.nummembers,self.photosphere,
									   self.abundances)
				if specfn == 'member':
					self.sinfo.from_member_abundances(**kwargs)

				elif specfn == 'center':
					self.sinfo.from_center_abundances(**kwargs)

				elif specfn == 'spectra':
					self.sinfo.from_center_spectrum(**kwargs)

