import os
import h5py
import copy
import inspect
import numpy as np
import pandas as pd
from tagspace import tagdir
from tagspace.data import gettimestr
#from tagspace.data.spectra import psmspectra
from tagspace.wrappers import getwrapperattrs
from tagspace import tagdir,ptnumdict,ptsymdict
from numpy.lib.recfunctions import append_fields
from tagspace.wrappers.genfns import normalgeneration,choice2dgeneration,lineargeneration
#from tagspace.data.spectra import spectra

class makeclusters(object):
	"""
	Used to generate synthetic cluster abundances and spectra.
	"""
	def __init__(self):
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
		return None
		# Read data and assign it to relevant class attributes (the latter
		# not yet implemented)
		'''
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
			self.synfilename = tagdir+'/synthetic_clusters/centergen_'+self.centergenfn.__name__
			if not os.path.isdir(self.synfilename):
				os.system('mkdir -p {0}'.format(self.synfilename))
			self.synfilename += '/clustering_data.hdf5'
			if os.path.isfile(self.synfilename):
				self.datafile = h5py.File(self.synfilename,'r+')
			elif not os.path.isfile(self.synfilename):
				self.datafile = h5py.File(self.synfilename,'w')
			self.datafile.close()
		return None
		'''

	def readh5py(self):
		return None

	def get_stellarparam(self,readdata=False,path=None,numstars=20,
					     params=[['TEFF','LOGG'],['VTURB']],
					     genfns=[choice2dgeneration,lineargeneration],
						 kwargdicts=[{'num':10,'numprop':2,
						   		     'readdata':'{0}/APOGEE/rc_teff_logg.npy'},
								     {'num':10,'numprop':1,'indeps':'LOGG',
								      'slope':0.32,
								      'intercept':2.48+0.1*np.random.randn(20)}]):
		"""
		Assign miscellaneous float attributes to star (e.g., effective 
		temperature, surface gravity, microturbulent velocity)

		readdata:		not implemented yet
		path:			not implemented yet
		numstars:		number of stars that need properties
		params:			list of parameters to assign. If parameters depend on
						each other, they should be grouped in sublists. 
						parameters independent of the others get their own
						sublist.
		genfns:			list of functions to generate the parameters.
		kwargdicts:		list of kwarg dictionaries to pass to generation
						functions

		Creates new attributes for the object, returns None
		"""
		if readdata:
			pass
		# Make nummembers into an interable
		self.numstars = numstars
		if isinstance(self.numstars,(int)):
			self.numstars = np.array([self.numstars]*self.numcluster)
		paramlabels =  [x.upper() for sublist in params for x in sublist]
		
		for param in range(len(params)):
			genfn = genfns[param]
			genresults = genfn(**kwargdicts[param])
			for p in range(len(params[param])):
				setattr(self,param[p],genresults[:,p])
		return None

	def polyindeps(self,degree=2,indeplabels=['TEFF']):
		"""
		Creates a matrix of independent values for a polynomial function.

		degree:			Degree of the polynomial
		indeplabels:	Names of parameters to use.

		Creates the matrix as the .indeps attributes, returns None.
		"""
		self.polynomial = PolynomialFeatures(degree=degree)
		indeps = np.zeros((self.numstars,len(indeplabels)))
		for i in range(len(indeplabels)):
			indeps[:,i] = getattr(self,indeplabels[i].lower())
			indeps[:,i] -= np.median(indeps[:,i])
		self.indeps = np.matrix(self.polynomial.fit_transform(indeps))
		return None

	def colpolyfit(self,i):
		"""
		Make a polynomial that is the best fit along each column of the
		data stored in .members.

		i:		index of the column to fit.

		Returns the best fit at column i of attribute .members.
		"""
		stars = np.matrix(self.members[:,i])
		covI = np.diag(1./self.uncertainties[:,i]**2)
		indeps = np.dot(self.indeps.T,np.dot(covI,self.indeps))
		stars = np.dot(self.indeps.T,np.dot(covI,stars.T))
		indepsI = np.linalg.inv(indeps)
		coeffs = np.dot(indepsI,stars)
		coeff_errs = np.array([np.sqrt(np.array(indepsI)[i][i]) for i in range(indeps.shape[1])])
		#if self.savecoeff:
		#	pass # save output somewhere
		return np.array(self.indeps*coeffs).flatten()

	def subpolyfit(self,maxcores=1,givencoeff=None,indeplabels=['TEFF'],
				   crossinds='all',degree=2,uncertainties=1,savecoeff=False):
		"""
		Subtract the best polynomial fit at each column of the data stored
		in .members.

		maxcores:		max number of parallel processes into which to split
						the fit on each column
		givencoeff:		not yet implemented
		indeplabels:	list of attribute names to use as independent variables
						in the polynomial fit.
		crossinds:		not yet implemented
		degree:			degree of polynomial to use in fit
		savecoeff:		not yet implemented.
		"""
		self.savecoeff = savecoeff
		#if isinstance(uncertainties,(int,float)):
		#	uncertainties = uncertainties*np.ones((self.nummembers,
		#										   self.specdim))
		#self.uncertainties = uncertainties
		self.polyindeps(degree=degree,indeplabels=indeplabels)
		if isinstance(crossinds,(list,np.ndarray)):
			pass # restrict self.indeps
		if not isinstance(givencoeff,(np.ndarray,str)):
			self.fitspec = np.array(ml.parallel_map(self.colpolyfit,
													range(self.members.shape[1]),
										  			numcores=maxcores))
			self.fitspec = self.fitspec.T
		elif isinstance(givencoeff,str):
			pass #read out coeff

		elif isinstance(givencoeff,np.ndarray):
			self.fitspec = indeps*givencoeff

		self.spectra -= self.fitspec
		return None


	def project(self):
		return None # project spectra onto array or read it from file to proj

	def addnoise(self):
		#self.uncertaintes = stuff
		return None



class abundances(makeclusters):

	def __init__(self):
		self.membertype='abundances'
		return None

	def get_centers(self,genfn=normalgeneration,instances=1,readdata=False,
				 filename=None,numcluster=20,maxcores=1,
				 elems=np.array([6,7,8,11,12,13,14,16,19,20]),
				 **kwargs):
		self.synfilename = 'test.hdf5'
		self.datafile = h5py.File(self.synfilename,'w')
		self.centergenfn = normalgeneration
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
		self.centers= np.zeros((self.instances,self.numcluster,self.numelem))
		# Update kwarg dictionary with parameters that are class attributes 
		kwargs.update({'num':self.numcluster,'numprop':self.numelem})

		# Create cluster center abundances for each requested instance
		for i in range(self.instances):
			# Get unique timestamp for this instance
			currenttime = gettimestr()
			self.timestamps[i] = currenttime
			# Generate cluster center abundances
			clustercenters = self.centergenfn(**kwargs)
			self.centers[i] = clustercenters

			# Assign attributes to HDF5 dataset
			dsetname = 'center_abundances_{0}'.format(currenttime)
			self.datafile[dsetname] = clustercenters
			centerinfo = self.datafile[dsetname]
			getwrapperattrs(centerinfo,self.centergenfn,kwargdict=kwargs)
			centerinfo.attrs['elemnames'] = self.elemnames
			centerinfo.attrs['atmnums'] = self.elems
		return None

	def get_members(self,genfn=normalgeneration,nummembers=20,
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
			self.numstars = np.sum(self.nummembers)

			# Create holder arrays
			self.abundances = np.zeros((self.instances,
										self.numstars,
										self.numelem))
			self.labels_true = -np.ones((self.instances,
										 self.numstars))
			
			# Deep copy kwargs so it can be updated with user specified kwargs
			kwargdict = copy.deepcopy(kwargs)
			for i in range(self.instances):

				# Get cluster centers for this instance
				centers = self.centers[i]

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



class spectra(makeclusters):

	def __init__(self):
		return None

	def get_centers(self,genfn=self.abundances().get_centers,**kwargs):
		return None

	def get_members(self,genfn=normalgeneration,**kwargs):
		return None

