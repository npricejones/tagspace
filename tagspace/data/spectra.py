import numpy as np 
from tagspace.wrappers.genfns import normalgeneration
from tagspace.clusters.makeclusters import makeclusters
#from apogee.modelspec.turbospec import synth as tsynth
#from apogee.modelatm import atlas9
from psm import generate_spectrum as psmgen
from psm import reference_point as psmref
from galpy.util import multi as ml
from sklearn.preprocessing import PolynomialFeatures

# undefined: specdim

class spectra(object):
	def __init__(self):
		return None

	def from_center_spectrum(self,centerspec,genfn=normalgeneration):

		return None

	def from_center_abundances(self,centerspec,genfn=normalgeneration):

		return None

	def from_member_abundances(self,maxcores):
		"""
		abundances:		h5py dataset with attribute atmnum specifying the
						atomic number of each column
		"""
		self.spectra = np.array(ml.parallel_map(self.genspec,
												range(self.nummembers),
										  		numcores=maxcores))
		return None


	def polyindeps(self,degree=2,indeplabels=['TEFF']):
		self.polynomial = PolynomialFeatures(degree=degree)
		indeps = np.zeros((self.nummembers,len(indeplabels)))
		for i in range(len(indeplabels)):
			indeps[:,i] = getattr(self,indeplabels[i].lower())
			indeps[:,i] -= np.median(indeps[:,i])
		self.indeps = np.matrix(self.polynomial.fit_transform(indeps))
		return None

	def pixelpolyfit(self,i):
		stars = np.matrix(self.spectra[:,i])
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
		self.savecoeff = savecoeff
		if isinstance(uncertainties,(int,float)):
			uncertainties = uncertainties*np.ones((self.nummembers,
												   self.specdim))
		self.uncertainties = uncertainties
		self.polyindeps(degree=degree,indeplabels=indeplabels)
		if isinstance(crossinds,(list,np.ndarray)):
			pass # restrict self.indeps
		if not isinstance(givencoeff,(np.ndarray,str)):
			self.fitspec = np.array(ml.parallel_map(self.pixelpolyfit,
													range(self.specdim),
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
		return None


# indexing and atomic number information for PSM
psmelems = np.array([6,7,8,11,12,13,14,16,19,20,22,23,25,26,28])
psmeleminds = np.array([3,4,5,6,7,8,9,10,11,12,13,14,15,16,17])
psmkeys = ['TEFF','LOGG','VTURB','C12C13']
psmkeyinds = np.array([0,1,2,18],dtype=int)

# maybe just want to assign all of the abundances psm knows about?
# what about cases with unknown abundances?
# then this framework makes less sense (i.e. not logical with from_center_spectrum)
class psmspectra(spectra):
	def __init__(self,members,photosphere,abundances):
		super(psmspectra,self).__init__()
		self.nummembers = members
		self.specdim = 7214
		self.abundances = abundances
		for k in range(len(psmkeys)):
			try:
				setattr(self,psmkeys[k].lower(),photosphere[psmkeys[k]])
			except ValueError:
				setattr(self,psmkeys[k].lower(),
					    np.ones(self.nummembers)*psmref[psmkeyinds[k]])

	def genspec(self,star):
		"""
		star:		index of abundances to use for spectra creation
					must match indexes of photosphere array in init

		returns 1 APOGEE-like (7214) spectrum 
		"""
		labels = np.copy(psmref) #defaults to reference point if no value given
		labels[:3] = np.array([self.teff[star]/1000.,self.logg[star],
							   self.vturb[star]])
		atmnums = self.abundances.attrs['atmnums']
		for a in range(len(atmnums)):
			if atmnums[a] in psmelems:
				ind = psmeleminds[psmelems==atmnums[a]]
				labels[ind] = self.abundances[:][star][a]
		labels[18]=self.c12c13[star]
		return psmgen(labels)