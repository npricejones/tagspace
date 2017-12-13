import numpy as np 
from tagspace.wrappers.genfns import normalgeneration
from apogee.modelspec.turbospec import synth as tsynth
from apogee.modelatm import atlas9
from psm import psm
from galpy.util import multi as ml
from sklearn.preprocessing import PolynomialFeatures
from tagspace.wrappers.specgenfns import turbospectrum,psm

# undefined: specdim

class spectra(object):
	def __init__(self):
		return None

	def from_center_spectrum(self,centerspec,genfn=normalgeneration):

		return None

	def from_center_abundances(self,centerspec,genfn=normalgeneration):

		return None

	def from_member_abundances(self,abundances):

		return None

	def pixelpolyfit(self,i):
		stars = self.spectra[:,i]
		covI = np.diag(1./self.uncertainties[:,i]**2)
		indeps = np.dot(self.indeps.T,np.dot(covI,self.indeps))
		stars = np.dot(self.indeps.T,np.dot(covI,stars))
		indepsI = np.linalg.inv(indeps)
		coeffs = np.dot(indepsI,stars)
		coeff_errs = np.array([np.sqrt(np.array(indepsI)[i][i]) for i in range(indeps.shape[1])])
		if self.savecoeff:
			pass # save output somewhere
		return indeps*coeffs

	def subpolyfit(self,maxcores=1,givencoeff=None,indeplabels=['TEFF'],
				   crossinds='all',degree=2,uncertainties=1,savecoeff=False):
		self.maxcores = maxcores
		self.savecoeff = savecoeff
		if isinstance(uncertainties,(int,float)):
			uncertainties = uncertainties*np.ones((self.nummembers,
												   self.specdim))
		self.uncertainties = uncertainties
		self.polynomial = PolynomialFeatures(degree=degree)
		indeps = np.zeros((self.nummembers,len(indeplabels)))
		for i in range(len(indeplabels)):
			indeps[:,i] = getattr(self,indeplabels[i])
			indeps[:,i] -= np.median(indeps[:,i])
		self.indeps = self.polynomial.fit_transform(indeps)
		if isinstance(crossinds,(list,np.ndarray)):
			pass # restrict self.indeps
		if not isinstance(givencoeff,(np.ndarray,str)):
			self.fitspec = np.array(ml.parallel_map(self.pixelpolyfit,
													range(self.specdim),
										  			numcores=self.maxcores))
		elif isinstance(givencoeff,str):
			pass #read out coeff

		elif isinstance(givencoeff,np.ndarray):
			self.fitspec = indeps*givencoeff

		self.spectra -= self.fitspec
		return None


	def project(self):
		return None # project spectra onto array or read it from file to proj


class psmspectra(spectra):
	def __init__(self,members,photosphere):
		super(psmspectra,self).__init__()
		self.nummembers = members
		self.TEFF = photosphere['TEFF']
		self.LOGG = photosphere['LOGG']
		self.VTURB = photosphere['VTURB']