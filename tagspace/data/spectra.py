import numpy as np 
from tagspace.wrappers.genfns import normalgeneration
from tagspace.clusters.makeclusters import cluster_spectra
#from apogee.modelspec.turbospec import synth as tsynth
#from apogee.modelatm import atlas9
from psm import generate_spectrum as psmgen
from psm import reference_point as psmref
from galpy.util import multi as ml
from sklearn.preprocessing import PolynomialFeatures


# indexing and atomic number information for PSM
psmelems = np.array([6,7,8,11,12,13,14,16,19,20,22,23,25,26,28])
psmeleminds = np.array([3,4,5,6,7,8,9,10,11,12,13,14,15,16,17])
psmkeys = ['TEFF','LOGG','VTURB','C12C13']
psmkeyinds = np.array([0,1,2,18],dtype=int)

# maybe just want to assign all of the abundances psm knows about?
# what about cases with unknown abundances?
# then this framework makes less sense (i.e. not logical with from_center_spectrum)
class psmspectra(cluster_spectra):
	def __init__(self,members,photosphere,abundances):
		super(psmspectra,self).__init__()
		self.nummembers = members
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