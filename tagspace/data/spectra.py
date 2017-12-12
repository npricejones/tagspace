import numpy as np 
from tagspace.wrappers.genfns import normalgeneration
from apogee.modelspec.turbospec import synth as tsynth
from apogee.modelatm import atlas9
from psm import psm
from tagspace.wrappers.specgenfns import turbospectrum,psm

class psmspectra(object):
	def __init__(self,members,photosphere):
		self.nummembers = members
		self.TEFF = photosphere['TEFF']
		self.LOGG = photosphere['LOGG']
		self.VTURB = photosphere['VTURB']

	def from_center_spectrum(self,centerspec,genfn=normalgeneration):

	def from_center_abundances(self,centerspec,genfn=normalgeneration):

	def from_member_abundances(self,abundances):