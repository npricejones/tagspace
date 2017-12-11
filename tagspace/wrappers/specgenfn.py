from apogee.modelspec.turbospec import synth as tsynth
from apogee.modelatm import atlas9
from psm import psm

def apoturbospectrum(elemnums,elemvals,teff=4750,logg=2.5,metals=-0.25,
					 am=0.25,cm=0.25,linelist='201404080919',lsf='all',
					 cont='cannon',vmacro=6.,isotopes='solar'):	
	atm = atlas9.Atlas9Atmosphere(teff=teff,logg=logg,metals=metals,am=am,cm=cm)
	elemproperties = []
	for e in range(len(elemnums)):
		elemproperties.append([elemnums[e],elemvals[e]])
	synspec= tsynth(*elemproperties,modelatm=atm,linelist=linelist,lsf=lsf,
					cont=cont,vmacro=vmacro,isotopes=isotopes)
	return synspec


def psm(elemnums,elemvals,teff=4750,logg=2.5,vturb=1.5,c12c13=0):
	labels = np.zeros(19)
	labels[:3] = np.array([teff/1000.,logg,vturb])
	elemvals = elemvals[atmsort]
	elemnums = elemnums[atmsort]
	psmelems = np.array([6,7,8,11,12,13,14,16,19,20,22,23,25,26,28])
	psminds = np.array([3,4,5,6,7,8,9,10,11,12,13,14,15,16,17])
	for e in range(len(elemnums)):
		if elemnums[e] in psmelems:
			ind = psminds[psmelems==elemnums[e]]
			labels[ind] = elemvals[e]
	labels[18]=c12c13
	return psm(labels)