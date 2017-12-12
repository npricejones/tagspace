from __future__ import print_function
import numpy as np		

def normalgeneration(num=10,numprop=10,centers=np.zeros(10),
					 stds=0.5*np.ones(10)):
	"""
	num: 			number of values to generate
	numprop:		number of properties per value
	centers:		center (mean value) of the numprop-dimensional Gaussian
	stds:			standard deviation of the numprop-dimensional Gaussian

	Returns an array of length num by numprop
	"""
	return np.random.randn(num,numprop)*stds + centers


def uniformgeneration(num=10,numprop=10,minvals=-0.1*np.ones(10),
					  maxvals=0.1*np.ones(10)):
	"""
	num:			number of values to generate
	numprop:		number of properties per value
	minvals:		min abundance value for each element
	maxvals:		max abundance value for each element

	Returns an array of length num by numprop
	"""
	return np.random.random(size=(num,numprop))*(maxvals-minvals) + minvals

def choice2dgeneration(num=10,numprop=2,choicearray=np.ones(100,3),
					   readdata=False,**kwargs):
	"""
	num:			number of values to generate
	numprop:		number of properties per value
	choicearray:	array of star properties to use to generate the histogram
					must have number of columns equal to numprop
	readdata:		option to specify a file location for reading the 
					choicearray
	**kwargs:		passed to numpy.histogramdd 

	Returns an array of length num by numprop
	"""
	if readdata:
		pass # do something to make choice array
	stars = np.zeros(num,numprop)
	hist,edges = np.histogram(choicearray[:,0],**kwargs)
	p1probs = hist[0]/np.sum(hist[0])
	p1choices =  (2*edges[0]-np.roll(edges[0],1))[1:]
	p1 = np.random.choice(p1choices,size=num,p=p1probs)
	stars[:,0] = p1
	starpos = 0
	for e in range(len(edges)-1):
		inds = np.where((edges[e] < choicearray[:,0]) & (choicearray[:,0] < edges[e+1]))
		hist,edges = np.histogram(choicearray[:,1][inds],**kwargs)
		p1matches = len(np.where((edges[e] < p1) & (p1 < edges[e+1]))[0])
		p2probs = hist[0]/np.sum(hist[0])
		p2choices =  (2*edges[0]-np.roll(edges[0],1))[1:]
		p2 = np.random.choice(p2choices,size=p1matches,p=p2probs)
		stars[:,1][starpos:starpos+p1matches] = p2
		starpos+=p1matches
	return stars


def lineargeneration(num=10,numprop=10,slope=1,intercept=0,
					 indeps=np.arange(10),indepfn=np.log10,**kwargs):
	"""
	num:			number of values to generate
	numprop:		number of properties per value
	slope:			slope of linear function
	intercept:		intercept of linear function
	indeps:			independent variable to use (must have len num or 
					be scalar)
	indepfn:		function to apply to independent variable
	**kwargs:		kwargs passed to indepfn

	"""
	if isinstance(indeps,(int,float)):
		indeps = np.ones(num)*indeps
	return slope*indepfn(indeps,**kwargs) + intercept