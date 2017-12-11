from __future__ import print_function
import numpy as np		

def normalgeneration(num=10,numelem=10,centers=np.zeros(10),
					 stds=0.5*np.ones(10)):
	"""
	num: 			number of values to generate
	numelem:		number of elements
	centers:		center (mean value) of the numelem-dimensional Gaussian
	stds:			standard deviation of the numelem-dimensional Gaussian
	"""
	return np.random.randn(num,numelem)*stds + centers


def uniformgeneration(num=10,numelem=10,minvals=-0.1*np.ones(10),
					  maxvals=0.1*np.ones(10)):
	"""
	num:			number of values to generate
	numelem:		number of elements
	minvals:		min abundance value for each element
	maxvals:		max abundance value for each element
	"""
	return np.random.random(size=(num,numelem))*(maxvals-minvals) + minvals