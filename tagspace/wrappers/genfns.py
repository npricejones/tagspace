import numpy as np

def normalgeneration(num=10,numelem=10,centers=np.zeros(10),
					 stds=0.5*np.ones(10)):
	"""
	num: 			number of values to generate
	numelem:		number of elements
	means:			means around which to generate values 
					if shape is numelem, then use same 
	"""
	return np.random.randn(num,numelem)*stds + centers


def uniformgeneration():
	return None