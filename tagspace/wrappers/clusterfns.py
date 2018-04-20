import numpy as np
# put user made clustering algorithms here to transform them into sklearn esque class objects

class RandomAssign(object):
   """
   Randomly assigns all stars to a cluster. If the number of clusters is not specified, 
   the algorithm chooses a random number less than the total number of stars.
   """
   def __init__(self,k = None):
      """
      Assigns the number of clusters.

      k:    The number of clusters to assign to (can be None but 
            this doesn't correspond to no clusters)
      """
      self.k = k
      return None

   def fit_predict(self,clusterdata):
      """
      Assigns each star in clusterdata to a random group. Will also update k
      to be a random number less than the total number of samples if k is 
      unspecified

      clusterdata:   array of abundance or spectra information to cluster

      Returns cluster assignments.
      """
      numstars = clusterdata.shape[0]
      if not self.k:
         self.k = np.random.randint(0,numstars,size=1)
      self.labels_pred = np.random.randint(0,self.k,size=numstars)
      return self.labels_pred