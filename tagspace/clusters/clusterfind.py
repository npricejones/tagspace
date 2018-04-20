import os
import h5py
import numpy as np
from galpy.util import multi as ml
from astropy.io import fits
from sklearn.metrics import homogeneity_score
from tagspace import tagdir
from tagspace.data import gettimestr
from tagspace.wrappers.genfns import normalgeneration
from tagspace.data.spectra import spectra
from tagspace.clusters import external_validation
from tagspace.clusters.makeclusters import makeclusters

class tag(makeclusters):
   """
   Using existing data, find clusters using a user-specified cluster finding algorithm
   structured after the template of scikit-learn (e.g. a class with a fit_predict 
   function)

   """
   def __init__(self,genfn=normalgeneration,instances=1,readdata=False,
             filename=None,numcluster=20,numelem=10,maxcores=1,
             elems=np.array([6,7,8,11,12,13,14,16,19,20]),
             **kwargs):
      """
      Create or read cluster data.

      genfn:      Function used to find cluster center abundances 
               (defaults to choosing from a normal distribution)
      instances:  The number of clusters to create with the same 
               parameters (defaults to one)
      readdata:   Boolean that looks for existing data if set True 
               (redundant with filename******)
      filename:   If readdata, read cluster center information from this 
               path. Unless the path starts are root '/' or home '~', 
               assumes data lives in the TAGSPACE_DATA environment 
               variable
      numcluster: number of clusters to generate (defaults to 20)
      maxcores:   Maxmimum number of cores to use for parallel processes
      elems:      List or array of atomic numbers or symbols for elements
               to be generated (defaults to: carbon, nitrogen, oxygen,
               sodum, magnesium, aluminum, silicon, sulfur, potassium,
               and calcium)
      **kwargs:   Passed to genfn

      Returns None
      """
      super(tag,self).__init__(self,genfn=normalgeneration,instances=instances,readdata=readdata,
                         filename=filename,numcluster=numcluster,numelem=numelem,maxcores=maxcores,
                         elems=elems,**kwargs)
      return None

   def cluster_wrapper(self,i):
      """
      Wrapper to find clusters in parallel.

      i:    Index to the cluster created instance to use.
      """
      numstars = self.clusterdata[i].shape[0]
      numproperties = self.clusterdata[i].shape[1]
      repeatpreds = np.zeros((numstars*self.repeats))
      startrack = 0
      for r in range(self.repeats):
         predict = self.clusterfn(**self.kwargs).fit_predict(self.clusterdata[i])
         repeatpreds[startrack:startrack+numstars] = predict
         startrack+=numstars
      self.instance['labels_pred_{0}'.format(self.timestamps[i])] = repeatpreds
      preds = self.instance['labels_pred_{0}'.format(self.timestamps[i])]
      preds.attrs['data'] = self.datapath[i]
      getwrapperattrs(preds,self.clusterfn,kwargdict=kwargs)
      return repeatpreds

   def cluster(self,datapath,clusterfn,repeats=1,**kwargs):
      """
      Identify clusters in data file

      datapath:      Path to cluster data in hdf5 path
      clusterfn:     Function used to identify cluster (in the style of scikit learn)
      repeats:    Number of times to apply cluster finding to a given data set.
      **kwargs:      Pass to clusterfn

      Returns None

      """
      # Open relevant datafile
      self.clsfilename = tagdir+'cluster_find/'+self.clusterfn.__name__+'/clusterfinding.hdf5'
      self.clsdatafile = h5py.File(self.clsfilename,'w')
      self.clusterpath = clusterfn.__name__

      # Create group if it doesn't already exist
      if self.clusterpath not in self.datafile:
         self.instance = self.datafile.create_group(self.clusterpath)
      elif self.clusterpath in self.datafile:
         self.instance = self.datafile[self.clusterpath]

      # Find cluster data
      self.datapath = datapath
      self.syndatafile = h5py.File(self.synfilename,'w')
      self.clusterdata = self.syndatafile[datapath][:]
      self.clusterfn = clusterfn
      self.repeats = repeats
      self.kwargs = kwargs
      self.labels_pred = np.array(ml.parallel_map(self.cluster_wrapper,
                                       range(self.instances),
                                       numcores=self.maxcores))
      self.datafile.close()
      return None

   def extval(self,metric=homogeneity_score):
      self.extscores = np.zeros((self.instances,self.repeats,self.labels_true.shape[1])) 
      for i in range(self.instances):
         for r in range(self.repeats):
            self.extscores[i][r] = external_validation(self.labels_true[i],self.labels_pred[i][r],metric=metric)
      return None

   def intval():
      return None

   def violinstats():
      return None

