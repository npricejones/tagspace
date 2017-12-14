import h5py
import numpy as np
from tagspace.clusters.makeclusters import makeclusters
from tagspace.wrappers.genfns import normalgeneration
from tagspace.data.spectra import psmspectra
from tagspace.clusters import external_validation
from sklearn.cluster import DBSCAN
from sklearn.metrics import homogeneity_completeness_v_measure
import matplotlib.pyplot as plt
from tqdm import tqdm
from galpy.util import multi as ml

# THEORETICAL INTRA CLUSTER SPREAD AT SNR=100                                                                                             
# From Ting et al 2016 arxiv:1602.06947                                                                                                   
ch_cls = 5e-3
nh_cls = 0.01
oh_cls = 0.01
nah_cls = 0.038
mgh_cls = 7.6e-3
alh_cls = 0.02
sih_cls = 8.2e-3
sh_cls = 0.024
kh_cls = 0.044
cah_cls = 0.016
tih_cls = 0.018
vh_cls = 0.06
mnh_cls = 0.013
nih_cls = 0.01
feh_cls = 4.3e-4
c12c13_cls = 0

tingspr = np.array([ch_cls,nh_cls,oh_cls,nah_cls,mgh_cls,alh_cls,sih_cls,sh_cls,kh_cls,
                    cah_cls,tih_cls,vh_cls,mnh_cls,nih_cls,feh_cls])

# INTRA CLUSTER SPREAD                                                                                                                    
# From 'global uncertainties' in Table 6 of Holtzmann et al 2015                                                                          

ch_cls = 0.035
nh_cls = 0.067
oh_cls = 0.050
nah_cls = 0.064
mgh_cls = 0.053
alh_cls = 0.067
sih_cls = 0.077
sh_cls = 0.063
kh_cls = 0.065
cah_cls = 0.059
tih_cls = 0.072
vh_cls = 0.088
mnh_cls = 0.061
feh_cls = 0.053
nih_cls = 0.060
c12c13_cls = 0

spreads = np.array([ch_cls,nh_cls,oh_cls,nah_cls,mgh_cls,alh_cls,sih_cls,sh_cls,kh_cls,
                    cah_cls,tih_cls,vh_cls,mnh_cls,nih_cls,feh_cls])


clusters = makeclusters(genfn=normalgeneration,instances=1,numcluster=50,maxcores=1,
                        elems=np.array([6,7,8,11,12,13,14,16,19,20,22,23,25,26,28]),
                        centers = np.zeros(15),stds = np.ones(15)*spreads*5)

clusters.create_abundances(genfn=normalgeneration,nummembers=100,
                           stds=tingspr)

mem = np.sum(clusters.nummembers)

teffdict = {'num':mem,'numprop':1,'centers':4500*np.ones(mem),'stds':100*np.ones(mem)}
loggdict = {'num':mem,'numprop':1,'centers':2.6*np.ones(mem),'stds':0.1*np.ones(mem)}
vturbdict = {'num':mem,'numprop':1,'centers':6*np.ones(mem),'stds':np.ones(mem)}
clusters.get_photosphere(nummembers=mem,params=[['TEFF'],['LOGG'],['VTURB']],
                         genfns=[normalgeneration,normalgeneration,
                                 normalgeneration],
                         kwargdicts=[teffdict,loggdict,vturbdict])
datafile = h5py.File(clusters.synfilename,'r+')
datakey = datafile['normalgeneration'].keys()[-1]
abundances = datafile['normalgeneration/'+datakey]
specinfo = psmspectra(mem,clusters.photosphere,abundances)
# Generation + fitting takes about 10 minutes with 5000 stars
specinfo.from_member_abundances(1) # < 1 min
print 'Done spectra form'
specinfo.subpolyfit(maxcores=4,indeplabels=['TEFF','LOGG'],
                   crossinds='all',degree=2,uncertainties=1)
print 'Done fit'
labels_true = abundances.attrs['labels_true']
metric = 'euclidean'
min_samples = 10
eps = np.arange(0.1,2.75,0.25)
labels_pred = -np.ones((len(eps),mem))
homscores = np.zeros(len(eps))
comscores = np.zeros(len(eps))
vscores = np.zeros(len(eps))

def DBSCANwrap(i):
    print 'starting {0} of {1}'.format(i+1,len(eps))
    return DBSCAN(min_samples=min_samples,metric=metric,
                  eps=eps[i]).fit_predict(specinfo.spectra)

labels_pred = ml.parallel_map(DBSCANwrap,range(len(eps)),numcores=4)

for e in tqdm(range(len(eps))):
    score = homogeneity_completeness_v_measure(labels_true,labels_pred[e])
    homscores[e],comscores[e],vscores[e] = score

plt.figure(figsize=(15,6))
plt.subplot(131)
plt.semilogx(eps,homscores)
plt.xlabel('epsilon')
plt.ylabel('homogeneity')
plt.subplot(132)
plt.semilogx(eps,comscores)
plt.xlabel('epsilon')
plt.ylabel('completeness')
plt.subplot(133)
plt.semilogx(eps,vscores)
plt.xlabel('epsilon')
plt.ylabel('harmonic mean')
plt.savefig('scores.pdf')

plt.figure(figsize=(15,10))
plt.subplot(111)
plt.imshow(specinfo.spectra,cmap='bwr',aspect=specinfo.spectra.shape[1]/float(specinfo.spectra.shape[0]),vmin=-0.5,vmax=0.5)
plt.colorbar()
plt.savefig('spectra2d.pdf')





