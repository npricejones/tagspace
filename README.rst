tagspace
-----------
Package to generate synthetic chemical space and identify clusters in a parallel way.

.. contents::

AUTHOR
======

Natalie Price-Jones

INSTALLATION
============

Navigate into the repository to install with:

`python setup.py install`

with possible option for local installation only:

`python setup.py install --prefix=<some local directory>`

This repository requires the following packages: [numpy](http://www.numpy.org/), [matplotlib](http://matplotlib.org/), [scipy](https://www.scipy.org/), [sklearn](http://scikit-learn.org/stable/), [astropy](http://www.astropy.org/), [galpy](https://github.com/jobovy/galpy), [apogee](https://github.com/jobovy/apogee), and [isodist](https://github.com/jobovy/isodist))

If you wish to try the example notebooks, you'll also need [jupyter](http://jupyter.org)

Please continue to the section below for instructions on how to set relevant environment variables

ENVIRONMENT VARIABLES AND FILE STRUCTURES
=========================================

This package will save the synthetic chemical space it creates (unless explicitly instructed otherwise during ```makeclusters``` object creation with ```save=False```). Files will be saved in the directory set to the environment variable **CHEMICAL_SPACE_DATA**.

PURPOSE
=======

This package is intended to test the limits 'chemical tagging', the process of identifying groups of stars from a large observational dataset by looking for clusters in stellar chemistry. This can be done either by looking at observed spectra directly (which encode chemical information through absorption lines) or by distilling observed spectra into a lower dimensional chemical space either by projecting them along relevent principal components or by using model spectra to derive chemical abundances.

CHEMICAL TAGGING
^^^^^^^^^^^^^^^^
Chemical tagging has a great deal of potential as a tool to understand the history of our Milky Way Galaxy. In its most general applications, it can be used to divide our Galaxy into broad populations that can be associated with features like the thin and thick disks or the halo (e.g. `Hawkins et. al. 2015 <https://arxiv.org/abs/1507.03604>`), or to constrian the chemical evolution of our galaxy (e.g. `Jofre et. al. 2017 <https://arxiv.org/abs/1611.02575>`). The technique has also been used to identify new members of a known chemical population (e.g. `Martell et. al. 2016 <https://arxiv.org/abs/1605.05792>`) or identify unusual new subpopulations (e.g. `Schiavon et. al. 2017 <https://arxiv.org/abs/1606.05651>`) With sufficient precision, chemical tagging may also be able to identify 'birth clusters' of stars that were born in the same gas cloud, even if their dynamical history has erased all evidence of their shared past. However, this requires a blind approach to chemical tagging (e.g. `Mitschang et. al. 2014 <https://arxiv.org/abs/1312.1759>` `Hogg et. al. 2016 <https://arxiv.org/abs/1601.05413>`, ), and may be limited by not only observational precision but by the intrinsic chemical homogeneity and uniqueness of the birth clusters (see e.g. `Blanco-Cuaresma et. al. 2015 <https://arxiv.org/abs/1503.02082>` and `Bovy 2016 <https://arxiv.org/abs/1510.06745>` for studies of instrinic birth cluster properties where open clusters are used as birth cluster proxies).

APPROACH
^^^^^^^^
Motivated by the studies listed above, I am developing this pacakge to test the efficacy of blind chemical tagging on various chemical spaces. I am primarily interested in identifying the conditions under which the groups of stars identified by a cluster-finding routine may be associated with stellar birth clusters. To that end, I am creating a flexible approach to cluster formation, allowing user specified functions to define how cluster members are chosen in chemical space. The code is built to create many instances of a particular chemical space and identify their clusters simultaneously. Success is measured by end result homogeneity (`scikit-learn implementation <http://scikit-learn.org/stable/modules/clustering.html#homogeneity-completeness-and-v-measure>`), which describes the level to which each of the algorithimically identified clusters contains only members of a single user-created cluster.

EXAMPLE USAGE
=============

This repository includes several notebooks in the ```examples``` folder that demonstrate more involved usage of the package.

BASIC WORKFLOW
^^^^^^^^^^^^^^

In general the workflow follows a few steps:

Making synthetic cluster data
+++++++++++++++++++++++++++++

Start by importing the repository's makecluster class object. You will also need to choose two generation functions: one to find the cluster centers and another to find members of a cluster. For this example, we'll use a normal distribution for both finding both cluster centers and members.

		import numpy as np
		from tagspace.clusters.makecluster import makecluster, normalgeneration

We'll use ```normalgeneration``` to find our cluster centers. This function takes three arguments: the number of clusters to identify, the mean of the normal distribution (i.e. the center of chemical space) and the standard deviation of the normal distribution. The latter two arguments may have dimensionality of your choosing. In this case we'll assume we're working with 10 chemical elements and want to input 20 clusters. We give the function and its kwargs to ```makeclusters```

		clusters = makeclusters(genfn=normalgeneration,num = 20, means = np.zeros(10), stds = 0.5*np.ones(10))

We have created our cluster centers. ```makeclusters``` has also automatically generated a directory associated with this data set, as well as a root string for saving individual cluster instances. We can overwrite these by passing the ```basepath``` and ```basename``` kwargs to change the directory and root name respectively.

We now have access to the function associated with ```makeclusters```, one of which is ```create_abundances```. This function will generate chemical abundances for members of the clusters given a function to use to find members and its kwargs. We'll use ```normalgeneration``` again, and give each cluster 15 members.

		clusters.create_abundances(genfn = normalgeneration, num = 15, means = cluster.centers, stds = 0.05*np.ones(10))

Since we're using ```normalgeneration``` and have given the ```means``` kwarg as an array with 20 rows (the number of clusters) and 10 columns (the number of chemical abundances), we will create 15 members for each of the 20 clusters. We could specify a different number of members for each cluster by changing our ```num``` kwarg to be an array with length 20.

With this we've created a very simple chemical space. Our abundances are in the array ```clusters.abundances```. We also have the array ```clusters.labels_true```, which tells us which original cluster each set of abundances (which correspond to a star) belong to.

Running cluster finding algorithms
++++++++++++++++++++++++++++++++++

Our next step is to call our cluster finding algorithm and apply it to our data. For this simple case, we'll use the wrapper for ```scikit-learn```'s KMeans algorithm. First we create a ```tag``` object, which takes a ```makeclusters``` object.

		from tagspace.clusters.clusterfind import tag
		tagclusters = tag(clusterdata=clusters)

Our ```tagclusters``` now has the properties of ```clusters``` as well as an array of zeros in ```tagclusters.labels_pred```. This is where we will store the indices that divide our stars into clusters according to the cluster finding algorithm we choose. We now run kmeans, which requires the number of clusters to find as input. We'll choose it to be 20, the true number of clusters.

		tagclusters.kmeans(n_clusters=20)

To see all of kmeans possible kwargs, run ```help(tagclusters.kmeans())```.

This function has now updated our ```tagclusters.labels_pred``` with the labels according to ```kmeans```. We could have used one of the other included wrappers or written our own by passing it through ```tagcluster.customfn(clusterfn = <name of function>,<kwargs>)```

Measuring success
+++++++++++++++++

Now that we have a prediction for how our data should be divided into clusters, we'd like to measure our level of success. We'll use the wrapper for ```sklearn.metric.homogeneity_score``` to compute this.

		tagclusters.homogeneity()

This function

MORE COMPLICATED CHEMICAL SPACES
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


SCALING UP
^^^^^^^^^^

In addition to using more complicated chemical spaces, we may also wish to scale up our analysis so we avoid relying on any individual cluster instance, which may be dominated by unusual cluster distributions. To achieve this, we simply ```makeclusters``` the ```instances``` kwarg. This is set to 1 by default. Choosing a higher number will create multiple cluster instances. Subsequent functions for cluster finding and success measurement know about the shape of the clusters and so can divide the resulting data appropriately.

The operations required to create and later find clusters in multiple instances of a data set automatically use all available cores. These can be constrained to a fixed value by setting the ```cores``` kwarg when creating a ```makeclusters``` object or by manually updating the variable in between function calls with ```<makeclusters object name>.cores = <integer>```. 

The cluster finding functions included in the ```tag``` object also support multiple cluster finding attempts through the ```repeats``` kwarg. Setting this to an integer will also automatically distribute processes to all possible cores.

USING EXISTING DATA
^^^^^^^^^^^^^^^^^^^