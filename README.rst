tagspace
-----------
Package to generate synthetic chemical space and identify clusters in a parallel way.

.. contents::
.. highlight:: python

AUTHOR
======

Natalie Price-Jones

INSTALLATION
============

Navigate into the repository to install with:

``python setup.py install``

with possible option for local installation only:

``python setup.py install --prefix=<some local directory>``

This repository requires the following packages: `numpy <http://www.numpy.org/>`__, `matplotlib <http://matplotlib.org/>`__, `scipy <https://www.scipy.org/>`__, `sklearn <http://scikit-learn.org/stable/>`, `astropy <http://www.astropy.org/>`__, `galpy <https://github.com/jobovy/galpy>`__, `apogee <https://github.com/jobovy/apogee>`__, and `isodist <https://github.com/jobovy/isodist>`__.

If you wish to try the example notebooks, you'll also need `jupyter <http://jupyter.org>`__

Please continue to the section below for instructions on how to set relevant environment variables

ENVIRONMENT VARIABLES AND FILE STRUCTURES
=========================================

ENVIRONMENT VARIABLES
^^^^^^^^^^^^^^^^^^^^^

This package will save the synthetic chemical space it creates (unless explicitly instructed otherwise during ``makeclusters`` object creation with ``save=False``). Files will be saved in the directory set to the environment variable **TAGSPACE_DATA**.

FILE STRUCTURES
^^^^^^^^^^^^^^^
Data created with the ``tagspace.clusters.makeclusters`` class for generating synthetic chemical spaces is saved in appropriately named directories, with the following structure
:: 
		$TAGSPACE_DATA/<cluster center generation function>/<cluster member generation function>/

Files are timestamped to avoid accidental overwriting.


OUTPUT DATA
^^^^^^^^^^^
Synthetic data created with ``tagspace.clusters.makeclusters`` is saved by default into .fits files. In their headers they contain the information about the functions used to create the cluster, including the parameters used for reproducability. The first file extension contains the cluster assignments, the second contains abundance data (if generated), the third contains spectra (if generated) and subsequent extensions contain information about possible changes to the spectra.

PURPOSE
=======

This package is intended to test the limits 'chemical tagging', the process of identifying groups of stars from a large observational dataset by looking for clusters in stellar chemistry. This can be done either by looking at observed spectra directly (which encode chemical information through absorption lines) or by distilling observed spectra into a lower dimensional chemical space either by projecting them along relevent principal components or by using model spectra to derive chemical abundances.

CHEMICAL TAGGING
^^^^^^^^^^^^^^^^
Chemical tagging has a great deal of potential as a tool to understand the history of our Milky Way Galaxy. In its most general applications, it can be used to divide our Galaxy into broad populations that can be associated with features like the thin and thick disks or the halo (e.g. `Hawkins et. al. 2015 <https://arxiv.org/abs/1507.03604>`__), or to constrian the chemical evolution of our galaxy (e.g. `Jofre et. al. 2017 <https://arxiv.org/abs/1611.02575>`__). The technique has also been used to identify new members of a known chemical population (e.g. `Martell et. al. 2016 <https://arxiv.org/abs/1605.05792>`__) or identify unusual new subpopulations (e.g. `Schiavon et. al. 2017 <https://arxiv.org/abs/1606.05651>`__) With sufficient precision, chemical tagging may also be able to identify 'birth clusters' of stars that were born in the same gas cloud, even if their dynamical history has erased all evidence of their shared past. However, this requires a blind approach to chemical tagging (e.g. `Mitschang et. al. 2014 <https://arxiv.org/abs/1312.1759>`__, `Hogg et. al. 2016 <https://arxiv.org/abs/1601.05413>`__), and may be limited by not only observational precision but by the intrinsic chemical homogeneity and uniqueness of the birth clusters (see e.g. `Blanco-Cuaresma et. al. 2015 <https://arxiv.org/abs/1503.02082>`__ and `Bovy 2016 <https://arxiv.org/abs/1510.06745>`__ for studies of instrinic birth cluster properties where open clusters are used as birth cluster proxies).

APPROACH
^^^^^^^^
Motivated by the studies listed above, I am developing this pacakge to test the efficacy of blind chemical tagging on various chemical spaces. I am primarily interested in identifying the conditions under which the groups of stars identified by a cluster-finding routine may be associated with stellar birth clusters. To that end, I am creating a flexible approach to cluster formation, allowing user specified functions to define how cluster members are chosen in chemical space. The code is built to create many instances of a particular chemical space and identify their clusters simultaneously. Success is measured by end result homogeneity (`scikit-learn implementation <http://scikit-learn.org/stable/modules/clustering.html#homogeneity-completeness-and-v-measure>`__), which describes the level to which each of the algorithimically identified clusters contains only members of a single user-created cluster.

EXAMPLE USAGE
=============

This repository includes several notebooks in the ``examples`` folder that demonstrate more involved usage of the package.

BASIC WORKFLOW - CLUSTERING SYNTHETIC DATA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In general the workflow follows a few steps:


Making synthetic cluster data
+++++++++++++++++++++++++++++

Start by importing the repository's makecluster class object. You will also need to choose two generation functions: one to find the cluster centers and another to find members of a cluster. For this example, we'll use a normal distribution for both finding both cluster centers and members.
::
		import numpy as np
		from tagspace.clusters.makeclusters import makeclusters, normalgeneration

We'll use ``normalgeneration`` to find our cluster centers. This function takes three arguments: the number of clusters to identify, the mean of the normal distribution (i.e. the center of chemical space) and the standard deviation of the normal distribution. The latter two arguments may have dimensionality of your choosing. In this case we'll assume we're working with 10 chemical elements and want to input 20 clusters. We give the function and its kwargs to ``makeclusters``
::
		clusters = makeclusters(genfn=normalgeneration,num = 20, means = np.zeros(10), stds = 0.5*np.ones(10))

We have created our cluster centers. ``makeclusters`` has also automatically generated a directory associated with this data set, as well as a root string for saving individual cluster instances. We can overwrite these by passing the ``basepath`` and ``basename`` kwargs to change the directory and root name respectively.

We now have access to the function associated with ``makeclusters``, one of which is ``create_abundances``. This function will generate chemical abundances for members of the clusters given a function to use to find members and its kwargs. We'll use ``normalgeneration`` again, and give each cluster 15 members.
::
		clusters.create_abundances(genfn = normalgeneration, num = 15, means = cluster.centers, stds = 0.05*np.ones(10))

Since we're using ``normalgeneration`` and have given the ``means`` kwarg as an array with 20 rows (the number of clusters) and 10 columns (the number of chemical abundances), we will create 15 members for each of the 20 clusters. We could specify a different number of members for each cluster by changing our ``num`` kwarg to be an array with length 20.

With this we've created a very simple chemical space. Our abundances are in the array ``clusters.abundances``. We also have the array ``clusters.labels_true``, which tells us which original cluster each set of abundances (which correspond to a star) belong to.

Running cluster finding algorithms
++++++++++++++++++++++++++++++++++

Our next step is to call our cluster finding algorithm and apply it to our data. For this simple case, we'll use the wrapper for ``scikit-learn``'s KMeans algorithm. First we create a ``tag`` object, which takes a ``makeclusters`` object.
::
		from tagspace.clusters.clusterfind import tag
		tagclusters = tag(data=clusters)

Our ``tagclusters`` now has the properties of ``clusters`` as well as an array of zeros in ``tagclusters.labels_pred``. This is where we will store the indices that divide our stars into clusters according to the cluster finding algorithm we choose. We now run kmeans, which requires the number of clusters to find as input. We'll choose it to be 20, the true number of clusters.
::
		tagclusters.kmeans(datatype='abundances',n_clusters=20)

To see all of kmeans possible kwargs, run ``help(tagclusters.kmeans())``.

This function has now updated our ``tagclusters.labels_pred`` with the labels according to ``kmeans``. We could have used one of the other included wrappers or written our own by passing it through ``tagcluster.customfn(clusterfn = <name of function>,<kwargs>)``

Measuring success
+++++++++++++++++

Now that we have a prediction for how our data should be divided into clusters, we'd like to measure our level of success. We'll use the wrapper for ``sklearn.metric.homogeneity_score`` to compute this.
::
		tagclusters.external.homogeneity()

This function measures the extent to which members of a cluster found by our chosen algorithm belonged to the same original cluster, so a value around 1 indicates successful clustering.

MORE COMPLICATED SYNTHETIC CHEMICAL SPACES
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Using chemical abundances as axes is the most common and straightforward approach to constructing a chemical space. However, we may wish to examine different versions of chemical space, and the construction of many of these is supported by ``tagspace``.

Spectra
+++++++

``tagspace`` supports two ways of generating spectra of stars to be members of a cluster. As in the previous section, we'll begin by creating a ``makeclusters`` object called ``clusters``. We will use the same simple case described above to find our cluster centers.
::		
		import numpy as np
		from tagspace.clusters.makeclusters import makeclusters, normalgeneration, uniformgeneration
		clusters = makeclusters(genfn=normalgeneration,num = 20, means = np.zeros(10), stds = 0.5*np.ones(10))


It is now necessary to specify photospheric parameters of our stars so we can generate the spectra. Unlike chemical abundances we do not expect these parameters to be similar for members of the same clusters. There are few options for how to get these parameters. If we already have this information, we can add it to our clusters information by passing it as a structured array to ``updateinfo``.
::
		clusters.updateinfo(<structured array of star information>)

The structured array must have labels ``'TEFF'``, ``'LOGG'`` and ``'VTURB'``.  The above is a convenience function; we could also have set ``clusters.spectra.teff``, ``clusters.spectra.logg`` and ``clusters.spectra.vturb`` to the appropriate values.

We have two other options to create the photospheric parameters - we can choose them in a related way with ``choosephotosphere`` or we can choose them independently using the function for each (e.g. ``chooseteff`` for selecting effective temperatures). For simplicity, let's choose the photospheric parameters simultaneously from a uniform distribution.
::
		clusters.choosephotosphere(genfn=uniformgeneration,bounds={'TEFF':[4000,5000],'LOGG':[2.0,4.0],'VTURB:[0.5,3.0]'})

Our first approach to generating member spectra begins by generating abundances and using those to create spectra. Start by identifying abundances.
::
		clusters.create_abundances(genfn = normalgeneration, num = 15, means = cluster.centers, stds = 0.05*np.ones(10),atmnum=[6,7,8,11,12,13,14,16,20,26])

We have added a new kwarg to ``create_abundances``; ``atmnum`` specifies which elements we are generating, since this is needed for spectra generation. We now have everything we need to create the spectra.
::
		clusters.create_spectra_abundances()

Alternatively, we can create a spectrum for each cluster center and vary it according to a generation function, in much the same way as we chose members in abundance space:
::
		cluster.create_spectra(genfn = normalgeneration, num = 15, means = cluster.centers, stds = 0.01*np.ones(10))

Finding clusters in spectral space is then a matter giving our cluster data to ``tag`` as before.


Fitting spectra
"""""""""""""""

Once spectra have been created, their use in chemical tagging can be improved by performing fits to remove differences between spectra due to differing photospheric parameters. To do this with ``tagspace``, use the function associated with the ``spectra`` object. If we assume we have created the ``clusters`` object from the previous section we can perform a fit in the following way. Let us assume we are interested in doing a second order polynomial fit in effective temperature, surface gravity and iron abundance with all cross terms included.
::
		clusters.spectra.fit(fitfn=polynomial,degree=2,variables=(clusters.spectra.teff,clusters.spectra.logg,clusters.spectra.abun['Fe']),crossterms=True)

This function updates the ``clusters.spectra.specs`` object and will save the new dataset.


Projecting spectra
""""""""""""""""""

We may wish to reduce the dimensionality of our spectra by projecting them along dimensions we think are important. We can supply a path to vectors describing these dimensions or provide them as an array. Either way we use ``project`` to do this in the following way.
::
		clusters.spectra.project(fname='<path to axis vectors>')

This function updates the ``clusters.spectra.specs`` object and will save the new dataset.

SCALING UP
^^^^^^^^^^

In addition to using more complicated chemical spaces, we may also wish to scale up our analysis so we avoid relying on any individual cluster instance, which may be dominated by unusual cluster distributions. To achieve this, we give ``makeclusters`` the ``instances`` kwarg. This is set to 1 by default. Choosing a higher number will create multiple cluster instances. Subsequent functions for cluster finding and success measurement know about the shape of the clusters and so can divide the resulting data appropriately.

The operations required to create and later find clusters in multiple instances of a data set automatically use all available cores. These can be constrained to a fixed value by setting the ``cores`` kwarg when creating a ``makeclusters`` object or by manually updating the variable in between function calls with ``<makeclusters object name>.cores = <integer>``. 

The cluster finding functions included in the ``tag`` object also support multiple cluster finding attempts through the ``repeats`` kwarg. Setting this to an integer will also automatically distribute processes to all possible cores.

USING EXISTING DATA
^^^^^^^^^^^^^^^^^^^

``tagspace`` is built to allow quick reproduction of previous results, as well as applications to non-synthetic datasets.

Data created with ``tagspace``
++++++++++++++++++++++++++++++

If we would like to work with previously created data in a new session, we will still need to create a ``makeclusters`` object and change its ``readdata`` kwarg from its default ``False`` to ``True``. We will also need to point ``makeclusters`` to the appropriate type of data. For example, if we wanted to use a specific file, we would give the ``fname`` kwarg with the path to the data (if this does not start from root ``/`` or from home ``~``, it is assumed to have the environment variable **CHEMICAL_SPACE_DATA** as its root). In this case our call would look like:
::
		from tagspace.clusters.makeclusters import makeclusters
		clusters = makeclusters(readdata=True,fname='<path>')

The ``fname`` kwarg also accepts a list or array of paths as input. If ``makeclusters``'s ``separate`` kwarg is set to ``False``, the stellar data are checked for shape and combined ,and initial clusters are appropriately reindexed.

Alternatively, if we wanted to use all data that was created with a particular generation function, our process takes an additional step. We will also need to specify what sort of stellar data we are looking for (either ``abundances``, ``spectra``, ``projspectra``, ``fitspectra`` or some list combining two or more of the proceeding), as well as the function used to generate the members for that data. Let's assume we are looking for all ``spectra`` and ``fitspectra`` data created with ``normalgeneration``. 
::
		from tagspace.clusters.makeclusters import makeclusters
		clusters = makeclusters(readdata=True,genfn=normalgeneration,separate=True)
		makeclusters.finddata(genfn=normalgeneration,datatype=['spectra','fitspectra'])

This will find all files that meet our criteria. The ``finddata`` function has additional options if, for example, we wanted to specify we were looking only for data where 10 members were created per cluster, or with particular standard deviation values.

External data
+++++++++++++

Data with known cluster assignments
"""""""""""""""""""""""""""""""""""

Data not created with ``tagspace``  but with known cluster assignments can be read in much the same way as previously created ``tagspace`` data, by using the ``fname`` kwarg of ``makeclusters`` to specify a path. Data should be in the form of a ``tagspace``-like .fits file (described in the Output Data subsection). The minimum requirements are a list of lists of data and a list of lists of cluster assignments with. The convenience function ``convert_to_TSfits`` in ``tagspace.data`` can easily convert the array (either from the current session or from file) into an appropriate fits file.
::
		from tagspace.data import convert_to_TSfits
		from tagspace.clusters.makeclusters import makeclusters
		convert_to_TSfits(<list of lists of star data>, <list of lists of cluster assignments>, datatype=<datatype>, fname='<path>')
		clusters = makeclusters(readdata=True,fname='<path>')

Here ``<datatype>`` refers to any of ``'abundances'``, ``'spectra'``, ``'projspectra'``, or ``'fitspectra'``

Data with unknown cluster assignments
"""""""""""""""""""""""""""""""""""""

Without known cluster assignments, we give our data directly to ``tag``,
::
		from tagspace.clusters.clusterfind import tag
		tagclusters = tag(data=<array of data to tag>)

and make use of the usual functions to create cluster assignments. Additional information about the stars (e.g. effective temperature, surface gravity), can be passed to ``tag`` as a structured array through the ``starinfo`` kwarg.


