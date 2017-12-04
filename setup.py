from setuptools import setup
import warnings

setup(name='tagspace',
      version='1.0.0',
      description='Testing cluster finding on synthetic and actual chemical space',
      author='Natalie Price-Jones',
      author_email='price-jones@astro.utoronto.ca',
      url='https://github.com/NatalieP-J/spectralspace',
      package_dir = {'tagspace/': ''},
      packages=['tagspace','tagspace/examples','tagspace/data','tagspace/clusters'],
      package_data={'tagspace/data':['distributions/apogee/rc_photosphere_params.npy',
                                     'distributions/apogee/rg_photosphere_params.npy']},
      dependency_links = ['https://github.com/jobovy/apogee/tarball/master#egg=apogee',
                          'https://github.com/jobovy/galpy/tarball/master#egg=galpy',
                          'https://github.com/jobovy/isodist/tarball/master#egg=isodist'],
      install_requires=['numpy','scipy','matplotlib','astropy','h5py',
                        'apogee','galpy','isodist','scikit-learn',
                        'periodictable']
      )

warnings.warn('''APOGEE installation requires environment variables to be set: SDSS_LOCAL_SAS_MIRROR, RESULTS_VERS, APOGEE_APOKASC_REDUX; see the package documentation for more details.''')
warnings.warn('''This package requires the TAG_SPACE environment variable to set to a storage directory''')