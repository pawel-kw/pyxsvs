from distutils.core import setup
from setuptools import setup
# Work around mbcs bug in distutils.
# http://bugs.python.org/issue10945
import codecs
try:
    codecs.lookup('mbcs')
except LookupError:
    ascii = codecs.lookup('ascii')
    func = lambda name, enc=ascii: {True: enc}.get(name=='mbcs')
    codecs.register(func)

setup(name='pyxsvs',
      version = '0.1.0',
      description = 'Python X-ray Speckle Visibility analysis tools',
      author = 'Pawel Kwasniewski',
      author_email = 'pawel.kw@gmail.com',
      url = 'https://github.com/pawel-kw/pyxsvs',
      install_requires = ['numpy>=1.7','matplotlib>=1.3','pyFAI'],
      packages = ['pyxsvs'],
      package_dir = {'pyxsvs': 'pyxsvs-src'},
      entry_points = {
          'console_scripts': [
              'pyxsvs = pyxsvs.pyxsvs:main',
              'makemask = pyxsvs.makemask:main',
              'createStatic = pyxsvs.createStatic:main',
              'reportResults = pyxsvs.reportResults:main'
              ]
          },
       package_data = {
        '': ['*.txt', '*.rst']
        }
      )

# ###############################################################################
# Check for Fabio to be present of the system
# ###############################################################################
try:
    import fabio
except ImportError:
    print("""pyFAI needs fabIO for all image reading and writing.
This python module can be found on:
http://sourceforge.net/projects/fable/files/fabio""")
