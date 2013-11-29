from distutils.core import setup
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
      version='0.1',
      description='Python X-ray Speckle Visibility analysis tools',
      author='Pawel Kwasniewski',
      author_email='pawel.kw@gmail.com',
      url='https://github.com/pawel-kw/pyxsvs',
      py_modules=['pyxsvs','pyxsvs.makemask','pyxsvs.createStatic'],
      )
