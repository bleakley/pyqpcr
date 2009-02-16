from distutils.core import setup
import glob
import sys

if sys.platform == 'darwin':
    extra_options = dict( \
                         setup_requires=['py2app'] \
                        )

elif sys.platform == 'win32':
    extra_options = dict( \
                         setup_requires=['py2exe'],
                         data_files = [(r'mpl-data', 
    glob.glob(r'C:\Python25\Lib\site-packages\matplotlib\mpl-data\*.*')),
                                      (r'mpl-data', 
    [r'C:\Python25\Lib\site-packages\matplotlib\mpl-data\matplotlibrc']),
                                      (r'mpl-data\images',
    glob.glob(r'C:\Python25\Lib\site-packages\matplotlib\mpl-data\images\*.*')),
                                      (r'mpl-data\fonts',
    glob.glob(r'C:\Python25\Lib\site-packages\matplotlib\mpl-data\fonts\*.*')) ] \
                        )
else:
    extra_options = dict( \
                         data_files=[('share/icons', ['pyQPCR-16.png', 
                             'pyQPCR-32.png']), ('share/applications',
                              ['pyQPCR.desktop'])] \
                         )

setup(name='pyQPCR',
      version='0.1',
      description='a free software to deal qPCR',
      long_description='a qt4 based interface to deal qPCR',
      author='Thomas Gastine',
      author_email='thomas.gastine@wanadoo.fr',
      url='http://sourceforge.net/',
      licence='GPLv3',
      packages=['pyQPCR', 'pyQPCR.dialogs', 'pyQPCR.utils'],
      scripts=['scripts/qpcr'],
      **extra_options
      )
