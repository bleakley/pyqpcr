from distutils.core import setup
import glob
import sys

if sys.platform == 'darwin':
    extra_options = dict(
                         setup_requires=['py2app']
                        )

elif sys.platform == 'win32':
    import py2exe
    extra_options = dict(
                         setup_requires=['py2exe'],
                         data_files = [(r'mpl-data', 
    glob.glob(r'C:\Python25\Lib\site-packages\matplotlib\mpl-data\*.*')),
                                      (r'mpl-data', 
    [r'C:\Python25\Lib\site-packages\matplotlib\mpl-data\matplotlibrc']),
                                      (r'mpl-data\images',
    glob.glob(r'C:\Python25\Lib\site-packages\matplotlib\mpl-data\images\*.*')),
                                      (r'mpl-data\fonts',
    glob.glob(r'C:\Python25\Lib\site-packages\matplotlib\mpl-data\fonts\*.*')) ]
                        )
    extra_options['options'] = \
           {
    'py2exe': { "includes" : ["sip", "matplotlib.backends",  
                               "matplotlib.backends.backend_qt4agg",
                               "matplotlib.figure","pylab", "numpy", 
                               "matplotlib.numerix.fft",
                               "matplotlib.numerix.linear_algebra", 
                               "matplotlib.numerix.random_array",
                               "matplotlib.backends.backend_tkagg"],
                'excludes': ['_gtkagg', '_tkagg', '_agg2', '_cairo', '_cocoaagg',
                             '_fltkagg', '_gtk', '_gtkcairo', ],
                'dll_excludes': ['libgdk-win32-2.0-0.dll',
                                 'libgobject-2.0-0.dll',
                                 'msvcp71.dll']
              }
          }
    extra_options['windows'] = [{'script' : 'scripts\qpcr',
                                 'icon_resources': [(1, 'logo.ico')]}]

else:
    import platform
    extra_options = dict(
          data_files=[('share/icons/hicolor/16x16/apps', ['pyQPCR-16.png']), 
                      ('share/icons/hicolor/32x32/apps', ['pyQPCR-32.png']), 
                      ('share/applications', ['pyQPCR.desktop'])] \
                        )
    if platform.dist()[0] == 'fedora':
        extra_options['options'] = \
                {
                'bdist_rpm': { 
                    'requires': ['python-matplotlib', 'PyQt4'],
                    'distribution_name': ['fedora']
                             },
                'install': {'optimize': '1', 'prefix' : ['/usr']}
                } 
    elif platform.dist()[0] == 'SuSe':
        extra_options['options'] = \
                {
                'bdist_rpm': { 
                    'requires': ['python-matplotlib', 'python-qt4'],
                    'distribution_name': ['opensuse']
                             },
                'install': {'prefix': ['/usr']}
                } 

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
