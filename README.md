# pyQPCR History

This project is a fork of pyQCPR, originally created by Thomas Gastine and Magali Hennion. It was originally hosted at http://pyqpcr.sourceforge.net/.

# pyQPCR Installation

## Dependencies

Before installing pyQPCR, you need to have these libraries:

 * PyQt 4.3 with Qt 4.2 or newer and SIP 4.19 or newer
   * https://www.riverbankcomputing.com/software/sip/download
   * https://www.riverbankcomputing.com/software/pyqt/download
   * https://download.qt.io/archive/qt/4.2/

 * python-matplotlib 0.98 or newer
   * http://matplotlib.sourceforge.net/

 * numpy 1.1 or newer
   * http://numpy.scipy.org/

 * scipy 0.7 or newer
   * http://www.scipy.org

## Running From the Source Tree

If you want to run pyQPCR from the source directory without installing, just run

  `direct-run.py`

## Installation

### Source Installation
After installing the dependencies, you can install pyQPCR by running:

  `sudo python setup.py install`

This will automatically build and install all required Python modules. To
start pyQPCR now you can use:

  `qpcr`

### RPM based distribution
You can build the rpm of pyQPCR by running:

  `python setup.py bdist_rpm`

This will automatically build the rpm in the bdist directory. You can then install it by running:

  `sudo rpm -ivh dist/pyQPCR-0.1-1.noarch.rpm`

This installation process is known to work on Fedora 7, 8 (i686 and x86_64) and
OpenSuSE 11.1 (i686)
