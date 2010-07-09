Download and install
********************

Download
========

The latest version of pyQPCR can be downloaded on the sourceforge webpage:
http://sourceforge.net/projects/pyqpcr/files/
You can also get the latest SVN [#f1]_ version with::

   svn co https://pyqpcr.svn.sourceforge.net/svnroot/pyqpcr



Run from a local directory
==========================

A script is given in the tarball to run directly pyQPCR without installing. 
Just open a console (or cmd on Windows) in the pyQPCR directory and simply run::

   direct-run.py


Install
=======

Windows
-------

Just download the .exe on the download page and install it. 
Everything should work flawlessly. You can also use the source code of pyQPCR
but you'll then need to install python explicitely on your windows install.
In such a case, `pythonxy <http://www.pythonxy.com>`_ is a good package to start.

Linux
-----

pyQPCR needs several python libraries as dependencies:

   * `python <http://www.python.org>`_ 2.5 or higher
   * `PyQt4 <http://www.riverbankcomputing.co.uk/software/pyqt/intro>`_ version 4.2 or higher
   * `python-matplotib <http://matplotlib.sourceforge.net/>`_ 0.98 or higher
   * `python-numpy <http://www.numpy.scipy.org>`_ 1.1 or higher
   * `scipy <http://www.scipy.org>`_ 0.6 or higher

Just install them on your favorite distribution. After installing these dependencies, you can install pyQPCR by running::

   sudo python setup.py install

If you work with a RPM-based distribution, you can also build the RPM 
with the following command::

   python setup.py bdist_rpm

and then install it with::

   sudo rpm -ivh dist/pyQPCR-[version]-1.noarch.rpm

t is known to work with Fedora and Opensuse Linux (don't forget to install python-devel and rpmbuild).


MacOS X
-------

Just download the .dmg on the download page and install it by dragging 
the .app file into your /Applications directory. Everything should work 
flawlessly. It is known to work with Mac OS 10.4 or later (Mac Intel only). 

.. warning:: For now, powerpc architectures are not supported.

Run
===

After your successful install, pyQPCR should now be in the menu entry (under
Linux and Windows) or in the Applications (under MacOS X)

You can also run it in the console with the command::

   qpcr

It may be useful for some debugging messages.

.. rubric:: Footnotes

.. [#f1] http://subversion.tigris.org/
