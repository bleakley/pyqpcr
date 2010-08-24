Quick start
***********

Import QPCR raw data / open existing file
=========================================

During this first step, you can:

* **Create a new project**: you give a project name, choose the PCR device (for now Eppendorf, Roche LightCycler, AB 7000 and AB StepOne ones are supported, but others can be easily added) and import your raw data (TXT or CSV files) of one or several plates. Some examples of these files are given with the source of pyQPCR.
* **Open an existing one**: pyQPCR has its own file format which is XML based. You can directly open these files (examples are in the source code of pyQPCR).

.. image:: img/newproj.png
   :scale: 50 %
   :align: center

.. image:: img/addplate.png
   :scale: 50 %
   :align: center

At any time, you can add or remove a plate from your project thanks to the
corresponding icons. It is also possible to merge two different projects.


Plate settings
==============

You can edit the data of each well separately or select and modify a group of
wells. You also can change the targets and samples properties (name, efficiency
of the primers), and remove or add new ones. You can disable wells in order to
not take them into account for calculations.

.. image:: img/settings.png
   :scale: 50 %
   :align: center

Standard curve calculation
==========================

You can define as "standards" the wells that contains dilutions of DNA in order
to calculate PCR efficiency. Then, you precise the amount of DNA (arbitrary
unit) in the different wells and the program will plot the standard curve and
calculate PCR efficiency for this set of primers. This efficiency will be taken
into account for subsequent relative quantifications.

.. image:: img/standard.png
   :scale: 50 %
   :align: center

Reference target and sample
===========================

For relative quantification calculations, you must define one or several
reference gene(s) and one reference sample. They can be either shared for all
plates or specific of each plate.

.. image:: img/refs.png
   :scale: 50 %
   :align: center


Relative quantification
=======================

The wells defined as "unknown" are used to calculate relative quantifications.
An improved :math:`\Delta\Delta Ct` method allows you to obtain reliable
quantifications and error. The confidence level is modifiable and can be either
gaussian or calculated using a T-test.

.. image:: img/rel1.png
   :scale: 50 %
   :align: center

The program plots results as histograms that are easy to customize (colors, 
legend, order, ...)

.. image:: img/rel2.png
   :scale: 50 %
   :align: center

.. image:: img/rel3.png
   :scale: 50 %
   :align: center

Results, export and save
========================

Results can be printed or exported in a pdf file containing a table with all
the data and plots for standard curves and/or relatives quantifications.

.. image:: img/res1.png
   :scale: 50 %
   :align: center

You can also save your project in the pyQPCR XML file format that allows you to
keep the entire project with the different plates and settings easily
recoverable.

Help
====

A Help menu is available and summarize the different functionalities of the
software.

.. image:: img/help.png
   :scale: 50 %
   :align: center
