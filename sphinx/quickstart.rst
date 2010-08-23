Quick start
***********

Import QPCR raw data / open existing file
=========================================

During this first step, you can:

* **Create a new project**: you give a project name, choose the PCR device 
(for now Eppendorf, Roche LightCycler, AB 7000 and AB StepOne ones are
supported, but others can be easily added) and import your raw data (TXT or CSV
files) of one or several plates. Some examples of these files are given with
the source of pyQPCR.
* **Open an existing one**: pyQPCR has its own file format which is XML based. 
You can directly open these files (examples are in the source code of pyQPCR).

.. image:: img/newproj.png
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
