# -*- coding: utf-8 -*-
#
# pyQPCR, an application to analyse qPCR raw data
# Copyright (C) 2008 Thomas Gastine
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

from PyQt4.QtCore import Qt, QString
from numpy import nan
import re

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Rev$"

class Ech:
    """
    An Ech object is used to define a sample. This object is constructed
    from the name of the sample.

    >>> Ech('A1', isRef=Qt.Unchecked)

    :attribute name: the sample name
    :attribute isRef: a boolean to determine if the sample is a reference
                      sample
    :attribute enabled: a boolean to determine if the sample is enabled or
                        disabled
    :attribute color: the sample color (for plotting purpose)
    """

    def __init__(self, nom, isRef=Qt.Unchecked):
        """
        Constructor of Ech object

        :param nom: the sample name
        :type nom: PyQt4.QtCore.QString
        :param isRef: a boolean to determine if the sample is enabled or
                      disabled
        :type isRef: PyQt4.QtCore.CheckState
        """
        self.name = QString(nom)
        self.isRef = isRef
        self.enabled = Qt.Checked

    def __str__(self):
        """
        Print function
        """
        st =  "%s" % self.name
        return unicode(st)

    def __repr__(self):
        """
        Print function
        """
        st =  "%s" % self.name
        return st

    def setName(self, name):
        """
        Set the sample name

        :param name: the new name of the sample
        :type name: PyQt4.QtCore.QString
        """
        self.name = name

    def setRef(self, checkBoxState):
        """
        Set if the sample is a reference sample or not.

        :param checkBoxState: a boolean which indicates wheter the sample
                              is a reference sample or not.
        :type checkBoxState: PyQt4.QtCore.CheckState
        """
        self.isRef = checkBoxState

    def setColor(self, color):
        """
        Set the sample color in the plot

        :param color: the sample color
        :type color: PyQt4.QtGui.QColor
        """
        self.color = color

    def setEnabled(self, ena):
        """
        enable/disable the current sample

        :param ena: a boolean which indicates whether the sample
                    is enabled or disabled
        :type ena: PyQt4.QtCore.CheckState
        """
        self.enabled = ena

class Gene:
    """
    A Gene object is used to define a target. This object is constructed
    from the name of the target and its efficiency.

    >>> g = Gene('foo', eff=95, pm=1)

    :attribute name: the name of the current target
    :attribute eff: the efficiency of the target
    :attribute pm: the standard error on efficiency
    :attribute isRef: a boolean to indicate if the target is a reference one
    :attribute enabled: a boolean to indicate if the current Gene is
                        enabled or not (useful for plots)

    """

    def __init__(self, nom, eff=100., pm=0., isRef=Qt.Unchecked):
        """
        Constructor of Gene object

        :param nom: the name of the target
        :type nom: PyQt4.QtCore.QString
        :param eff: the efficiency of the target
        :type eff: float
        :param pm: the standard error on the efficiency
        :type eff: float
        :param isRef: a boolean to indicate if the gene is a reference one
        :type isRef: PyQt4.QtCore.CheckState
        """
        self.name = QString(nom)
        self.eff = eff
        self.pm = pm
        # Pour dire si le gene est un gene de reference
        self.isRef = isRef
        self.ctref = nan
        # Pour dire si on veut tracer un gene
        self.enabled = Qt.Checked

    def __str__(self):
        """
        Printing method
        """
        st =  "%s" % self.name
        return unicode(st)

    def __repr__(self):
        """
        Printing method
        """
        st =  "%s (%.2f %% +/- %.2f)" % (self.name, self.eff, self.pm) 
        return st

    def setName(self, name):
        """
        Set the target name

        :param name: the new name of the target
        :type name: PyQt4.QtCore.QString
        """
        self.name = name

    def setPm(self, pm):
        """
        Set the target efficiency standard error

        :param pm: the standard error of the efficiency
        :type pm: float
        """
        self.pm = pm

    def setEff(self, eff):
        """
        Set the target efficiency

        :param eff: the target efficiency
        :type eff: float
        """
        self.eff = eff

    def setRef(self, checkBoxState):
        """
        Set if the target is a reference sample or not.

        :param checkBoxState: a boolean which indicates wheter the target
                              is a reference one or not.
        :type checkBoxState: PyQt4.QtCore.CheckState
        """
        self.isRef = checkBoxState

    def setColor(self, color):
        """
        Set the color of the target color (for the plot of quantification)

        :param color: the target color
        :type color: PyQt4.QtGui.QColor
        """
        self.color = color

    def setEnabled(self, ena):
        """
        Enable/disable the current target

        :param ena: a boolean which indicates whether the target
                    is enabled or disabled
        :type ena: PyQt4.QtCore.CheckState
        """
        self.enabled = ena

    def setR2(self, R):
        """
        Set the value of R2 for the current target

        :param R: the value of R2
        :type R: float
        """
        self.R2 = R

    def calcCtRef(self, listePuits):
        """
        This methods calculates the mean ct of every wells of a given
        target.

        .. math:: {c_t} = \dfrac{1}{n_w}\sum_{w} {c_t}_w

        :param listePuits: the list of wells of the current gene
        :type listePuits: pyQPCR.wellGeneSample.Puits
        """
        qt = 0
        k = 0
        brokenWells = []
        for well in listePuits:
            try:
                if well.enabled and well.type == 'unknown':
                    qt += well.ct
                    k += 1
            except TypeError:
                well.setWarning(True)
                brokenWells.append(well.name)
                continue
        if len(brokenWells) != 0:
            raise WellError(brokenWells)
        elif k == 0:
            self.ctref = 0
        else:
            self.ctref = float(qt/k)

class Puits:
    """
    The Puits object is used to define a well. It is constructed from the
    name of the well, its type and the value of the ct.

    >>> a = Puits('A1', ct=23.4, ech='ech', gene='gene', type='unknown')

    :attribute name: the name of the well
    :attribute ech: the Ech object of the well
    :attribute gene: the Gene object of the well
    :attribute ct: the value of the ct
    :attribute ctmean: the mean value of ct in a replicate
    :attribute ctdev: the standard deviation of ct in a replicate
    :attribute amount: the quantity for a standard-type well
    :attribute NRQ: the value of NRQ after quantification
    :attribute NRQerror: the standard-error on NRQ
    :attribute enabled: a boolean to indicate if the well is enabled or not
    :attribute warning: a boolean to indicate if there is a problem
                        with this well
    """

    def __init__(self, name, ech=QString(''), ct=nan, ctmean=nan, 
            ctdev=nan, gene=QString('')):
        """
        Constructor of the Puits object

        :param name: the name of the well
        :type name: PyQt4.QtCore.QString
        :param ech: the name of the sample of the well
        :type ech: PyQt4.QtCore.QString
        :param ct: the value of ct
        :type ct: float
        :param ctmean: the mean value of ct in a replicate
        :type ctmean: float
        :param ctdev: the standard deviation of ct in a replicate
        :type ctdev: float
        :param gene: the name of the target of the well
        :type gene: PyQt4.QtCore.QString
        """
        self.name = name
        self.ech = Ech(ech)
        self.gene = Gene(gene)
        self.ct = ct
        self.ctmean = ctmean
        self.ctdev = ctdev
        self.amount = ''
        self.type = QString("unknown")
        self.NRQ = ''
        self.NRQerror = ''
        self.getPosition()
        self.enabled = True
        self.warning = False

    def __str__(self):
        """
        Print method
        """
        st = '\nPuit ' + self.name + "\n" + '(' + str(self.ech) + ', ' + \
              str(self.gene) + ')' + "\n" + \
              "ct = %.2f, ctmean = %.2f, ctdev = %.2f"%( \
                self.ct, self.ctmean, self.ctdev)
        return st

    def __repr__(self):
        """
        Print method
        """
        st = "%s: %s, %s, %s" % (self.name, self.gene, self.ech, self.type)
        return st

    def writeWellXml(self):
        """
        This method is used to represent the Well object (with its
        attribute) in a XML string. This is used to save the project
        in an XML file.
        """
        amount = str(self.amount)
        if str(self.ct) != '':
            try:
                ct = "%.2f" % self.ct
            except TypeError:
                ct = str(self.ct)
        else:
            ct = ''
        if str(self.ctmean) != '':
            try:
                ctmean = "%.2f" % self.ctmean
            except TypeError:
                ctmean = str(self.ctmean)
        else:
            ctmean = ''
        if str(self.ctdev) != '':
            try:
                ctdev = "%.2f" % self.ctdev
            except TypeError:
                ctdev = str(self.ctdev)
        else:
            ctdev = ''
        if str(self.NRQ) != '':
            NRQ = "%.2f" % self.NRQ
            NRQerror = "%.2f" % self.NRQerror
        else:
            NRQ = ''
            NRQerror = ''

        st ="<WELL CT='%s' CTMEAN='%s' CTDEV='%s' " % (ct, ctmean, ctdev)
        st += "AMOUNT='%s' NRQ='%s' NRQERROR='%s' " % (amount, NRQ, NRQerror)
        st += "ENABLED='%i' >\n" % int(self.enabled)
        st += "<NAME>%s</NAME>\n" % self.name
        st += "<TYPE>%s</TYPE>\n" % self.type
        st += "<TARGET EFF='%.2f' PM='%.2f'>%s</TARGET>\n" % \
            (self.gene.eff, self.gene.pm, self.gene.name)
        st += "<SAMPLE>%s</SAMPLE>\n" % self.ech.name
        st += "</WELL>\n"
        return st

    def writeHtml(self, ctMin=35, ectMax=0.3):
        """
        This method is used to display a Puits object in an HTML string.
        It is called when one wants print the result in a PDF file.

        :param ctMin: the minimum value for a ct (negative)
        :type ctMin: float
        :param ectMax: the maximum value for the standard deviation of the
                       ct in a replicate.
        :type ectMax: float
        """
        try:
            amount = "%.2f" % self.amount
        except TypeError:
            amount =  str(self.amount)
        try:
            if self.ct <= ctMin and self.type == 'negative':
                ct = "<img src=':/flag.png' width=8> "
            else:
                ct = ''
            ct += "%.2f" % self.ct
        except TypeError:
            ct = ''
        try:
            ctmean = "%.2f" % self.ctmean
        except TypeError:
            ctmean = ''
        try:
            if self.ctdev >= ectMax:
                ctdev = "<img src=':/flag.png' width=8> "
            else:
                ctdev = ''
            ctdev += "%.2f" % self.ctdev
        except TypeError:
            ctdev = ''
        try:
            NRQ = "%.2f" % self.NRQ
        except TypeError:
            NRQ = ''
        try:
            NRQerror = "%.2f" % self.NRQerror
        except TypeError:
            NRQerror = ''
        eff = "%.2f%s%.2f" % (self.gene.eff, unichr(177), self.gene.pm)

        if self.enabled:
            name = "<img src=':/enable.png' width=8> <b>%s</b>" % self.name
        else:
            name = "<img src=':/disable.png' width=8> <b>%s</b>" % self.name
        if self.type == 'unknown':
            bgcolor = '#e6e6fa'
        elif self.type == 'standard':
            bgcolor = '#ffe4e1'
        elif self.type == 'negative':
            bgcolor = '#fff8d6'

        st = ("<tr><td align=center><b>%s</b></td>\n" # name
              "<td bgcolor=%s align=center>%s</td>\n" # type
              "<td align=center>%s</td>\n" # gene
              "<td align=center>%s</td>\n" # sample
              "<td align=center>%s</td>\n" # ct
              "<td align=center>%s</td>\n" # ctmean
              "<td align=center>%s</td>\n" # ctdev
              "<td align=center>%s</td>\n" # amount
              "<td align=center>%s</td>\n" # eff
              "<td align=center>%s</td>\n" # NRQ
              "<td align=center>%s</td></tr>\n") % (name, bgcolor,
                    self.type, self.gene.name, self.ech.name, ct, 
                    ctmean, ctdev, amount, eff, NRQ, NRQerror)
        return st

    def getPosition(self):
        """
        This method gives the well position (x,y) of a given well
        according to its name.

        ex. A11 xpos=0 ypos=10
        """
        motif = re.compile(r"([A-P])(\d+)")
        groups = motif.search(self.name).groups()
        self.ypos = int(groups[1])-1
        lettres = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
                   'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
        dict = {}
        for xpos, l in enumerate(lettres):
            dict[l] = xpos
        self.xpos = dict[self.name[0]]

    def setGene(self, gene):
        """
        Set the target of the well

        :param gene: the target (object)
        :type gene: Gene
        """
        self.gene = gene

    def setEch(self, ech):
        """
        Set the sample of the well

        :param ech: the sample (object)
        :type ech: Ech
        """
        self.ech = ech

    def setType(self, name):
        """
        Set the type of the well

        :param name: the type of the well (standard, unknown or negative)
        :type name: PyQt4.QtCore.QString
        """
        if self.type == 'standard' and name == 'unknown':
            self.amount = ''
        self.type = QString(name)

    def setAmount(self, qte):
        """
        Set the amount of standard-type well

        :param qte: the value of the amount
        :type qte: float
        """
        if self.type != 'unknown':
            self.amount = qte
        else:
            self.amount = ''

    def setCt(self, ct):
        """
        Set the ct value

        :param ct: the value of ct
        :type ct: float
        """
        self.ct = ct

    def setCtmean(self, ctmean):
        """
        Set the ct mean value (in a replicate)

        :param ctmean: the value of ctmean
        :type ctmean: float
        """
        self.ctmean = ctmean

    def setCtdev(self, ctdev):
        """
        Set the standard deviation of ct (in a replicate)

        :param ctdev: the value of ctdev
        :type ctdev: float
        """
        self.ctdev = ctdev

    def setEnabled(self, ena):
        """
        Enable/disable the current well

        :param ena: a boolean
        :type ena: PyQt4.QtCore.CheckState
        """
        self.enabled = ena

    def setNRQ(self, nrq):
        """
        Set the quantification NRQ of the well

        :param nrq: the quantification value
        :type nrq: float
        """
        try:
            self.NRQ = float(nrq)
        except ValueError:
            self.NRQ = nrq

    def setNRQerror(self, err):
        """
        Set the standard error of the quantification of the well

        :param err: the standard error of the quantification
        :type err: float
        """
        self.NRQerror = err

    def setWarning(self, warn):
        """
        Put a warning flag on a broken well

        :param warn: a boolean to indicate if a well is broken or not
        :type warn: PyQt4.QtCore.CheckState
        """
        self.warning = warn


class WellError(Exception):
    """
    This exception is raised if some wells on a plate have a warning state
    """

    def __init__(self, brokenWells):
        """
        Constructor of WellError

        :param brokenWells: a list of the brokenWells
        :type brokenWells: list
        """
        self.brokenWells = brokenWells

    def __str__(self):
        """
        Print error
        """
        return repr(self.brokenWells)



if __name__=="__main__":
    a = 'toto'
    eff = 90
    pm = 0.1
    g = Gene(a, eff, pm)
    ech = Ech('si')

