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
    This class defines the sample.

    @ivar name: the sample name
    @type name: PyQt4.QtCore.QString
    @ivar isRef: a boolean to determine if the sample is a reference
                 sample
    @type isRef: PyQt4.QtCore.CheckState
    @ivar enabled: a boolean to determine if the sample is enabled or
                   disabled
    @type enabled: PyQt4.QtCore.CheckState
    @ivar color: the sample color (for plotting purpose)
    @type color: PyQt4.QtGui.QColor
    """

    def __init__(self, nom, isRef=Qt.Unchecked):
        """
        Ech constructor:

            >>> Ech('A1', isRef=Qt.Unchecked)

        @param nom: the sample name
        @type nom: PyQt4.QtCore.QString
        @param isRef: a boolean to determine if the sample is enabled or
                      disabled
        @type isRef: PyQt4.QtCore.CheckState
        """
        self.name = QString(nom)
        self.isRef = isRef
        self.enabled = Qt.Checked

    def __str__(self):
        st =  "%s" % self.name
        return unicode(st)

    def __repr__(self):
        st =  "%s" % self.name
        return st

    def setName(self, name):
        """
        Set the sample name

        @param name: the new name of the sample
        @type name: PyQt4.QtCore.QString
        """
        self.name = name

    def setRef(self, checkBoxState):
        """
        Set if the sample is a reference sample or not.

        @param checkBoxState: a boolean which indicates wheter the sample
                              is a reference sample or not.
        @type checkBoxState: PyQt4.QtCore.CheckState
        """
        self.isRef = checkBoxState

    def setColor(self, color):
        """
        Set the sample color in the plot

        @param color: the sample color
        @type color: PyQt4.QtGui.QColor
        """
        self.color = color

    def setEnabled(self, ena):
        """
        enable/disable the current sample

        @param ena: a boolean which indicates whether the sample
                    is enabled or disabled
        @type ena: PyQt4.QtCore.CheckState
        """
        self.enabled = ena

class Gene:

    def __init__(self, nom, eff=100., pm=0., isRef=Qt.Unchecked):
        self.name = QString(nom)
        self.eff = eff
        self.pm = pm
# Pour dire si le gene est un gene de reference
        self.isRef = isRef
        self.ctref = nan
# Pour dire si on veut tracer un gene
        self.enabled = Qt.Checked

    def __str__(self):
        #st =  "%s (%.2f %% +/- %.2f)" % (self.name, self.eff, self.pm) 
        st =  "%s" % self.name
        return unicode(st)

    def __repr__(self):
        st =  "%s (%.2f %% +/- %.2f)" % (self.name, self.eff, self.pm) 
        return st

    def setPos(self, pos):
        self.pos = pos

    def setName(self, name):
        """
        Set the target name

        @param name: the new name of the target
        @type name: PyQt4.QtCore.QString
        """
        self.name = name

    def setPm(self, pm):
        """
        Set the target efficiency relative-error

        @param pm: the target efficiency relative-error
        @type pm: float
        """
        self.pm = pm

    def setEff(self, eff):
        """
        Set the target efficiency

        @param eff: the target efficiency
        @type eff: float
        """
        self.eff = eff

    def setRef(self, checkBoxState):
        """
        Set if the target is a reference sample or not.

        @param checkBoxState: a boolean which indicates wheter the target
                              is a reference one or not.
        @type checkBoxState: PyQt4.QtCore.CheckState
        """
        self.isRef = checkBoxState

    def setColor(self, color):
        """
        set the target color in the plot

        @param color: the target color
        @type color: PyQt4.QtGui.QColor
        """
        self.color = color

    def setEnabled(self, ena):
        """
        enable/disable the current target

        @param ena: a boolean which indicates whether the target
                    is enabled or disabled
        @type ena: PyQt4.QtCore.CheckState
        """
        self.enabled = ena

    def setR2(self, R):
        self.R2 = R

    def calcCtRef(self, listePuits):
        """
        This methods calculates the mean ct of every wells of a given
        target.

        @param listePuits: the list of wells of the current gene
        @type listePuits: pyQPCR.wellGeneSample.Puits
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

    def __init__(self, name, ech=QString(''), ct=nan, ctmean=nan, 
            ctdev=nan, gene=QString('')):
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
        st = '\nPuit ' + self.name + "\n" + '(' + str(self.ech) + ', ' + \
              str(self.gene) + ')' + "\n" + \
              "ct = %.2f, ctmean = %.2f, ctdev = %.2f"%( \
                self.ct, self.ctmean, self.ctdev)
        return st

    def __repr__(self):
        st = "%s: %s, %s, %s" % (self.name, self.gene, self.ech, self.type)
        return st

    def writePuits(self, html=False):
        if str(self.amount) == '':
            amount = '"-"'
        else:
            amount = str(self.amount)
        if str(self.ct) == '':
            ct = '""'
        else:
            ct = str(self.ct)
        if str(self.ctmean) == '':
            ctmean = '""'
        else:
            ctmean = "%.2f" % self.ctmean
        if str(self.ctdev) == '':
            ctdev = '""'
        else:
            ctdev = "%.2f" % self.ctdev
        if str(self.NRQ) == '':
            NRQ = '""'
        else:
            NRQ = "%.2f" % self.NRQ
        if str(self.NRQ) == '':
            NRQerror = '""'
        else:
            NRQerror = "%.2f" % self.NRQerror

        st = '"%s"\t"%s"\t%s\t%s\t%s\t%s\t"%s"\t"%s"\t%s\t%s\n' % \
            (self.name, str(self.ech), ct, ctmean, 
             ctdev, amount, self.gene.name, self.type,
             NRQ, NRQerror)
        return st

    def writeWellXml(self):
        amount = str(self.amount)
        if str(self.ct) != '':
            try:
                ct = "%.2f" % self.ct
            except TypeError:
                ct = str(self.ct)
        else:
            ct = ''
        if str(self.ctmean) != '':
            ctmean = "%.2f" % self.ctmean
        else:
            ctmean = ''
        if str(self.ctdev) != '':
            ctdev = "%.2f" % self.ctdev
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
        motif = re.compile(r"([A-H])(\d+)")
        groups = motif.search(self.name).groups()
        self.ypos = int(groups[1])-1
        lettres = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        dict = {}
        for xpos, l in enumerate(lettres):
            dict[l] = xpos
        self.xpos = dict[self.name[0]]

    def setGene(self, gene):
        self.gene = gene

    def setEch(self, ech):
        self.ech = ech

    def setType(self, name):
        if self.type == 'standard' and name == 'unknown':
            self.amount = ''
        self.type = QString(name)

    def setAmount(self, qte):
        if self.type != 'unknown':
            self.amount = qte
        else:
            self.amount = ''

    def setCt(self, ct):
        self.ct = ct

    def setCtmean(self, ctmean):
        self.ctmean = ctmean

    def setCtdev(self, ctdev):
        self.ctdev = ctdev

    def setEnabled(self, ena):
        self.enabled = ena

    def setNRQ(self, nrq):
        try:
            self.NRQ = float(nrq)
        except ValueError:
            self.NRQ = nrq

    def setNRQerror(self, err):
        self.NRQerror = err

    def setWarning(self, warn):
        """Si un puit est casse, on met un flag warning dessus"""
        self.warning = warn


class WellError(Exception):

    def __init__(self, brokenWells):
        self.brokenWells = brokenWells

    def __str__(self):
        return repr(self.brokenWells)



if __name__=="__main__":
    a = 'toto'
    eff = 90
    pm = 0.1
    g = Gene(a, eff, pm)
    ech = Ech('si')

