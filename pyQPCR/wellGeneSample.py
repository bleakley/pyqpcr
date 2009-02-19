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

from PyQt4.QtCore import *
from PyQt4.QtGui import *
from numpy import nan
import re

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Rev$"

class Ech:

    def __init__(self, nom, isRef=Qt.Unchecked):
        self.name = str(nom)
# Pour dire si l'echantillon est un echantillon de reference
        self.isRef = isRef
# Pour dire si on veut tracer un echantillon
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
        @type name: string
        """
        self.name = name

    def setRef(self, checkBoxState):
        self.isRef = checkBoxState

    def setColor(self, color):
        """
        Set the sample color in the plot

        @param color: the sample color
        @type color: PyQt4.QtGui.QColor
        """
        self.color = color

    def setEnabled(self, ena):
        self.enabled = ena

class Gene:

    def __init__(self, nom, eff=100., pm=0., isRef=Qt.Unchecked):
        self.name = str(nom)
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
        self.name = name

    def setPm(self, pm):
        self.pm = pm

    def setEff(self, eff):
        self.eff = eff

    def setRef(self, checkBoxState):
        self.isRef = checkBoxState

    def setColor(self, color):
        """
        set the target color in the plot

        @param color: the target color
        @type color: PyQt4.QtGui.QColor
        """
        self.color = color

    def setEnabled(self, ena):
        self.enabled = ena

    def calcCtRef(self, listePuits):
        qt = 0
        k = 0
        for well in listePuits:
            qt += well.ct
            k += 1
        self.ctref = float(qt/k)

class Puits:

    def __init__(self, name, ech=None, ct=nan, ctmean=nan, 
            ctdev=nan, gene=None):
        self.name = name
        self.ech = Ech(ech)
        self.gene = Gene(gene)
        self.ct = ct
        self.ctmean = ctmean
        self.ctdev = ctdev
        self.amount = ''
        self.type = "unknown"
        self.NRQ = ''
        self.NRQerror = ''
        self.getPosition()
        self.enabled = True

    def __str__(self):
        st = '\nPuit ' + str(self.name) + "\n" + '(' + str(self.ech) + ', ' + \
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
            (str(self.name), str(self.ech), ct, ctmean, 
             ctdev, amount, str(self.gene.name), str(self.type),
             NRQ, NRQerror)
        return st

    def writeHtml(self):
        amount = str(self.amount)
        try:
            ct = "%.2f" % self.ct
        except TypeError:
            ct = ''
        try:
            ctmean = "%.2f" % self.ctmean
        except TypeError:
            ctmean = ''
        try:
            ctdev = "%.2f" % self.ctdev
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

        st = ("<tr><td align=center><b>%s</b></td>\n" # name
              "<td align=center>%s</td>\n" # type
              "<td align=center>%s</td>\n" # gene
              "<td align=center>%s</td>\n" # sample
              "<td align=center>%s</td>\n" # ct
              "<td align=center>%s</td>\n" # ctmean
              "<td align=center>%s</td>\n" # ctdev
              "<td align=center>%s</td>\n" # amount
              "<td align=center>%s</td>\n" # NRQ
              "<td align=center>%s</td></tr>\n") % (str(self.name), str(self.type), 
                    str(self.gene.name), str(self.ech.name), ct, 
                    ctmean, ctdev, amount, NRQ, NRQerror)
        return st

    def getPosition(self):
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

    def setEch(self, name):
        self.ech = name

    def setType(self, name):
        self.type = name

    def setAmount(self, qte):
        self.amount = qte

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


if __name__=="__main__":
    a = 'toto'
    eff = 90
    pm = 0.1
    g = Gene(a, eff, pm)
    print g
    ech = Ech('si')
    print ech

