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

import re
import csv
from pyQPCR.wellGeneSample import Ech, Gene, Puits, WellError
from scipy.stats import t, norm
from PyQt4.QtCore import Qt, QString, QFileInfo
from numpy import mean, std, sqrt, log, log10, polyval, polyfit, sum, \
array, append
from pyQPCR.utils.odict import OrderedDict
from pyQPCR.utils.ragged import RaggedArray2D

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Rev$"


class Plaque:
    
    def __init__(self, filename=None, machine='Eppendorf'):
        self.unsaved = False
        self.filename = filename

        self.listePuits = []
        self.listGene = [Gene('')]
        self.listEch = [Ech('')]
        self.listAmount = ['']

        self.geneRef = ''
        self.echRef = ''
 
        if self.filename is not None:
            self.determineFileType(self.filename)
            if machine == 'Eppendorf':
                self.parseEppendorf()
            elif machine == 'Applied':
                self.parseApplied()
            # Raise exception if no well are detected
            if len(self.listePuits) == 0:
                raise PlateError(self.filename, machine)
# Permet eventuellement de connaitre les genes/ech de ref 
        #self.getRefsFromFile()

    def determineFileType(self, filename):
        extent = filename[-3:]
        if extent in ["txt", "TXT"]:
            self.fileType = "txt"
        elif extent in ["csv", "CSV"]:
            self.fileType = "csv"
        else:
            raise IOError

    def parseEppendorf(self):
        """
        This method allows to parse Eppendorf raw data.
        """
        file = open(self.filename, "r")
        motif = re.compile(r"[\w\s]*")
        amountMotif = re.compile(r"Amount SYBR ?\[(.*)\]")
        if self.fileType == "txt":
            splitter = re.compile(r'("[\w\s.,\-\(\)\[\]\+\\/]*"|\d+[.,]?\d*)')
            iterator = file.readlines()
        if self.fileType == "csv":
            iterator = csv.reader(file, delimiter=";")
        for ind, line in enumerate(iterator):
            if self.fileType == "txt": line = splitter.findall(line)

            if ind == 0:
                self.header = OrderedDict()
                for i, field in enumerate(line):
                    st = field.strip('"')
                    if amountMotif.match(st):
                        self.stdUnit = amountMotif.match(st).group(1)
                        self.header['Amount SYBR'] = i
                    else:
                        self.header[st] = i
                ncol = len(self.header.keys())

            if len(line) == ncol and ind != 0:
                champs = []
                for field in line:
                    try:
                        dat = float(field.replace(',', '.'))
                    except ValueError:
                        dat = field.strip('"')
                    champs.append(dat)
                if self.header.has_key('Pos'):
                    name = champs[self.header['Pos']]
                    x = Puits(name)
                else:
                    raise KeyError
                if self.header.has_key('Name'):
                    echName = champs[self.header['Name']]
                    x.setEch(Ech(echName))
                if self.header.has_key('Ct SYBR'):
                    ct = champs[self.header['Ct SYBR']]
                    x.setCt(ct)
                if self.header.has_key('Ct Mean SYBR'):
                    ctmean = champs[self.header['Ct Mean SYBR']]
                    x.setCtmean(ctmean)
                if self.header.has_key('Ct Dev. SYBR'):
                    ctdev = champs[self.header['Ct Dev. SYBR']]
                    x.setCtdev(ctdev)
                if self.header.has_key('Amount SYBR'):
                    amount = champs[self.header['Amount SYBR']]
                    if amount != '-':
                        x.setAmount(amount)
                    else:
                        x.setAmount('')
                if self.header.has_key('Target SYBR'):
                    geneName = champs[self.header['Target SYBR']]
                    x.setGene(Gene(geneName))
                if self.header.has_key('Type'):
                    type = champs[self.header['Type']]
                    x.setType(type)
                if self.header.has_key('NRQ'):
                    nrq = champs[self.header['NRQ']]
                    x.setNRQ(nrq)
                if self.header.has_key('NRQerror'):
                    nrqerror = champs[self.header['NRQerror']]
                    x.setNRQerror(nrqerror)
 
                setattr(self, x.name, x)
                self.listePuits.append(x)

            if len(line) == 2:
                name =  motif.findall(line[0].strip('"'))[0].replace(' ', '')
                value = line[1]
                try:
                    value = float(value)
                except ValueError:
                    pass
                setattr(self, name, value)
        file.close()

    def parseApplied(self):
        """
        This method allows to parse Applied raw data.
        """
        file = open(self.filename, 'r')
        iterator = file.readlines()
        fileencoding = "utf-8"
        splitter = re.compile(r'([\w .,\-\(\)\[\]\+\\/]*|\d+[.,]?\d*)\t', re.UNICODE)
        result = re.compile(r'\[Results\]')
        motifSample = re.compile(r'Reference Sample = (.*)')
        motifTarget = re.compile(r'Endogenous Control = (.*)')
        hasHeader = False
        for ind, rawline in enumerate(iterator):
            if len(result.findall(rawline)) != 0:
                hasHeader = True
                initTab = ind + 1
                continue
            rawline = rawline.decode(fileencoding)
            line = splitter.findall(rawline)

            if hasHeader:
                if ind == initTab:
                    self.header = OrderedDict()
                    for i, field in enumerate(line):
                        self.header[field] = i
                    ncol = len(self.header.keys())
            if hasHeader:
                if ind != initTab and len(line) == ncol:
                    champs = []
                    for k, field in enumerate(line):
                        try:
                            if self.header.keys()[k] not in ('Sample Name', 'Target Name'):
                                dat = float(field.replace(',', '.'))
                            else:
                                dat = field
                        except ValueError:
                            dat = field
                        champs.append(dat)
                    if self.header.has_key('Well'):
                        name = champs[self.header['Well']]
                        x = Puits(name)
                    else:
                        raise KeyError
                    if self.header.has_key('Sample Name'):
                        echName = champs[self.header['Sample Name']]
                        x.setEch(Ech(echName))
                    if self.header.has_key(u'C\u0442'):
                        ct = champs[self.header[u'C\u0442']]
                        x.setCt(ct)
                    if self.header.has_key(u'C\u0442 Mean'):
                        ctmean = champs[self.header[u'C\u0442 Mean']]
                        x.setCtmean(ctmean)
                    if self.header.has_key(u'C\u0442 SD'):
                        ctdev = champs[self.header[u'C\u0442 SD']]
                        x.setCtdev(ctdev)
                    if self.header.has_key('Quantity'):
                        amount = champs[self.header['Quantity']]
                    if self.header.has_key('Target Name'):
                        geneName = champs[self.header['Target Name']]
                        x.setGene(Gene(geneName))
                    #if self.header.has_key('Task'):
                        #type = champs[self.header['Task']]
                        #x.setType(type)
                    if self.header.has_key(u'\u0394\u0394C\u0442'):
                        nrq = champs[self.header[u'\u0394\u0394C\u0442']]
                        x.setNRQ(nrq)
                    #if self.header.has_key('NRQerror'):
                        #nrqerror = champs[self.header['NRQerror']]
                        #x.setNRQerror(nrqerror)

                    setattr(self, x.name, x)
                    self.listePuits.append(x)
            if motifSample.match(rawline):
                self.echRef = QString(motifSample.findall(rawline)[0])
            if motifTarget.match(rawline):
                self.geneRef = QString(motifTarget.findall(rawline)[0])


        file.close()

    def write(self, filename):
        f = open(filename, 'w')
        self.determineFileType(filename)
        header=['Pos', 'Name', 'Ct SYBR', 'Ct Mean SYBR', 'Ct Dev. SYBR',
                'Amount SYBR', 'Target SYBR', 'Type', 'NRQ', 'NRQerror']
        if hasattr(self, 'stdUnit'):
            header[5] = 'Amount SYBR [%s]' % self.stdUnit
# -----------TXT------------- 
        if self.fileType == "txt":
# Ecriture du header
            for field in header:
                f.write('"%s"'% field)
                f.write("\t")
            f.write("\n")
# Ecriture des puits
            for well in self.listePuits:
                st = well.writePuits()
                f.write(st)
# Ecriture du `Analysis Parameters`
            f.write("\n")
            f.write('"Analysis Parameters"\n')
            if hasattr(self, 'geneRef'):
                f.write('"refTarget"\t"%s"\n' % self.geneRef.name)
            if hasattr(self, 'echRef'):
                f.write('"refSample"\t"%s"\n' % self.echRef.name)
# -----------CSV------------- 
        elif self.fileType == "csv":
            writer = csv.writer(f, quoting=csv.QUOTE_NONNUMERIC, delimiter=";")
# Ecriture du header
            writer.writerow(header)
# Ecriture des puits
            for well in self.listePuits:
                writer.writerow([well.name, well.ech.name, well.ct, well.ctmean,
                    well.ctdev, well.amount, well.gene.name, well.type, well.NRQ,
                    well.NRQerror])
# Ecriture du `Analysis Parameters`
            writer.writerow("")
            writer.writerow(['Analysis Parameters'])
            if hasattr(self, 'geneRef'):
                writer.writerow(['refTarget', self.geneRef.name])
            if hasattr(self, 'echRef'):
                writer.writerow(['refSample', self.echRef.name])
        f.close()

    def writeHtml(self):
        html = u""
        html += "<table cellpadding=2 cellspacing=0 border=1 width=100%>\n"
        html += ("<tr>\n"
                 "<th align=center>Well</th>\n"
                 "<th align=center>Type</th>\n"
                 "<th align=center>Target</th>\n"
                 "<th align=center>Sample</th>\n"
                 "<th align=center>Ct</th>\n"
                 "<th align=center>Ctmean</th>\n"
                 "<th align=center>Ctdev</th>\n"
                 "<th align=center>Amount</th>\n"
                 "<th align=center>Efficiency</th>\n"
                 "<th align=center>NRQ</th>\n"
                 "<th align=center>NRQerror</th>\n"
                 "</tr>\n")
        for well in self.listePuits:
            html += well.writeHtml()
        html += "</table>"
        return html

    def setDicoGene(self):
        self.dicoGene = OrderedDict()
        for well in self.listePuits:
            nomgene = well.gene.name
            if nomgene != '':
                if self.dicoGene.has_key(nomgene):
                    self.dicoGene[nomgene].append(well)
                else:
                    self.dicoGene[nomgene] = [well]

    def setDicoEch(self):
        self.dicoEch = OrderedDict()
        for well in self.listePuits:
            nomech = well.ech.name
            if nomech != '':
                if self.dicoEch.has_key(nomech):
                    self.dicoEch[nomech].append(well)
                else:
                    self.dicoEch[nomech] = [well]

    def getRefsFromFile(self):
        """
        This methods allows to determine both the Gene of Interest and
        the reference-sample thanks to the input data file
        """
        if hasattr(self, 'refTarget'):
            goi = self.refTarget.strip('"')
            ind = self.adresseGene[QString(goi)]
            self.geneRef = self.listGene[ind]
            self.listGene[ind].setRef(Qt.Checked)
        if hasattr(self, 'refSample'):
            soi = self.refSample.strip('"')
            ind = self.adresseEch[QString(soi)]
            self.echRef = self.listEch[ind]
            self.listEch[ind].setRef(Qt.Checked)


class Replicate:

    def __init__(self, listePuits, type=QString('unknown'), 
                 confidence=0.9, errtype="normal"):
        self.confidence = confidence
        self.errtype = errtype
        self.type = type
        self.listePuits = listePuits
        if len(self.listePuits) != 0:
            self.gene = self.listePuits[0].gene
        else:
            self.gene = Gene('')
        self.ctList =array([])
        for well in self.listePuits:
            self.ctList = append(self.ctList, well.ct)
        if self.type == QString('unknown'):
            self.ech = self.listePuits[0].ech
        elif self.type == QString('standard'):
            self.amList = array([])
            for well in self.listePuits:
                self.amList = append(self.amList, well.amount)
        self.calcMeanDev()

    def __str__(self):
        st = '{%s:[' % self.type
        for well in self.listePuits:
            st = st + well.name + ','
        st += ']'
        if self.type == QString('unknown'):
            st += ' %s, %s}' % (self.gene, self.ech)
        else:
            st += ' %s}' % self.gene
        return st

    def __repr__(self):
        st = '['
        for well in self.listePuits:
            st = st + well.name + ','
        st += ']'
        return st

    def setNRQ(self, NRQ):
        self.NRQ = NRQ

    def setNRQerror(self, NRQerr):
        self.NRQerror = NRQerr

    def calcMeanDev(self):
        """
        Compute the mean ct of a replicate
        Formule 7
        """
        try:
            self.ctmean = self.ctList.mean()
        except TypeError:
            brokenWells = []
            for well in self.listePuits:
                try:
                    f = float(well.ct)
                except ValueError:
                    brokenWells.append(well.name)
                    well.setWarning(True)
            raise WellError(brokenWells)

        if len(self.ctList) > 1:
            # Formule 8
            stdctList = self.ctList.std()*sqrt(1./(len(self.ctList)-1.))
            if self.errtype == "student":
                # coeff Student
                talpha = t.ppf(1.-(1.-self.confidence)/2., len(self.ctList)-1) 
            elif self.errtype == "normal":
                talpha = norm.ppf(1.-(1.-self.confidence)/2.) # Gaussian
            self.ctdev = stdctList*sqrt(len(self.ctList)-1.)
            self.ctdevtalpha = talpha * stdctList
        else:
            self.ctdev = 0.
            self.ctdevalpha = 0.

        for well in self.listePuits:
            well.setCtmean(self.ctmean)
            well.setCtdev(self.ctdev)

    def calcDCt(self):
        self.dct = self.gene.ctref - self.ctmean # Formule 10
        self.RQ = (1.+self.gene.eff/100.)**(self.dct) # Formule 11

    def calcRQerror(self):
        # Formule 12
        err = sqrt( self.RQ**2 * ((self.dct*(self.gene.pm/100.) \
                /(1.+self.gene.eff/100.))**2 \
                + (log(1.+self.gene.eff/100.)*self.ctdevtalpha)**2 \
                ))
        self.RQerror = err


class ReplicateError(Exception):

    def __init__(self, listRep):
        self.listRep = listRep

    def __str__(self):
        st = "<ul>"
        for trip in self.listRep:
            st += "<li>(<b>%s, %s</b>) : E(ct)=%.2f </li>" % (trip.gene, trip.ech, trip.ctdev)
        st += "</ul>"

        return st

class PlateError(Exception):

    def __init__(self, filename, machine):
        self.filename = filename
        self.machine = machine

    def __str__(self):
        st = "<b>Warning</b> : The file <b>%s </b> does not contain any well at the right format. " %  \
              QFileInfo(self.filename).fileName()
        st += "It probably comes from your raw data file. Your current PCR device is"
        st += " <b>%s</b>, check your file corresponds to this machine !" % self.machine
        st += " If the error continues to occur, post a message at "
        st += r' <a href="http://sourceforge.net/projects/pyqpcr/forums/forum/918935">'
        st += r'http://sourceforge.net/projects/pyqpcr/forums/forum/918935</a>'
        return st


if __name__ == '__main__':
    pl = Plaque('../samples/raw_std.txt')
    print str(pl.A1.NRQ)
