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
from pyQPCR.wellGeneSample import Ech, Gene, Puits
from scipy.stats import t
from PyQt4.QtCore import Qt
from PyQt4.QtGui import *
from numpy import mean, std, sqrt, log, log10, polyval, polyfit, sum, \
array, append
from pyQPCR.utils.odict import OrderedDict
from pyQPCR.utils.ragged import RaggedArray2D

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Rev$"

class Plaque:
    
    def __init__(self, filename):
        self.filename = filename

        self.listePuits = []
        self.listGene = [Gene('')]
        self.listEch = [Ech('')]
        self.listAmount = ['']
        self.adresseGene = OrderedDict()
        self.adresseEch = OrderedDict()
        self.adresseAmount = OrderedDict()
        self.adresseEch[''] = 0
        self.adresseGene[''] = 0
        self.adresseAmount[''] = 0
 
        self.determineFileType(self.filename)
        self.read()
        self.setDicoGene()
        self.setDicoEch()
        self.setDicoAm()
# Permet eventuellement de connaitre les genes/ech de ref 
        self.getRefsFromFile()

    def determineFileType(self, filename):
        extent = filename[-3:]
        if extent in ["txt", "TXT"]:
            self.fileType = "txt"
        elif extent in ["csv", "CSV"]:
            self.fileType = "csv"
        else:
            raise IOError

    def read(self):
        file = open(self.filename, "r")
        motif = re.compile(r"[\w\s]*")
        amountMotif = re.compile(r"Amount SYBR ?\[(.*)\]")
        if self.fileType == "txt":
            splitter = re.compile(r'("[\w\s.,\-\(\)\[\]\+]*"|\d+[.,]?\d*)')
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
                    if not self.adresseEch.has_key(echName):
                        self.listEch.append(Ech(echName))
                        self.adresseEch[echName] = len(self.listEch)-1
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
                        if not self.adresseAmount.has_key(str(amount)):
                            self.listAmount.append(str(amount))
                            self.adresseAmount[str(amount)] = len(self.listAmount)-1
                    else:
                        x.setAmount('')
                if self.header.has_key('Target SYBR'):
                    geneName = champs[self.header['Target SYBR']]
                    if not self.adresseGene.has_key(geneName):
                        self.listGene.append(Gene(geneName))
                        self.adresseGene[geneName] = len(self.listGene)-1 
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
#
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
                 "<th align=center>NRQ</th>\n"
                 "<th align=center>NRQerror</th>\n"
                 "</tr>\n")
        for well in self.listePuits:
            html += well.writeHtml()
        html += "</table>"
        return html

    def setDicoGene(self):
# Mise a jour de dicoGene
        self.dicoGene = OrderedDict()
        for well in self.listePuits:
            nomgene = well.gene.name
            if self.dicoGene.has_key(nomgene):
                self.dicoGene[nomgene].append(well)
            else:
                self.dicoGene[nomgene] = [well]
# Mise a jour de adresseGene
        self.adresseGene = OrderedDict()
        self.adresseGene[''] = 0
        ind = 1
        for gene in self.listGene:
            if not self.adresseGene.has_key(str(gene)):
                self.adresseGene[str(gene)] = ind
                ind += 1

    def setDicoEch(self):
# Mise a jour de dicoEch
        self.dicoEch = OrderedDict()
        for well in self.listePuits:
            nomech = well.ech.name
            if self.dicoEch.has_key(str(nomech)):
                self.dicoEch[str(nomech)].append(well)
            else:
                self.dicoEch[str(nomech)] = [well]
# Mise a jour de adresseEch
        self.adresseEch = OrderedDict()
        self.adresseEch[''] = 0
        ind = 1
        for ech in self.listEch:
            if not self.adresseEch.has_key(str(ech.name)):
                self.adresseEch[str(ech.name)] = ind
                ind += 1

    def setDicoAm(self):
# Mise a jour de dicoAmount
        self.dicoAmount = OrderedDict()
        for well in self.listePuits:
            if self.dicoAmount.has_key(str(well.amount)):
                self.dicoAmount[str(well.amount)].append(well)
            else:
                self.dicoAmount[str(well.amount)] = [well]
# Mise a jour de adresseAmount
        self.adresseAmount = OrderedDict()
        self.adresseAmount[''] = 0
        ind = 1
        for am in self.listAmount:
            if not self.adresseAmount.has_key(str(am)):
                self.adresseAmount[str(am)] = ind
                ind += 1

    def getRefsFromFile(self):
        """
        This methods allows to determine both the Gene of Interest and
        the reference-sample thanks to the input data file
        """
        if hasattr(self, 'refTarget'):
            goi = self.refTarget.strip('"')
            ind = self.adresseGene[str(goi)]
            self.geneRef = self.listGene[ind]
            self.listGene[ind].setRef(Qt.Checked)
        if hasattr(self, 'refSample'):
            soi = self.refSample.strip('"')
            ind = self.adresseEch[str(soi)]
            self.echRef = self.listEch[ind]
            self.listEch[ind].setRef(Qt.Checked)

    def findTriplicat(self, ectMax):
# Elimination du gene avec la chaine vide puis calcul du ctref
        for g in self.listGene[1:]:
            if self.dicoGene.has_key(g.name):
                g.calcCtRef(self.dicoGene[g.name])
                for well in self.dicoGene[g.name]:
                    well.setGene(g)
        self.dicoTrip = RaggedArray2D()
        for key in self.dicoGene.keys():
            dicoEch = RaggedArray2D()
            for well in self.dicoGene[str(key)]:
                if str(well.type) == 'unknown' and well.enabled == True:
                    if dicoEch.has_key(str(well.ech.name)):
                        dicoEch[str(well.ech.name)].append(well)
                    else:
                        dicoEch[str(well.ech.name)] = [well]
# Suppression de la chaine vide
            if dicoEch.has_key(""):
                dicoEch.pop("")
            for ech in dicoEch.keys():
                trip = Replicate(dicoEch[ech], ectMax=ectMax)
                if not hasattr(trip, "ctdev"):
                    # if ctdev undefined raise an exception
                    raise ValueError
                trip.calcDCt()
                dicoEch[ech] = trip
                self.dicoTrip[key] = dicoEch
            
    def findStd(self, ectMax):
        """
        This method allows to build a dictionnary for standard
        wells.
        """
        self.dicoStd = RaggedArray2D()
        for key in self.dicoGene.keys():
            dicoAmount = RaggedArray2D()
            for well in self.dicoGene[str(key)]:
                if str(well.type) == 'standard' and well.enabled == True:
                    if dicoAmount.has_key(str(well.amount)):
                        dicoAmount[str(well.amount)].append(well)
                    else:
                        dicoAmount[str(well.amount)] = [well]
            if dicoAmount.has_key(""):
                dicoAmount.pop("")
            for amount in dicoAmount.keys():
                trip = Replicate(dicoAmount[amount], type="standard",
                                 ectMax=ectMax)
                if not hasattr(trip, "ctdev"):
                    # if ctdev undefined raise an exception
                    raise ValueError
                dicoAmount[amount] = trip
                self.dicoStd[str(key)] = dicoAmount

    def calcNRQ(self):
        for g in self.dicoTrip.keys():
            for ech in self.dicoTrip[g].keys():
# Calcul de NRQ et rajout comme argument a chaque triplicat
                NRQ = self.dicoTrip[g][ech].RQ/ \
                    self.dicoTrip[g][self.echRef.name].RQ* \
                    self.dicoTrip[self.geneRef.name][self.echRef.name].RQ/ \
                    self.dicoTrip[self.geneRef.name][ech].RQ
                self.dicoTrip[g][ech].setNRQ(NRQ)
                self.dicoTrip[g][ech].calcRQerror()
                for well in self.dicoTrip[g][ech].listePuits:
                    well.setNRQ(NRQ)
# Calcul de NRQerror et rajout comme argument a chaque triplicat
        for g in self.dicoTrip.keys():
            for ech in self.dicoTrip[g].keys():
                NRQerror = self.dicoTrip[g][ech].NRQ  \
                        * sqrt((self.dicoTrip[self.geneRef.name][ech].RQerror \
                        / self.dicoTrip[self.geneRef.name][ech].RQ)**2 \
                        + (self.dicoTrip[g][ech].RQerror \
                        / self.dicoTrip[g][ech].RQ)**2 
# Rajout de 2 termes supplementaires avec notre definition de NRQ
                        + (self.dicoTrip[g][self.echRef.name].RQerror \
                        / self.dicoTrip[g][self.echRef.name].RQ)**2
                        + (self.dicoTrip[self.geneRef.name][self.echRef.name].RQerror \
                        / self.dicoTrip[self.geneRef.name][self.echRef.name].RQ)**2)
                self.dicoTrip[g][ech].setNRQerror(NRQerror)
                for well in self.dicoTrip[g][ech].listePuits:
                    well.setNRQerror(NRQerror)

    def calcStd(self, confidence):
        for geneName in self.dicoStd.keys():
            x = array([])
            y = array([])
            for trip in self.dicoStd[geneName].values():
                x = append(x, trip.amList)
                y = append(y, trip.ctList)
            x = log10(x)
            slope, orig = polyfit(x, y, 1)
            yest = polyval([slope, orig], x)
            seps = sqrt(sum((yest-y)**2)/(len(y)-2)) # Formule 2
            sx = sqrt(sum((x-x.mean())**2)/(len(x))) # Formule 3
            stderr = seps / (sx*sqrt(len(x))) # Formule 4 corrigee
            talpha = t.ppf(1.-(1.-confidence)/2., len(x)-2) # Student
            slopeerr = talpha * stderr
            eff = (10**(-1./slope)-1)*100 # Formule 5 adaptee
            # Erreur(Eff) = (Eff+100) * slopeerr / slope**2
            stdeff = (eff+100)*log(10)*slopeerr/slope**2 # Formule 6 adaptee
            # Coefficient de Pearsson de correlation
            R2 = 1 - sum((y-yest)**2)/sum((y-mean(y))**2)
            print eff, stdeff, R2
            # Mise a jour de l'efficacite des puits
            for well in self.dicoGene[geneName]:
                well.gene.setEff(eff)
                well.gene.setPm(stdeff)
            # il faut aussi mettre a jour les genes de listGene
            # qui servent a remplir les comboBox
            ind = self.adresseGene[geneName]
            self.listGene[ind].setEff(eff)
            self.listGene[ind].setPm(stdeff)


class Replicate(QDialog):

    def __init__(self, listePuits, type='unknown', parent=None, ectMax=0.3):
        self.parent = parent
        self.ectMax = ectMax
        QDialog.__init__(self, parent)
        self.type = type
        self.listePuits = listePuits
        if len(self.listePuits) != 0:
            self.gene = self.listePuits[0].gene
        else:
            self.gene = Gene('')
        self.ctList =array([])
        for well in self.listePuits:
            self.ctList = append(self.ctList, well.ct)
        if self.type == 'unknown':
            self.ech = self.listePuits[0].ech
        elif self.type == 'standard':
            self.amList = array([])
            for well in self.listePuits:
                self.amList = append(self.amList, well.amount)
        self.calcMeanDev()

    def __str__(self):
        st = '{%s:[' % self.type
        for well in self.listePuits:
            st = st + well.name + ','
        st += ']'
        if self.type == 'unknown':
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

            QMessageBox.warning(self, "Problem in calculation !",
                "A problem occured in the calculations. It seems to come from the \
                 well %s. Check whether ct are correctly defined." \
                 % (brokenWells))
            return

        if len(self.ctList) > 1:
            self.ctdev = self.ctList.std() * sqrt(len(self.ctList)/ \
                         (len(self.ctList)-1.))
        else:
            self.ctdev = 0.

        for well in self.listePuits:
            well.setCtmean(self.ctmean)
            well.setCtdev(self.ctdev)

        if self.ctdev >= self.ectMax:
            if self.type == 'unknown':
                QMessageBox.warning(self, "Warning Replicates",
                    "Warning: E(ct) of replicate (%s, %s) greater than %.2f" \
                    % (self.gene, self.ech, self.ectMax))
            elif self.type == 'standard':
                QMessageBox.warning(self, "Warning Replicates",
                    "Warning: E(ct) of replicate (%s, %s) greater than %.2f" \
                    % (self.gene, self.amList[0], self.ectMax))

    def calcDCt(self):
        self.dct = self.gene.ctref - self.ctmean # Formule 10
        self.RQ = (1.+self.gene.eff/100.)**(self.dct) # Formule 11

    def calcRQerror(self):
        # Formule 12
        err = sqrt( self.RQ**2 * ((self.dct*(self.gene.pm/100.) \
                /(1.+self.gene.eff/100.))**2 \
                + (log(1.+self.gene.eff/100.)*self.ctdev)**2 \
                ))
        self.RQerror = err


if __name__ == '__main__':
    pl = Plaque('../samples/raw_std.txt')
    print str(pl.A1.NRQ)
