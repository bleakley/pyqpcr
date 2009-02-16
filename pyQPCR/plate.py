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
from utils import *
from PyQt4.QtCore import Qt
from numpy import mean, std, sqrt, log, asarray
from pyQPCR.utils.odict import OrderedDict
from pyQPCR.utils.ragged import RaggedArray2D

__author__ = "$Author$:"
__date__ = "$Date$:"
__version__ = "$Rev$:"

class Plaque:
    
    def __init__(self, filename):
        self.filename = filename
#
        self.listePuits = []
        self.listGene = [Gene('')]
        self.listEch = [Ech('')]
        self.adresseGene = OrderedDict()
        self.adresseEch = OrderedDict()
        self.adresseEch[''] = 0
        self.adresseGene[''] = 0
# 
        self.determineFileType(self.filename)
        self.read()
        self.setDicoGene()
        self.setDicoEch()

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
        if self.fileType == "txt":
            splitter = re.compile(r'("[\w\s.,\-\(\)\[\]\+]*"|\d+[.,]?\d*)')
            iterator = file.readlines()
        if self.fileType == "csv":
            iterator = csv.reader(file, delimiter=";")
        for ind, line in enumerate(iterator):
            if self.fileType == "txt": line = splitter.findall(line)
            if ind == 0:
                self.header = []
                for field in line:
                    self.header.append(field.strip('"'))
                ncol = len(self.header)
            if len(line) == ncol and ind != 0:
                champs=[]
                for field in line:
                    try:
                        dat = float(field.replace(',', '.'))
                    except ValueError:
                        dat = field.strip('"')
                    champs.append(dat)
                name, ech, ct, ctmean, ctdev, amount, gene = champs[0:7]
                x = Puits(name, ech, ct, ctmean, ctdev, amount, gene)
# Dans le cas ou on a plus de 7 colonnes:
# le type est donne dans la colonne 7
# NRQ en 8 et NRQerror en 9
                if len(champs) > 7:
                    x.setType(champs[7])
                    x.setNRQ(champs[8])
                    x.setNRQerror(champs[9])
                setattr(self, x.name, x)
                self.listePuits.append(x)
                if not self.adresseEch.has_key(ech):
                    self.listEch.append(Ech(ech))
                    self.adresseEch[ech] = len(self.listEch)-1
# Pour les genes on fait un dico + une liste de genes
                if not self.adresseGene.has_key(gene):
                    self.listGene.append(Gene(gene))
                    self.adresseGene[gene] = len(self.listGene)-1 
            if len(line) == 2:
                name =  motif.findall(line[0].strip('"'))[0].replace(' ', '')
                value = line[1]
                try:
                    value = float(value)
                except ValueError:
                    pass
                setattr(self, name, value)
# Permet eventuellement de connaitre les genes/ech de ref 
        self.getRefsFromFile()
        file.close()

    def write(self, filename):
        f = open(filename, 'w')
        self.determineFileType(filename)
# Ajout eventuel dans le header
        if len(self.header) < 8:
            self.header.append('Type')
            self.header.append('NRQ')
            self.header.append('NRQerror')
# -----------TXT------------- 
        if self.fileType == "txt":
# Ecriture du header
            for field in self.header:
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
            f.write('"GOI"\t"%s"\n' % self.geneRef.name)
            f.write('"SOI"\t"%s"\n' % self.echRef.name)
# -----------CSV------------- 
        elif self.fileType == "csv":
            writer = csv.writer(f, quoting=csv.QUOTE_NONNUMERIC, delimiter=";")
# Ecriture du header
            writer.writerow(self.header)
# Ecriture des puits
            for well in self.listePuits:
                writer.writerow([well.name, well.ech.name, well.ct, well.ctmean,
                    well.ctdev, well.amount, well.gene.name, well.type, well.NRQ,
                    well.NRQerror])
# Ecriture du `Analysis Parameters`
            writer.writerow("")
            writer.writerow(['Analysis Parameters'])
            writer.writerow(['GOI', self.geneRef.name])
            writer.writerow(['SOI', self.echRef.name])
        f.close()

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
            if not self.adresseGene.has_key(str(gene.name)):
                self.adresseGene[str(gene.name)] = ind
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

    def setRefs(self):
        """
        Determine the reference target and sample
        """
        for g in self.listGene:
            if g.isRef == 2:
                self.geneRef = g
        for e in self.listEch:
            if e.isRef == 2:
                self.echRef = e

    def getRefsFromFile(self):
        """
        This methods allows to determine both the Gene of Interest and
        the reference-sample thanks to the input data file
        """
        if hasattr(self, 'GOI'):
            goi = self.GOI.strip('"')
            ind = self.adresseGene[str(goi)]
            self.geneRef = self.listGene[ind]
            self.listGene[ind].setRef(Qt.Checked)
        if hasattr(self, 'SOI'):
            soi = self.SOI.strip('"')
            ind = self.adresseEch[str(soi)]
            self.echRef = self.listEch[ind]
            self.listEch[ind].setRef(Qt.Checked)

    def findTriplicat(self):
# Elimination du gene avec la chaine vide
        for g in self.listGene[1:]:
            if self.dicoGene.has_key(g.name):
                g.calcCtRef(self.dicoGene[g.name])
                for well in self.dicoGene[g.name]:
                    well.setGene(g)
        self.dicoTrip = OrderedDict()
        self.dicoTrip = RaggedArray2D()
        for key in self.dicoGene.keys():
            #dicoEch = OrderedDict()
            dicoEch = RaggedArray2D()
            dicoAmount = OrderedDict()
            for well in self.dicoGene[str(key)]:
                if str(well.type) == 'unknown' and well.enabled == True:
                    if dicoEch.has_key(str(well.ech.name)):
                        dicoEch[str(well.ech.name)].append(well)
                    else:
                        dicoEch[str(well.ech.name)] = [well]
                if str(well.type) == 'standard' and well.enabled == False:
                    if dicoAmount.has_key(well.amount):
                        dicoAmount[well.amount].append(well)
                    else:
                        dicoAmount[well.amount] = [well]
# Suppression de la chaine vide
            try:
                dicoEch.pop("")
            except KeyError:
                pass
            for ech in dicoEch.keys():
                trip = Triplicat(dicoEch[ech], str(dicoEch[ech][0].type))
                trip.calcDCt()
                dicoEch[ech] = trip
                self.dicoTrip[key] = dicoEch
            
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

class Triplicat:

    def __init__(self, listePuits, type):
        self.listePuits = listePuits
        self.type = type
        self.gene = self.listePuits[0].gene
        self.ech = self.listePuits[0].ech
        self.amount = self.listePuits[0].amount
        #
        self.ctList =[]
        for well in self.listePuits:
            self.ctList.append(well.ct)
        self.calcMeanDev()

    def __str__(self):
        st = '{%s:[' % self.type
        for well in self.listePuits:
            st = st + well.name + ','
        st += ']'
        if self.type == 'unknown':
            st += ' %s, %s}' % (self.gene, self.ech)
        elif self.type == 'standard':
            st += ' %s, %s}' % (self.gene, str(self.amount))
        return st

    def __repr__(self):
        st = '{%s:[' % self.type
        for well in self.listePuits:
            st = st + well.name + ','
        st += ']'
        if self.type == 'unknown':
            st += ' %s, %s}' % (self.gene.name, self.ech)
        elif self.type == 'standard':
            st += ' %s, %s}' % (self.gene.name, str(self.amount))
        return st

    def setNRQ(self, NRQ):
        self.NRQ = NRQ

    def setNRQerror(self, NRQerr):
        self.NRQerror = NRQerr

    def calcMeanDev(self):
        """
        Description:
        -----------
           compute the mean ct of a replicate

        """
        self.ctmean = mean(self.ctList)
        #self.ctdev = std(self.ctList, ddof=1)
        self.ctdev = asarray(self.ctList).std()* sqrt(len(self.ctList)/ \
                     (len(self.ctList)-1.))

        for well in self.listePuits:
            well.setCtmean(self.ctmean)
            well.setCtdev(self.ctdev)

    def calcDCt(self):
        self.dct = self.gene.ctref - self.ctmean
        self.RQ = (1.+self.gene.eff/100.)**(self.dct)

    def calcRQerror(self):
        err = sqrt( self.RQ**2 * ((self.dct*(self.gene.pm/100.) \
                /(1.+self.gene.eff/100.))**2 \
                + (log(1.+self.gene.eff/100.)*self.ctdev)**2 \
                ))
        self.RQerror = err


if __name__ == '__main__':
    pl = Plaque('sortiesrealplex/ref-machine.txt')
    print str(pl.A1.NRQ)
    print pl.GOI
    pl.write('toto.csv')

