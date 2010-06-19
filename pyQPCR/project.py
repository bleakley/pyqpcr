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
#
from pyQPCR.saxProjecthandler import *
from pyQPCR.utils.odict import OrderedDict
from pyQPCR.utils.ragged import RaggedArray2D
from pyQPCR.plate import StdObject
from PyQt4.QtCore import *
from PyQt4.QtXml import QXmlSimpleReader, QXmlInputSource
from numpy import mean, std, sqrt, log, log10, polyval, polyfit, sum, \
array, append

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Rev$"

class Project:

    def __init__(self, fname=None):
        self.dicoPlates = OrderedDict()
#
        self.hashGene = OrderedDict()
        self.hashGene[''] = Gene('')
#
        self.hashEch = OrderedDict()
        self.hashEch[''] = Ech('')
#
        self.hashAmount = OrderedDict()
        self.hashAmount[''] = ''
        self.unsaved = False
        self.filename = fname
        if self.filename is not None:
            self.openProject(self.filename)

    def openProject(self, fname):
        error = None
        fh = None
        try:
            handler = SaxProjectHandler(self)
            parser = QXmlSimpleReader()
            parser.setContentHandler(handler)
            parser.setErrorHandler(handler)
            fh = QFile(fname)
            input = QXmlInputSource(fh)
            if not parser.parse(input):
                raise ValueError, handler.error
# For each plate, adress the Gene location
            self.initLocGene()
            self.initLocEch()
            self.initLocAm()
# A deplacer :
            for pl in self.dicoPlates.values():
                pl.setDicoGene()
                pl.setDicoEch()
                if len(pl.geneRef) != 0:
                    for geneName in pl.geneRef:
                        if self.hashGene.has_key(geneName):
                            self.hashGene[geneName].setRef(Qt.Checked)
                        else:
                            pl.geneRef.remove(geneName)
                if pl.echRef != '' and self.hashEch.has_key(pl.echRef):
                    self.hashEch[pl.echRef].setRef(Qt.Checked)
            self.setDicoAm()
        except (IOError, OSError, ValueError), e:
            error = "Failed to import: %s" % e
        finally:
            if fh is not None:
                fh.close()
            if error is not None:
                return False, error

    def exportXml(self, fname):
        error = None
        CODEC = "UTF-8"
        try:
            fh = QFile(fname)
            if not fh.open(QIODevice.WriteOnly):
                raise IOError, unicode(fh.errorString())
            stream = QTextStream(fh)
            stream.setCodec(CODEC)
            stream << ("<?xml version='1.0' encoding='%s'?>\n"
                       "<!DOCTYPE QPCR>\n"
                       "<QPCR VERSION='1.0'>\n" % CODEC)
            for key in self.dicoPlates.keys():
                stream << ("<PLATE NAME='%s'>\n" % key)
                for geneName in self.dicoPlates[key].geneRef:
                    if geneName != '':
                        stream << ("<REFTARGET NAME='%s'></REFTARGET>\n") % \
                                   geneName
                stream << ("<REFSAMPLE NAME='%s'></REFSAMPLE>\n") % \
                        self.dicoPlates[key].echRef
                for well in self.dicoPlates[key].listePuits:
                    stream << well.writeWellXml()
                stream << "</PLATE>\n"
            stream << "</QPCR>\n"
        except (IOError, OSError), e:
            error = "Failed to export: %s" % e
        finally:
            if fh is not None:
                fh.close()
            if error is not None:
                return False, error
            
    def addPlate(self, plate, key=None):
        """
        This method allows to add a plate to an existing Project.

        @param plate: the plate we want to add
        @type plate: pyQPCR.Plaque
        @param key: the key associated with the plate
        @type key: PyQt4.QtCore.QString
        """
        if not key:
            key = QFileInfo(plate.filename).fileName()
        self.dicoPlates[key] = plate
# For each plate, adress the Gene location
        self.initLocGene(plate)
        self.initLocEch(plate)
        self.initLocAm(plate)
# A deplacer :
        plate.setDicoGene()
        plate.setDicoEch()
        self.setDicoAm()

    def removePlate(self, plateName):
        """
        This method allows to delete a specific plate of an existing
        project.

        @param plateName: the name of the plate we want to delete
        @type plateName: string
        """
        oldGenes = self.dicoPlates[plateName].dicoGene.keys()
        oldEchs = self.dicoPlates[plateName].dicoEch.keys()
        oldGenesRef = self.dicoPlates[plateName].geneRef
        oldEchRef = self.dicoPlates[plateName].echRef
        self.dicoPlates.__delitem__(plateName)

        listGeneRef = []
        for oldg in oldGenes:
            delete = True
            for pl in self.dicoPlates.values():
                if pl.dicoGene.has_key(oldg) or oldg == '':
                    delete = False
                    if oldGenesRef.__contains__(oldg):
                        listGeneRef.append(oldg)
            if delete:
                self.hashGene.__delitem__(oldg)

        if len(listGeneRef) > 0:
            uncheck = True
            for geneName in listGeneRef:
                for pl in self.dicoPlates.values():
                    if pl.geneRef.__contains__(geneName):
                        uncheck = False
                if uncheck:
                    self.hashGene[geneName].setRef(Qt.Unchecked)

        listEchRef = []
        for oldech in oldEchs:
            delete = True
            for pl in self.dicoPlates.values():
                if pl.dicoEch.has_key(oldech) or oldech == '':
                    delete = False
                    if oldEchRef == oldech:
                        listEchRef.append(oldech)
            if delete:
                self.hashEch.__delitem__(oldech)

        if len(listEchRef) > 0:
            uncheck = True
            for echName in listEchRef:
                for pl in self.dicoPlates.values():
                    if pl.echRef == echName:
                        uncheck = False
                if uncheck:
                    self.hashEch[echName].setRef(Qt.Unchecked)

                
# A deplacer :
        self.setDicoAm()

    def initLocGene(self, plate=None):
        if plate is None:
            for pl in self.dicoPlates:
                for well in self.dicoPlates[pl].listePuits:
                    nomgene = well.gene.name
                    if not self.hashGene.has_key(nomgene):
                        self.hashGene[nomgene] = well.gene
        else:
            for well in plate.listePuits:
                nomgene = well.gene.name
                if not self.hashGene.has_key(nomgene):
                    self.hashGene[nomgene] = well.gene
                    if hasattr(plate, 'geneRef'):
                        for geneName in plate.geneRef:
                            if geneName == nomgene:
                                self.hashGene[nomgene].setRef(Qt.Checked)

    def initLocEch(self, plate=None):
        if plate is None:
            for pl in self.dicoPlates:
                for well in self.dicoPlates[pl].listePuits:
                    nomech = well.ech.name
                    if not self.hashEch.has_key(nomech):
                        self.hashEch[nomech] = well.ech
        else:
            for well in plate.listePuits:
                nomech = well.ech.name
                if not self.hashEch.has_key(nomech):
                    self.hashEch[nomech] = well.ech
                    if hasattr(plate, 'echRef'):
                        if plate.echRef == nomech:
                            self.hashEch[nomech].setRef(Qt.Checked)

    def initLocAm(self, plate=None):
        if plate is None:
            for pl in self.dicoPlates:
                for well in self.dicoPlates[pl].listePuits:
                    try:
                        key = QString("%.2f" % well.amount)
                    except TypeError:
                        key = QString(well.amount)
                    if not self.hashAmount.has_key(key):
                        self.hashAmount[key] = well.amount
        else:
            for well in plate.listePuits:
                try:
                    key = QString("%.2f" % well.amount)
                except TypeError:
                    key = QString(well.amount)
                if not self.hashAmount.has_key(key):
                    self.hashAmount[key] = well.amount

    def setDicoAm(self):
        self.dicoAmount = OrderedDict()
        for pl in self.dicoPlates:
            for well in self.dicoPlates[pl].listePuits:
                try:
                    key = QString("%.2f" % well.amount)
                except TypeError:
                    key = QString(well.amount)
                if self.dicoAmount.has_key(key):
                    self.dicoAmount[key].append(well)
                else:
                    self.dicoAmount[key] = [well]

    def findTrip(self, ectMax, confidence, errtype):
        self.dicoTriplicat = OrderedDict()
        largeCtTrip = []
        for plate in self.dicoPlates.keys():
            pl = self.dicoPlates[plate]
            for g in self.hashGene.values()[1:]:
                if pl.dicoGene.has_key(g.name):
                    g.calcCtRef(pl.dicoGene[g.name])
                    for well in pl.dicoGene[g.name]:
                        well.setGene(g)
            dicoTrip = RaggedArray2D()
            for key in pl.dicoGene.keys():
                dicoEch = RaggedArray2D()
                for well in pl.dicoGene[key]:
                    if well.type == QString('unknown') and well.enabled == True:
                        if dicoEch.has_key(well.ech.name):
                            dicoEch[well.ech.name].append(well)
                        else:
                            dicoEch[well.ech.name] = [well]
# Suppression de la chaine vide
                if dicoEch.has_key(""):
                    dicoEch.pop("")
                for ech in dicoEch.keys():
                    trip = Replicate(dicoEch[ech], confidence=confidence, 
                                     errtype=errtype)
                    if not hasattr(trip, "ctdev"):
                        # if ctdev undefined raise an exception
                        raise ValueError
                    if trip.ctdev >= ectMax:
                        largeCtTrip.append(trip)
                    trip.calcDCt()
                    trip.calcRQerror()
                    dicoEch[ech] = trip
                    dicoTrip[key] = dicoEch
            self.dicoTriplicat[plate] = dicoTrip
        if len(largeCtTrip) != 0:
            raise ReplicateError(largeCtTrip)

    def calcCF(self):
        self.CF = OrderedDict()
        self.CFerror = OrderedDict()
        for pl in self.dicoTriplicat.keys():
            for g in self.dicoTriplicat[pl].keys():
                for pl1 in self.dicoTriplicat.keys():
                    if pl1 != pl and self.dicoTriplicat[pl1].has_key(g):
                        for ech in self.dicoTriplicat[pl][g].keys():
                            if self.dicoTriplicat[pl1][g].has_key(ech):
                                NRQ = self.dicoTriplicat[pl][g][ech].NRQ
                                NRQerror = self.dicoTriplicat[pl][g][ech].NRQerror
                                if self.CF.has_key(pl):
                                    self.CF[pl] = append(self.CF[pl], NRQ)
                                    self.CFerror[pl] = append(self.CFerror[pl], 
                                            (NRQerror/NRQ)**2)
                                else:
                                    self.CF[pl] = array([NRQ])
                                    self.CFerror[pl] = array([(NRQerror/NRQ)**2])
            self.CF[pl] = pow(self.CF[pl].prod(),1./len(self.CF[pl]))
            self.CFerror[pl] = self.CF[pl] * sqrt(self.CFerror[pl].sum())

    def calcNRQ(self):
        broken = []
        NF = {} ; NFerror = {}
        for pl in self.dicoTriplicat.keys():
            plate = self.dicoPlates[pl]
# Boucle de calcul des NF : plusieurs genes de references
            for ech in self.hashEch.keys():
                tabRQ = array([]) ; tabRQerror = array([])
                for g in plate.geneRef:
                    try:
                        if self.dicoTriplicat[pl][g].has_key(ech):
                            tmp = self.dicoTriplicat[pl][g][ech].RQ
                            tmp2 = (self.dicoTriplicat[pl][g][ech].RQerror / \
                                   self.dicoTriplicat[pl][g][ech].RQ)**2
                            tabRQ = append(tabRQ, tmp)
                            tabRQerror = append(tabRQerror, tmp2)
                    except KeyError,e:
                        if ech != '':
                            broken.append((g, ech))
                        continue
                if len(tabRQ) != 0:
                    NF[ech] = (tabRQ.prod())**(1./len(tabRQ))
                    NFerror[ech] = NF[ech]*sqrt(tabRQerror.sum())/len(tabRQ)
            if len(broken) != 0:
                raise NRQError(broken)
            for g in self.dicoTriplicat[pl].keys():
                for ech in self.dicoTriplicat[pl][g].keys():
# Calcul de NRQ et rajout comme argument a chaque triplicat
                    try:
                        NRQ = self.dicoTriplicat[pl][g][ech].RQ/ \
                              self.dicoTriplicat[pl][g][plate.echRef].RQ* \
                              NF[plate.echRef]/ NF[ech]
                        self.dicoTriplicat[pl][g][ech].setNRQ(NRQ)
                        for well in self.dicoTriplicat[pl][g][ech].listePuits:
                            #print pl, well.name, NRQ
                            well.setNRQ(NRQ)
                    except KeyError, ke:
                        broken.append((g,ech))
                        continue
# Calcul de NRQerror et rajout comme argument a chaque triplicat
            for g in self.dicoTriplicat[pl].keys():
                for ech in self.dicoTriplicat[pl][g].keys():
                    try:
                        NRQerror = self.dicoTriplicat[pl][g][ech].NRQ  \
                             * sqrt((NFerror[ech] / NF[ech])**2 \
                             + (self.dicoTriplicat[pl][g][ech].RQerror \
                             / self.dicoTriplicat[pl][g][ech].RQ)**2  \
                             + (self.dicoTriplicat[pl][g][plate.echRef].RQerror \
                             / self.dicoTriplicat[pl][g][plate.echRef].RQ)**2 \
                             + (NFerror[plate.echRef] / NF[plate.echRef])**2)
                        self.dicoTriplicat[pl][g][ech].setNRQerror(NRQerror)
                        for well in self.dicoTriplicat[pl][g][ech].listePuits:
                            well.setNRQerror(NRQerror)
                    except (KeyError, AttributeError):
                        continue
        if len(broken) != 0:
            raise NRQError(broken)

    def findStd(self, ectMax, confidence, errtype):
        """
        This method allows to build a dictionnary for standard
        wells.
        """
        self.dicoStd = RaggedArray2D()
        largeCtTrip = []
        for pl in self.dicoPlates.values():
            for key in pl.dicoGene.keys():
                dicoAmount = RaggedArray2D()
                for well in pl.dicoGene[key]:
                    if well.type == QString('standard') and well.enabled == True:
                        if dicoAmount.has_key(str(well.amount)):
                            dicoAmount[str(well.amount)].append(well)
                        else:
                            dicoAmount[str(well.amount)] = [well]
                if dicoAmount.has_key(""):
                    dicoAmount.pop("")
                for amount in dicoAmount.keys():
                    trip = Replicate(dicoAmount[amount], type=QString('standard'),
                                     confidence=confidence, errtype=errtype)
                    if not hasattr(trip, "ctdev"):
                        # if ctdev undefined raise an exception
                        raise ValueError
                    if trip.ctdev >= ectMax:
                        largeCtTrip.append(trip)
                    dicoAmount[amount] = trip
                    self.dicoStd[key] = dicoAmount
        if len(largeCtTrip) != 0:
            raise ReplicateError(largeCtTrip)

    def calcStd(self, confidence, errtype):
        self.dicoPlotStd = OrderedDict()
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
            if errtype == "student":
                talpha = t.ppf(1.-(1.-confidence)/2., len(x)-2) # Student
            elif errtype == "normal":
                talpha = norm.ppf(1.-(1.-confidence)/2.) # Gaussian
            slopeerr = talpha * stderr
            eff = (10**(-1./slope)-1)*100 # Formule 5 adaptee
            # Erreur(Eff) = (Eff+100) * slopeerr / slope**2
            stdeff = (eff+100)*log(10)*slopeerr/slope**2 # Formule 6 adaptee
            # Coefficient de Pearsson de correlation
            R2 = 1 - sum((y-yest)**2)/sum((y-mean(y))**2)
            self.dicoPlotStd[geneName] = StdObject(x, y, yest, slope, orig, R2, eff, stdeff)
            # output for debugging stuff:
            # print eff, stdeff, R2
            # Mise a jour de l'efficacite des puits
            for pl in self.dicoPlates.values():
                if pl.dicoGene.has_key(geneName):
                    for well in pl.dicoGene[geneName]:
                        well.gene.setEff(eff)
                        well.gene.setPm(stdeff)
            # il faut aussi mettre a jour les genes de listGene
            # qui servent a remplir les comboBox
            self.hashGene[geneName].setEff(eff)
            self.hashGene[geneName].setPm(stdeff)

    def findBars(self, width, spacing, sens='geneEch', plates=None):
        """
        Fonction qui determine le nb de barres dans chaque categorie
        (gene, ech) en vue du trace. Le dictionnaire barWidth est d'abord
        rempli du nombre de barre:
           ex. barWidth[ctrl] = 8
        Puis a partir de la, on trouve les abscisses correspondantes.
        """
        leftMargin = 0.1
        if plates is None:
            plates = self.dicoPlates.keys() 
        self.barWidth = OrderedDict()
        self.barXticks = OrderedDict()
        for g in self.hashGene.keys()[1:]:
            for ech in self.hashEch.keys()[1:]:
                for pl in plates:
                    if self.dicoTriplicat[pl].has_key(g) and \
                            self.hashGene[g].enabled == Qt.Checked and \
                            self.dicoTriplicat[pl][g].has_key(ech) and \
                            self.hashEch[ech].enabled == Qt.Checked:
                        if sens == 'geneEch':
                            if self.barWidth.has_key(ech):
                                self.barWidth[ech] += 1
                            else:
                                self.barWidth[ech] = 1
                        elif sens == 'echGene':
                            if self.barWidth.has_key(g):
                                self.barWidth[g] += 1
                            else:
                                self.barWidth[g] = 1
        nbar = array(self.barWidth.values())
        largeur = nbar*width+spacing
        i = 0
        for e in self.barWidth.keys():
            self.barWidth[e] = largeur[:i].sum() + leftMargin
# matplotlib veut des str et non des QString
            self.barXticks[str(e)] = largeur[:i].sum() + \
                                     leftMargin +(nbar[i]-1.)*width/2.
            i += 1

    def findUnknown(self):
        """
        This methods computes a list of length nplate which contains booleans
        which indicate wheter a plate contains 'unknown' wells or not.
        """
        for pl in self.dicoPlates:
            liste = []
            for well in self.dicoPlates[pl].listePuits:
                if well.type == QString('unknown'):
                    liste.append(well)

            if len(liste) == 0:
                self.dicoPlates[pl].setUkn(False)
            else:
                self.dicoPlates[pl].setUkn(True)


class NRQError(Exception):

    def __init__(self, broken):
        self.broken = broken

    def __str__(self):
        st = "<ul>"
        for line in self.broken:
            st += "<li>(<b>%s, %s</b>)</li>" % (line[0], line[1])
        st += "</ul>"
        return st


if __name__ == "__main__":
    proj = Project('toto.xml')
    print proj.dicoPlates[QString('mh101109-1m.TXT')].echRef
    proj.findTrip(0.3, 0.9, 'student')
    print proj.dicoTriplicat[QString('mh101109-1m.TXT')]
    proj.calcNRQ()

