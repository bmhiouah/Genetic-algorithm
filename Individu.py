## Classe individu : RotTable sans les écarts-types

import mathutils
import math
import random

from Traj3D import Traj3D


class Individu:
    

    couplage={"AA":"TT","AC":"GT","AG":"CT","CA":"TG","CC":"GG","CT":"AG","GA":"TC","GG":"CC","GT":"AC","TC":"GA","TG":"CA","TT":"AA"}
    
    __ORIGINAL_ROT_TABLE = {\
        "AA": [35.62, 7.2, -154, 0.06, 0.6, 0],\
        "AC": [34.4, 1.1, 143, 1.3, 5, 0],\
        "AG": [27.7, 8.4, 2, 1.5, 3, 0],\
        "AT": [31.5, 2.6, 0, 1.1, 2, 0],\
        "CA": [34.5, 3.5, -64, 0.9, 34, 0],\
        "CC": [33.67, 2.1, -57, 0.07, 2.1, 0],\
        "CG": [29.8, 6.7, 0, 1.1, 1.5, 0],\
        "CT": [27.7, 8.4, -2, 1.5, 3, 0],\
        "GA": [36.9, 5.3, 120, 0.9, 6, 0],\
        "GC": [40, 5, 180, 1.2, 1.275, 0],\
        "GG": [33.67, 2.1, 57, 0.07, 2.1, 0],\
        "GT": [34.4, 1.1, -143, 1.3, 5, 0],\
        "TA": [36, 0.9, 0, 1.1, 2, 0],\
        "TC": [36.9, 5.3, -120, 0.9, 6, 0],\
        "TG": [34.5, 3.5, 64, 0.9, 34, 0],\
        "TT": [35.62, 7.2, 154, 0.06, 0.6, 0]\
        }

    dinu_to_modif = ["AA", "AC", "AG", "AT", "CA", "CC", "GA", "CG", "GC", "TA"] #liste des dinucleotides suffisantes à modifier pour modifier toute la rot table (avec les symetries)
    
    def __init__(self,Table = __ORIGINAL_ROT_TABLE):
        self.Table = {}
        for dinucleotide in Table:
            self.Table[dinucleotide] = Table[dinucleotide][:]
        

    def mutation(self,p,coef = 1): #p: probabilité de mutation #coef = à quel point on s'autorise à sécarter de la moyenne, il est superdicielle ici
        if p>random.random():
            dinucleotide=self.dinu_to_modif[random.randint(0,9)] #choix équiprobable du dinucléotide à modifier
            n=random.randint(0,1)#choix de l'angle à modifier
            ecart=self.Table[dinucleotide][n+3] #récupération de l'écart à ne pas dépasser
            alea=random.triangular(-1,1)
            variation=ecart*coef*alea
            valeur=self.Table[dinucleotide][n]+variation #nouvelle valeur aprés mutation
            self.modif(dinucleotide,n,valeur) #application du changement
    

    #la fonction suivante permet de rester cohérent lors des modifications et de conserver les symétries.
    def modif(self,dinucleotide,n,val):#modifie la valeur n du dinucleotide et la remplace par val
        mini=self.__ORIGINAL_ROT_TABLE[dinucleotide][n]-self.__ORIGINAL_ROT_TABLE[dinucleotide][n+3]
        maxi=self.__ORIGINAL_ROT_TABLE[dinucleotide][n]+self.__ORIGINAL_ROT_TABLE[dinucleotide][n+3]
        if val>=mini and val<=maxi:#val doit rester entre une valeur inimale et une valeur maximale
            self.Table[dinucleotide][n]=val
            if dinucleotide in self.couplage:#chaque dinucléotide est couplé à une autre et donc on doit aussi modifié cette autre valeur
                couple=self.couplage[dinucleotide]
                if n==2:
                    self.Table[couple][n]=-val  

                else:
                    self.Table[couple][n]=val                 


    def getTwist(self, dinucleotide):
        return self.Table[dinucleotide][0]

    def getWedge(self, dinucleotide):
        return self.Table[dinucleotide][1]

    def getDirection(self, dinucleotide):
        return self.Table[dinucleotide][2]

    
    def fitness(self,seq): #donne pour l'individu et la dna_seq = seq l'écart entre le premier et le dernier élément de la chaîne
        traj = Traj3D()
        traj.compute(seq, self)
        a = traj.getTraj()
        b = len(a)
       
        return ((a[0][0] - a[b-1][0])**2 + (a[0][1] - a[b-1][1])**2 + (a[0][2] - a[b-1][2])**2)**0.5



