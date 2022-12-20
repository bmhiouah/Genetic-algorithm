from random import uniform, sample, shuffle, randint
from Individu import Individu
from copy import deepcopy

# classe composée d'une liste de n individus

class ListIndividu:

    

    def __init__(self,nb_indiv,seq):
        self.seq = seq
        self.nb_indiv = nb_indiv
        self.list_indiv = [None for i in range(nb_indiv)]
    


    def initial_list(self):

        # on génère des individus sur la plage [valeur - ecart, valeur + ecart]
        for i in range(self.nb_indiv):
            indiv = Individu()
            table = indiv.Table
            for dinucleotide in indiv.dinu_to_modif:
                for j in range(3):
                    new_value = table[dinucleotide][j] + uniform(-1,1)*table[dinucleotide][j+3]
                    indiv.modif(dinucleotide,j,new_value)
            self.change_indiv(i,indiv)


    def change_indiv(self,i,indiv): #permet de changer un élément de la liste par un autre individu
        self.list_indiv[i] = deepcopy(indiv)
    

    def select(self): #c'est la fonction selection par tournoi, elle prend une liste d'individu et retourne une autre liste d'individu
        shuffle(self.list_indiv) #on commence par mélanger pour créer de la diversité 
        for i in range(0,self.nb_indiv-1,2): #on les fais combattre 2 à 2
            comb1, comb2 = self.get_indiv(i), self.get_indiv(i+1)
            if comb1.fitness(self.seq) < comb2.fitness(self.seq): #et enfin on remplace le perdant par le gagnant on a ainsi une liste de doublon de gagnant
                self.change_indiv(i+1,comb1)
            else:
                self.change_indiv(i,comb2)
        


    def fusion(self,i,j,k): #fusionne les individus d'indice i et j et mets les enfants à leur place
        indiv1, indiv2 = self.get_indiv(i), self.get_indiv(j)
        choosen1 = sample(indiv1.dinu_to_modif, k) # k indique le nombre d'élements qui seront interchangé.

        for dinu in indiv1.dinu_to_modif:
            
            if dinu in choosen1:
                for n in range(3):
                    x = indiv1.Table[dinu][n]
                    y = indiv2.Table[dinu][n]
                    indiv2.modif(dinu,n,x)
                    indiv1.modif(dinu,n,y)
            

    def fusion_total(self): #c'est la fusion appliqué à toute la liste, elle prend une liste et change la liste
        for i in range(0,self.nb_indiv,4): #la fonction est constitué de manière à laisser le premier doublon intact et se servir du deuxième doublon pour faire une fusion
            k = randint(1,8) 
            self.fusion(i,i+2,k)


    def get_indiv(self,k): # renvoie l'individu au rang k de self
        return self.list_indiv[k]


    def mutate(self, p, coef = 1): #applique une mutation à toute la liste sauf le meilleur
        best_indiv=self.get_best()
        for i in range(self.nb_indiv):
            if self.get_indiv(i) != best_indiv:
                self.get_indiv(i).mutation(p,coef)
    

    def get_best(self): #permet de trouver l'individu de la liste avec le meilleur résultat
        best = self.get_indiv(0)
        for i in range(self.nb_indiv):
            if self.get_indiv(i).fitness(self.seq) < best.fitness(self.seq):
                best = self.get_indiv(i)
        return best

        
