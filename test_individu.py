from Individu import Individu
import random

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



def test_init():
    Tab=__ORIGINAL_ROT_TABLE
    ind=Individu(Tab)
    assert ind.Table==Tab
    



def test_mutation():
    Tab=__ORIGINAL_ROT_TABLE
    nbre_mutation=0
    for i in range(10000):
        ind=Individu(Tab)
        ind.mutation(0.1,coef=1)
        for dinucleotide in Tab:
            for n in range(2):
                if Tab[dinucleotide][n]!=ind.Table[dinucleotide][n]:
                    nbre_mutation+=1
                    mini=(__ORIGINAL_ROT_TABLE[dinucleotide][n])-(__ORIGINAL_ROT_TABLE[dinucleotide][n+3])
                    maxi=(__ORIGINAL_ROT_TABLE[dinucleotide][n])+(__ORIGINAL_ROT_TABLE[dinucleotide][n+3])
                    assert (ind.Table[dinucleotide][n]<=maxi and ind.Table[dinucleotide][n]>=mini)
    print(nbre_mutation/(10000*1.6))
    assert abs((nbre_mutation/(10000*1.6))-0.1)<0.05




def test_modif():
    Tab=__ORIGINAL_ROT_TABLE
    ind=Individu(Tab)
    for dinucleotide in ind.dinu_to_modif:
        for n in range(2):
            rand=random.uniform(-2,2)
            valeur=Tab[dinucleotide][n]+Tab[dinucleotide][n+3]*rand
            ind.modif(dinucleotide,n,valeur)
            if abs(rand)<=1:
                assert ind.Table[dinucleotide][n]==valeur
            else:
                assert ind.Table[dinucleotide][n]==Tab[dinucleotide][n]
