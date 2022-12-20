from ListIndividu import ListIndividu
from Individu import Individu
from copy import deepcopy
def test_init():
    init = ListIndividu(5,"ACTTCGGTAACACGTGAC")
    assert init.list_indiv == [None for i in range(5)]
    init = ListIndividu(10,"ACTTCGGTAACACGTGAC")
    assert init.list_indiv == [None for i in range(10)]

def test_initial_list():
    L = ListIndividu(5,"ACTTCGGTAACACGTGAC")
    L.initial_list()
    assert L.nb_indiv == 5
    assert len(L.list_indiv) == 5
    indiv = Individu()
    
    for dinucleotide in indiv.Table:
        assert dinucleotide in L.list_indiv[0].Table
    
        assert indiv.getDirection(dinucleotide) != L.list_indiv[0].getDirection(dinucleotide) or indiv.Table[dinucleotide][5] == 0
        assert indiv.getWedge(dinucleotide) != L.list_indiv[0].getWedge(dinucleotide) or indiv.Table[dinucleotide][4] == 0
        assert indiv.getTwist(dinucleotide) != L.list_indiv[0].getTwist(dinucleotide) or indiv.Table[dinucleotide][3] == 0

        assert indiv.Table[dinucleotide][3] == L.list_indiv[0].Table[dinucleotide][3]
        assert indiv.Table[dinucleotide][4] == L.list_indiv[0].Table[dinucleotide][4]
        assert indiv.Table[dinucleotide][5] == L.list_indiv[0].Table[dinucleotide][5]



def test_change_indiv():
    L = ListIndividu(12,"ACTTCGGTAACACGTGAC")
    L.change_indiv(2,"abc")
    assert L.list_indiv[2] == "abc"

def test_select():
    lis=ListIndividu(40,"AAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAG")
    lis.initial_list()
    best_elt=lis.list_indiv[0]
    dist=lis.list_indiv[0].fitness("AAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAG")
    for elt in lis.list_indiv:
        if elt.fitness("AAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAG")<dist:
            dist=elt.fitness("AAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAG")
            best_elt=elt
    lis.select()
    assert best_elt in lis.list_indiv
    print(lis.list_indiv)
    for i in range(20):
        print(i)
        assert lis.list_indiv[2*i].Table==lis.list_indiv[2*i+1].Table

def test_fusion():
    lis=ListIndividu(40,"AAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAG")
    lis.initial_list()
    ind=Individu()
    tab1=deepcopy(ind.Table)
    lis.change_indiv(0,ind)
    ind2=lis.list_indiv[1]
    tab2=deepcopy(ind2.Table)
    lis.fusion(0,1,0)
    assert lis.get_indiv(0).Table==ind.Table
    assert lis.get_indiv(1).Table==ind2.Table
    lis.fusion(0,1,10)
    for clé in tab2:
        for i in range(6):
            assert lis.get_indiv(0).Table[clé][i]==tab2[clé][i]
            assert lis.get_indiv(1).Table[clé][i]==tab1[clé][i]

def test_indiv():
    l_ind=ListIndividu(12,"AAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAG")
    ind=Individu()
    l_ind.change_indiv(0,ind)
    assert l_ind.get_indiv(0).Table==ind.Table


def test_get_best():
    seq="AAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAG"
    l_ind=ListIndividu(12,seq)
    l_ind.initial_list()
    best=l_ind.get_best()
    for ind in l_ind.list_indiv:
        assert best.fitness(seq)<=ind.fitness(seq)