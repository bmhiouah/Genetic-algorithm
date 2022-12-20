from ListIndividu import ListIndividu
from Individu import Individu

from Traj3D import Traj3D

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--filename", help="input filename of DNA sequence", default= "plasmid_8k.fasta")
parser.add_argument("--nb_gen", help="input number of generations")
parser.add_argument("--nb_indiv", help="input size of the population")
parser.parse_args()
args = parser.parse_args()

def main(nb_gen = 60, nb_indiv = 100):
    
    if args.filename:
        traj = Traj3D()
        lineList = [line.rstrip('\n') for line in open(args.filename)]
        seq = ''.join(lineList[1:])
    elif args.filename == "small": seq = "AAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCCAGTAAACGAAAAAACCGCCTGGGGAGGCGGTTTAGTCGAAGGTTAAGTCAG"

    list_individu = ListIndividu(nb_indiv,seq)

    list_individu.initial_list()
    
    
    print("ca part")

    for i in range(nb_gen):
        list_individu.select()


        list_individu.fusion_total()


        list_individu.mutate(0.1)
        
        best_indiv = list_individu.get_best()
        print("generation = ",i+1)
        print(best_indiv.fitness(seq))
        print(" ")
    
    print(best_indiv.Table)
    traj = Traj3D()
    traj.compute(seq, best_indiv)
    traj.draw("sample.png")

    

if __name__ == "__main__" :
    if args.nb_gen and args.nb_indiv:
        main(int(args.nb_gen),(int(args.nb_indiv)//4)*4) #la construction de notre algorithme impose un nb_indiv divisible par 4
    else:
        main()