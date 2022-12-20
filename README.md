# groupe7 EI-Algorithmique Génétique

## Membres
Badre Mhiouah, Paul Casteras, Baptiste Carbillet, Guillaume Malle

## Contexte
Les séquences d'ADN sont aujourd'hui très largement étudiées à travers leur séquence textuelle (succession de A, C, G et T), il est aussi très instructif de les étudier à partir de leur trajectoire tri-dimensionnelle. En 1993, des biopysiciens ont établi un modèle de conformation 3D qui permet de transformer une suite de nucléotides (sous forme de lettres) en une trajectoire tri-dimensionnelle. Dès lors, il est possible de représenter toute séquence textuelle d'ADN en une trajectoire 3D.
si on observe un chromosome bactérien (longue séquence d'ADN constituant une bactérie) ou un plasmide (petite séquence présente au sein des bactéries), on s'aperçoit que ce chromosome ou ce plasmide est circulaire, i.e. les deux extrémités ont été "collées" l'une à l'autre. Le modèle pré-cité ne rend pas compte de ce phénomène lorsque l'on représente la trajectoire 3D d'un chromosome bactérien ou d'un plasmide.

## Projet
L'objectif de ce projet est de modifier le modèle de conformation 3D donné afin de rendre un plasmide circulaire. Pour cela, un algorithme génétique sera développé en Python et structuré en classes (programmation orientée objet).

## How to use
- Download the whole repository
- Install *mathutils*, *numpy* and *matplotlip* : run `pip install mathutils numpy matplotlib`
- Run Main.py with Python : `python Main.py` or `python3 Main.py`
    - There are 2 DNA sequence files named `plasmid_8k.fasta` and `plasmid_180k.fasta`.
    - Set the DNA sequence with argument `--filename` (*default = plasmid_8k.fasta*).
    - Set the number of indivuals in the population with argument `--nb_indiv`.
    - Set the number of generations with argument `--nb_gen`.


## Structure du code
Il y a 3 classes principales :
- 'Individu' : Permet de représenter un individu de la population. Un individu est une table de rotation, i.e un dictionnaire qui à chaque dinucléotide (paire de nucléotides A-A, A-C, ...) associe trois valeurs d'angles de rotation.
- 'Traj3D' : Permet de calculer et dessiner une trajectoire dans l'espace à partir d'une séquence d'ADN et d'une table de rotation.
- 'ListIndividu' : Permet de représenter la population entière, sous forme d'une liste d'individu. C'est cette liste qu'on modifie au fur et à mesure de l'évolution de la population.

Le fichier Main.py organise et structure les différentes parties du code (initalise une population, lance la sélection d'une partie de cette population, puis la reproduction et enfin la mutation). C'est depuis Main.py qu'on peut lancer des simulations.

