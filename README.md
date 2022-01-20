# cPEPmatch


cPEPmatch is computational approach to rationally design cyclic peptides that potentially bind at desired regions of the interface of protein-protein complexes. The methodology is based on comparing the protein backbone structure of short peptide segments (epitopes) at the protein-protein interface with a collection of cyclic peptide backbone structures. A cyclic peptide that matches the backbone structure of the segment is used as a template for a binder by adapting the amino acid side chains to the side chains found in the target complex. 


Created by Brianda L. Santini - Physics Department T38, Technical University of Munich, Garching, Germany

Citation: Santini BL, Zacharias M. Rapid in silico Design of Potential Cyclic Peptide Binders Targeting Protein-Protein Interfaces. Front Chem. 2020  doi: 10.3389/fchem.2020.573259

## Environment

- Download and place this directory in your working directory (eg. inside /home/user/Documents/cPEPmatch/).
      
- The overall program should be ran under Anaconda environment (https://docs.anaconda.com/anaconda/install/linux/). 
      
      
- Modules to install: 

            Numpy and Scipy - pip install numpy scipy matplotlib ipython jupyter pandas sympy nose

            BioPython - pip install biopython (https://biopython.org/wiki/Download)

            VMD - conda install -c conda-forge vmd-python (https://vmd.robinbetz.com/)

            Modeller - This one is a bit tricky try: 

                  conda config --add channels salilab

                  conda install modeller

            You will be prompted after installation to edit a file to add your Modeller license key. Alternatively, set the KEY_MODELLER environment variable to your license key before you run 'conda install'. (https://salilab.org/modeller/9.17/release.html)


            NGL Viewer - pip install nglview



## Program Sections

### 1. Cyclic Peptide Database (Database folder: cyclo_database)

- There is no need to run this program constatntly unless you want to add a new cyclic peptide to the library. 
      
- The database folder contains all cyclic peptides, clean and renumbered. Name should be ‘pdb’-cp.pdb.
      
- There are two different lists of cyclic peptides, (1) contains only structures with standard amino acids, and (2) contains a list of all amino acids including those with staples or NON-standard amino acids. The reason for this division is that our mutation step with Modeller does not work if there are unrecognized structures in it. This is still to be fixed, but note that the matching does work for all of them.
      
- The output of this program is a database.pkl/ database_nonstandard.pkl binary file. This file is read by the Protein Matching program




### 2. Protein Matching (Program: cPEP-match.ipynb)

This program has three sections: 

a.  Interface Characterization - This selects the all the protein residues that are within an interface_cutoff of the GAG. It outputs an interface.pdb file which will be read by the next step.

b.  Motif construction -  This calculates the distance motifs for the interface residues contained in the interface.pdb file.

c.  Match + Mutate -  This section matches the interface motifs with those from the database.pkl file. The matches read and those that overlap with the GAG are deleted. Then, mutation + minimization of the adapted side chains is done using Modeller. 


## Running instructions:

- Before running this program, make sure to create a new directory inside your working folder (eg. inside /home/user/Documents/cPEPmatch/), the name for this directory should be the pdb name (eg. 1qqp). Place inside the clean protein-gag complex file pdb (eg. 1qqp.pdb) that you want to test. 
      
- Variables that need to be adjusted: 


            interface_cutoff (units are Angstroms, eg. 6-10)

            frmsd_threshold (units are Angstroms, eg. 0.5-1)

            protein1 residues (eg. ‘1 to 669’)

            gag residues (eg. ‘670 to 675’) 

            pdb name (eg. ‘1qqp’)

            my_location (eg, '/home/user/Documents/cPEPmatch/')

      

- The output of this program is a match_list.txt file, along with separate PDB files of all the matched and adapted cyclic peptides.


### Recent updates:

- Supports all motif size selections (4-6 work best).
- Allows a selection of non-consecutive amino acids (set "conscutive = False").
- The database program creator is now attached within the same program (if you want to run it, set "CREATE_CYCLOLIB=1")
