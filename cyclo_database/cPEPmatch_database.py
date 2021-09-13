#!/usr/bin/env python
# coding: utf-8

# ## Cyclo Library

# Characterization the cyclic peptide motifs.

# In[1]:


import sys,os,re, glob
import operator #set of efficient functions
import pickle #store data on disc
import time
import warnings
import numpy as np
import shutil #import copyfile
import csv

from scipy import spatial
from glob import glob
from collections import defaultdict
from itertools import product
from sys import exit


import Bio
from Bio.PDB import PDBParser
from Bio.PDB.parse_pdb_header import parse_pdb_header

import vmd
from vmd import molecule, atomsel
cwd = os.getcwd()


# In[2]:


CREATE_CYCLOLIB=1
motif_size=4


# In[3]:


def int_list(lyst, sep=' ', quote=None):
    """Return the list of numbers in 'lyst' as a string separated by 'sep'."""
    f = str if quote is None else lambda i: quote+str(i)+quote
    return sep.join(map(f, lyst))

def cyclo_distance_matrix(nres, R, fname):
    M = dict()
    
    for i,res in enumerate(R[:-nres+1]):
        #print(res, res+nres, R[i:i+nres])
        resids = R[i:i+nres]
        if not np.all(np.diff(resids)==1): continue
        rlist = int_list(resids, quote='"')
        #if not rnbrs: continue
        QA = '''(chain M and name CA and (altloc "A" or altloc "") and resid %(rlist)s)'''
        QA = ' '.join(str(QA).split())
        
        resid_query = atomsel(QA % dict(rlist=rlist))
        r = np.array(resid_query.centerperresidue())
        
        #
        # Ignore sequences that are either near the ends of a chain, or are near ends of a broken chain
        #
        if len(r) != nres: continue
        assert len(r) == nres, (len(r), r, nres, res, rnbrs, resid_query, fname)
        
        dr = spatial.distance_matrix(r,r)
        m = dr[np.triu_indices(nres,k=2)]
        M[resids[0]] = dict(
            m=m,
#            r=r,
            resid=resids,
        )
    return M



def cyclo_motifs(molid):
    C = atomsel('chain M and protein and backbone and (altloc "A" or altloc "") ')
    C_res = list(sorted(set(C.resid)))
    MC = cyclo_distance_matrix(motif_size, C_res, cyclo)
    
    return (MC)


# In[7]:


if CREATE_CYCLOLIB:
    os.chdir(cwd)

    parser = PDBParser()
    cyclo_mtfs = []
    
    #All cyclic peptides, including those with non-standard amino acids.
    #lib = ['2c9t', '1cnl', '1tk2', '5eoc', '3zwz', '4z0f', '4zjx', '4ib5', '3oe0', '6c1r', '1urc', '4w4z', 
          # '1ebp', '3pqz', '3avb', '4n7y', '5xco', '4k8y', '1jbl', '3p8f', '4i80', '5afg', '3wbn', '1npo', 
          # '2mgo', '3wmg', '5kgn', '5b4w', '5lso', '5xn3', '1sld', '5bxo', '4zqw', '4ktu', '4kts', '4twt', 
          # '3qn7', '3av9', '3ava', '3avi', '3wnf', '1sle', '1vwb', '5n99', '4os1', '4mnw', '4jk5', '4gly', 
          # '4x1n', '5glh', '2m2g', '2lwv', '6pip', '2lwt', '1foz', '2ajw', '5lff', '2ox2', '6axi', '2lwu', 
          # '3wne', '2ak0', '2lws', '1qx9', '6pin', '2n07', '1l5g', '6pio', '2jrw', '6awm', '5vav', '6awk' ]
    
    #Only cyclic peptides with standard amino acids.
    lib = ['1cnl', '1ebp', '1foz', '1jbl', '1l5g', '1npo', '1qx9', '1sld', '1sle', '1urc', '1vwb', '2ajw', 
           '2ak0', '2c9t', '2jrw', '2lws', '2lwt', '2lwu', '2lwv', '2m2g', '2mgo', '2n07', '2ox2', '3av9', 
           '3ava', '3avb', '3avi', '3p8f', '3wbn', '3wne', '3wnf', '3zwz', '4gly', '4ib5', '4jk5', 
           '4k8y', '4kts', '4ktu', '4mnw', '4mq9', '4twt', '4x1n', '5bxo', '5eoc', '5glh', '5lff', '5lso', 
           '5vav', '5xco', '5xn3', '6awk', '6awm', '6axi', '6pin', '6pio', '6pip' ]

   
    
    for c in lib:
        motifs = dict()
        cyclo = ('%s-cp.pdb' % c)
        molid = molecule.load("pdb", cyclo )
        c_mtfs = cyclo_motifs(molid)
        if c_mtfs is not None:
            #save_interacting_sidechains(os.path.splitext(f)[0], isc)
            motifs.update({c: c_mtfs})
            cyclo_mtfs.append(motifs)
            
            
database = open("database.pkl", "wb")
pickle.dump(cyclo_mtfs, database)
database.close()


# In[6]:


if CREATE_CYCLOLIB:
    os.chdir(cwd)

    parser = PDBParser()
    cyclo_mtfs = []
    
    #All cyclic peptides, including those with non-standard amino acids.
    lib = ['2c9t', '1cnl', '1tk2', '5eoc', '3zwz', '4z0f', '4zjx', '4ib5', '3oe0', '6c1r', '1urc', '4w4z', 
           '1ebp', '3pqz', '3avb', '4n7y', '5xco', '4k8y', '1jbl', '3p8f', '4i80', '5afg', '3wbn', '1npo', 
           '2mgo', '3wmg', '5kgn', '5b4w', '5lso', '5xn3', '1sld', '5bxo', '4zqw', '4ktu', '4kts', '4twt', 
           '3qn7', '3av9', '3ava', '3avi', '3wnf', '1sle', '1vwb', '5n99', '4os1', '4mnw', '4jk5', '4gly', 
           '4x1n', '5glh', '2m2g', '2lwv', '6pip', '2lwt', '1foz', '2ajw', '5lff', '2ox2', '6axi', '2lwu', 
           '3wne', '2ak0', '2lws', '1qx9', '6pin', '2n07', '1l5g', '6pio', '2jrw', '6awm', '5vav', '6awk' ]

    
    for c in lib:
        motifs = dict()
        cyclo = ('%s-cp.pdb' % c)
        molid = molecule.load("pdb", cyclo )
        c_mtfs = cyclo_motifs(molid)
        if c_mtfs is not None:
            #save_interacting_sidechains(os.path.splitext(f)[0], isc)
            motifs.update({c: c_mtfs})
            cyclo_mtfs.append(motifs)
            
            
database = open("database_nonstandard.pkl", "wb")
pickle.dump(cyclo_mtfs, database)
database.close()


# In[ ]:




