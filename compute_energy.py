import pandas as pd
import itertools
import numpy as np
from scipy import stats

residue_mapping = {
  "ALA":'A',
  "CYS":'C',
  "ASP":'D',
  "GLU":'E',
  "PHE":'F',
  "GLY":'G',
  "HIS":'H',
  "ILE":'I',
  "LYS":'K',
  "LEU":'L',
  "MET":'M',
  "ASN":'N',
  "PRO":'P',
  "GLN":'Q',
  "ARG":'R',
  "SER":'S',
  "THR":'T',
  "VAL":'V',
  "TRP":'W',
  "TYR":'Y'}

def read_starting_comments(f, comment_char = "#"):
    comments = []
    with open(f_MJ) as FIN:
        for line in FIN:
            line = line.strip()
            if line[0] == comment_char:
                comments.append(line[1:])
            else: 
                break
    return comments
    
f_MJ = "base_interactions/MJ.txt"
comments = read_starting_comments(f_MJ)
residue_rows = comments[-1].upper().split()
residue_rows_letters = [residue_mapping[item] for item in residue_rows]


df = pd.read_csv(f_MJ, sep='\s+', comment="#", 
                 skiprows=len(comments), header=None)
df.columns = residue_rows_letters
df.index   = residue_rows_letters

def sub_block(rowL,colL):
    k1,k2 = map(len, [rowL, colL])
    A = np.zeros((k1,k2))
    
    for i,j in itertools.product(range(k1),range(k2)):
        A[i,j] = df[rowL[i]][colL[j]]

    return A

def representative_weight(rowL, colL):
    return sub_block(rowL, colL).mean()

def sub_matrix(letter_blocks):
    k = len(letter_blocks)
    A = np.zeros((k,k))
    
    for i,j in itertools.product(range(k),repeat=2):
        A[i,j] = representative_weight(letter_blocks[i], letter_blocks[j])
    return A

def compute_errors(scheme,B):
    ''' 
    This is supposed to match with Table 3 and this gets the right 
    rank order, but the magnitudes are not correct.
    '''
    
    ''' B is a submatrix and scheme is list of list of letters '''
    #print scheme, '\n', B
    reduced, exact = [],[]
    X = np.zeros((20,20))
    for i,j in itertools.product(range(len(scheme)),repeat=2):
        r = B[i,j]
        for l1,l2 in itertools.product(scheme[i], scheme[j]):
            reduced.append(r)
            exact.append(df[l1][l2])

    exact = np.array(exact)
    reduced = np.array(reduced)
    rms_error = np.sqrt(((exact-reduced)**2).mean())
    _,_,corr,_,_, = stats.linregress(exact,reduced)
    return rms_error, corr

                                    

'''
Five-bead schemes in Table 2 of Luthra et. al. 
'''

five_bead_schemes = {}
five_bead_schemes["PAM_Koisol_2004"] = "AGTSNQDEHRKP,W,YF,MIVL,C"
five_bead_schemes["WAG_Koisol_2004"] = "AGTSNQDEHRKP,CV,IML,FY,W"
five_bead_schemes["Wang_Wang_1999"]  = "CMIFLYWV,AHT,GP,QNRSK,DE"
five_bead_schemes["Wang_Wang_2002"]  = "CMFI,LVWY,AGTS,NQDE,HPRK"
five_bead_schemes["Li_2003"]         = "CFYW,MLIV,G,PATS,NHQEDRK"
five_bead_schemes["Chemical_prop"]   = "IVL,FYWH,KRDE,GACS,TMNQP"
five_bead_schemes["Cieplak_2001"]    = "LFI,MVWCY,HA,TGPRQSNED,K"

for key,scheme in five_bead_schemes.items():
    five_bead_schemes[key] = [x for x in scheme.split(',')]

for key,scheme in five_bead_schemes.items():
    B = sub_matrix(scheme)
    print key, '\n', B, '\n', compute_errors(scheme,B)
    print

