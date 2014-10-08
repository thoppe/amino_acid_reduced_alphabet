import pandas as pd
import itertools
import numpy as np
from scipy import stats

'''
Computes the interaction potential from the MJ matrix and a reduced
alphabet of amino acids.

@article{luthra2007method,
  title={A method for computing the inter-residue interaction 
         potentials for reduced amino acid alphabet},
  author={Luthra, Abhinav and Jha, Anupam Nath and 
          Ananthasuresh, GK and Vishveswara, Saraswathi},
  journal={Journal of biosciences},
  volume={32},
  number={1},
  pages={883--889},
  year={2007},
  publisher={Springer}
}
'''

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

def read_starting_comments(f_MJ, comment_char = "#"):
    comments = []
    with open(f_MJ) as FIN:
        for line in FIN:
            line = line.strip()
            if line[0] == comment_char:
                comments.append(line[1:])
            else: 
                break
    return comments

def read_MJ_matrix(f_MJ = "base_interactions/MJ.txt"):

    df = pd.read_csv(f_MJ, sep='\s+', comment="#", 
                     skiprows=len(comments), header=None)
    df.columns = residue_rows_letters
    df.index   = residue_rows_letters

    return df

def sub_block(rowL,colL, df):
    k1,k2 = map(len, [rowL, colL])
    A = np.zeros((k1,k2))
    
    for i,j in itertools.product(range(k1),range(k2)):
        A[i,j] = df[rowL[i]][colL[j]]

    return A

def representative_weight(rowL, colL, df):
    return sub_block(rowL, colL, df).mean()

def sub_matrix(letter_blocks,df):
    k = len(letter_blocks)
    A = np.zeros((k,k))
    
    for i,j in itertools.product(range(k),repeat=2):
        w = representative_weight(letter_blocks[i], letter_blocks[j],df)
        A[i,j] = w
    return A

def compute_errors(scheme,df, B=None):
    ''' 
    This is supposed to match with Table 3 and this gets the right 
    rank order, but the magnitudes are not correct.
    '''

    if B == None:
        B = sub_matrix(scheme, df)
    
    ''' B is a submatrix and scheme is list of list of letters '''
    #print scheme, '\n', B
    reduced, exact = [],[]

    scheme_mapping = {}
    for k, block in enumerate(scheme):
        for letter in block:
            scheme_mapping[letter] = k

    # We may not be passing a full scheme here
    letter_subset = scheme_mapping.keys()

    for l1,l2 in itertools.product(letter_subset,repeat=2):
        if l1>=l2:
            i,j = scheme_mapping[l1], scheme_mapping[l2]
            exact.append( df[l1][l2] )
            reduced.append( B[i,j] )

    exact = np.array(exact)
    reduced = np.array(reduced)
    total_L2_error = ((exact-reduced)**2).sum()
    return total_L2_error
    #rms_error = np.sqrt(((exact-reduced)**2).mean())
    #_,_,corr,_,_, = stats.linregress(exact,reduced)
    #return rms_error, corr


comments = read_starting_comments("base_interactions/MJ.txt")
residue_rows = comments[-1].upper().split()
residue_rows_letters = [residue_mapping[item] for item in residue_rows]
                                   

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

if __name__ == "__main__":
    df = read_MJ_matrix()

    for key,scheme in five_bead_schemes.items():
        B = sub_matrix(scheme, df)
        print key, '\n', B, '\n', compute_errors(scheme,df)
        print

