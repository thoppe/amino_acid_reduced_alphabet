import pandas as pd
import itertools
import numpy as np
from scipy import stats
#from numba import jit

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

'''
Amino acid propensities are taken from the publication:

@article{costantini2006amino,
  title={Amino acid propensities for secondary structures are influenced by the protein structural class},
  author={Costantini, Susan and Colonna, Giovanni and Facchiano, Angelo M},
  journal={Biochemical and biophysical research communications},
  volume={342},
  number={2},
  pages={441--451},
  year={2006},
  publisher={Elsevier}
}
'''

residue_propensities = {
    "A" : .0773,
    "C" : .0184,
    "D" : .0582,
    "E" : .0661,
    "F" : .0405,
    "G" : .0711,
    "H" : .0235,
    "I" : .0566,
    "K" : .0627,
    "L" : .0883,
    "M" : .0208,
    "N" : .0450,
    "P" : .0452,
    "Q" : .0394,
    "R" : .0503,
    "S" : .0613,
    "T" : .0553,
    "V" : .0691,
    "W" : .0151,
    "Y" : .0354,
}


def read_starting_comments(f, comment_char = "#"):
    comments = []
    with open(f) as FIN:
        for line in FIN:
            line = line.strip()
            if line[0] == comment_char:
                comments.append(line[1:])
            else: 
                break
    return comments

_valid_interaction_names = ["MJ96", "BT", "SJKG"]

def load_interaction_matrix(name):

    if name not in _valid_interaction_names:
        raise KeyError("Matrix {} not defined".format(matrix_name))

    if name == "MJ96":
        f_matrix = "base_interactions/MJ96.txt"

    if name == "BT":
        f_matrix = "base_interactions/BT.txt"

    if name == "SJKG":
        f_matrix = "base_interactions/SJKG.txt"

    comments = read_starting_comments(f_matrix)
    residue_rows = comments[-1].upper().split()
    residue_rows_letters = [residue_mapping[item] for item in residue_rows]

    df = read_interaction_matrix(f_matrix,comments,residue_rows_letters)

    return df, residue_rows_letters


def read_interaction_matrix(f_interaction,comments,residue_rows_letters):
    comments = read_starting_comments(f_interaction)

    df = pd.read_csv(f_interaction, sep='\s+', comment="#", 
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

def compute_errors(scheme,df, L, B=None, residue_weighted=False):
    ''' 
    This is supposed to match with Table 3 and this gets the right 
    rank order, but the magnitudes are not correct.
    '''

    A = np.array(df).ravel()

    if type(B) == type(None):
        B = sub_matrix(scheme, df)        
    
    ''' B is a submatrix and scheme is list of list of letters '''
    #print scheme, '\n', B

    scheme_mapping = {}
    for k, block in enumerate(scheme):
        for letter in block:
            scheme_mapping[letter] = k

    X = np.zeros((20,20))
    W = np.zeros((20,20))
    for i,j in zip(*np.triu_indices(20,k=0)):

        l1 = L[i]
        l2 = L[j]

        W[i,j] = residue_propensities[l1]*residue_propensities[l2]

        if l1 in scheme_mapping and l2 in scheme_mapping:
            bi = scheme_mapping[l1]
            bj  = scheme_mapping[l2]
            X[i,j] = X[j,i] = B[bi,bj]

        else:
            X[i,j] = X[j,i] = df[l1][l2]

    X = X.ravel()
    W = W.ravel()
    
    # Remove weight if not asked for
    if not residue_weighted:  W = np.ones(W.shape)
    squares = (A-X)**2
    
    epsilon = np.sqrt(np.average(squares,weights=W))
    return epsilon


'''
Five-bead schemes in Table 2 of Luthra et. al. 
'''

import collections
published_schemes = collections.defaultdict(dict)

published_schemes[5]["PAM_Koisol_2004"] = "AGTSNQDEHRKP,W,YF,MIVL,C"
published_schemes[5]["WAG_Koisol_2004"] = "AGTSNQDEHRKP,CV,IML,FY,W"

published_schemes[2]["Wang_Wang_1999"]  = "LFIMVWCY,HATGPRQSNEDK"
published_schemes[5]["Wang_Wang_1999"]  = "CMIFLYWV,AHT,GP,QNRSK,DE"

published_schemes[5]["Wang_Wang_2002"]  = "CMFI,LVWY,AGTS,NQDE,HPRK"
published_schemes[5]["Li_2003"]         = "CFYW,MLIV,G,PATS,NHQEDRK"
published_schemes[5]["Chemical_prop"]   = "IVL,FYWH,KRDE,GACS,TMNQP"
published_schemes[5]["Cieplak_2001"]    = "LFI,MVWCY,HA,TGPRQSNED,K"

published_schemes[2]["Cieplak_2001"]    = "LFIMVWCY,HATGPRQSNEDK"

published_schemes[4]["Shepherd_2007"]   = "LFI,MVWCY,HATGP,RQSNEDK"
published_schemes[8]["Shepherd_2007"]   = "LF,I,MVW,CY,HAT,GP,RQS,NEDK"

published_schemes[5]["Rakshit_2007"]   = "C,KR,MFILVWYA,DE,GTSNQHP"



for parts in published_schemes:
    for key,scheme in published_schemes[parts].items():
        published_schemes[parts][key] = [x for x in scheme.split(',')]

def pretty_string(scheme):
    return ' '.join([''.join(x) for x in scheme])

if __name__ == "__main__":

    df, L = load_interaction_matrix("MJ96")

    np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

    for key,scheme in published_schemes[5].items():
        B = sub_matrix(scheme, df)
        print key
        print pretty_string(scheme), '\n', B, '\n', 
        print compute_errors(scheme,df,L,B)
        print
        break

