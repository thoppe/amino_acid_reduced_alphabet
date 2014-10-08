import numpy as np
import compute_energy
from stirling import stirling_set
import sys

from compute_energy import sub_matrix

np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

# Use command line to input number of beads
bead_target = int(sys.argv[1])

L = compute_energy.residue_rows_letters
df = compute_energy.read_MJ_matrix()

global_min = np.inf

target_set = compute_energy.five_bead_schemes["Cieplak_2001"]
target_error = compute_energy.compute_errors(target_set, df, 
                                             sub_matrix(target_set, df))

def branch_bound(scheme):

    # Always discard set if it has too many members
    if len(scheme)>bead_target: return True
    
    if scheme:
        B = compute_energy.sub_matrix(scheme, df)
        epsilon = compute_energy.compute_errors(scheme,df,B)
        if epsilon > global_min: 
            return True

    return False

def pretty_string(scheme):
    return ' '.join([''.join(x) for x in scheme])

best_scheme = None

for scheme in stirling_set(L, early_break=branch_bound):
    B = compute_energy.sub_matrix(scheme, df)
    epsilon = compute_energy.compute_errors(scheme,df,B)
    if epsilon < global_min and len(scheme)==bead_target:
        global_min = epsilon
        vals = pretty_string(scheme), epsilon, target_error 
        print "{} {:.4f} {:.4f}".format(*vals)
        best_scheme = scheme


        
B = compute_energy.sub_matrix(best_scheme, df)
print
print B



