import numpy as np
import compute_energy
from stirling import stirling_set
import sys
from compute_energy import sub_matrix

import logging
logging.basicConfig(level=logging.DEBUG)

np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

# Use command line to input number of beads
bead_target = int(sys.argv[1])
L = compute_energy.residue_rows_letters

# Generate intermediate steps
starting_points = []
index_cutoff = 5
def find_starting_points(A,**kwargs): 
    if kwargs["index"]==index_cutoff:
        starting_points.append(A)
        return True
    return False

for index in stirling_set(L, early_break=find_starting_points): pass

import multiprocessing 
global_min = multiprocessing.Value('f', 1.0)

df = compute_energy.read_MJ_matrix()

target_set = compute_energy.five_bead_schemes["Cieplak_2001"]
target_error = compute_energy.compute_errors(target_set, df, 
                                             sub_matrix(target_set, df))

def branch_bound(scheme,**kwargs):

    # Always discard set if it has too many members
    if len(scheme)>bead_target: return True

    if scheme:
        B = compute_energy.sub_matrix(scheme, df)
        epsilon = compute_energy.compute_errors(scheme,df,B)
        if epsilon > global_min.value: 
            return True

    return False

def pretty_string(scheme):
    return ' '.join([''.join(x) for x in scheme])

def solve_scheme(start_pos=None,start_index=None,residue_letters=None):
    global global_min

    best_scheme = None

    if start_pos == None:
        scheme_iter = stirling_set(L, early_break=branch_bound)
    else:
        scheme_iter = stirling_set(L, 
                                   items=start_pos,
                                   index=start_index,
                                   early_break=branch_bound)

    for scheme in scheme_iter:

        B = compute_energy.sub_matrix(scheme, df)
        epsilon = compute_energy.compute_errors(scheme,df,B)
        if epsilon < global_min.value and len(scheme)==bead_target:

            with global_min.get_lock():
                global_min.value = epsilon

            vals = pretty_string(scheme), epsilon, target_error 
            logging.info("{} {:.4f} {:.4f}".format(*vals))
            best_scheme = scheme

    return best_scheme


P = multiprocessing.Pool()
PROCS = []

for pt in starting_points:
    kwargs = {"start_index":index_cutoff,
              "start_pos"  :pt,
              "residue_letters":L}
    sol = P.apply_async(solve_scheme, kwds=kwargs)
    PROCS.append(sol)

# Read the results from each proc
RESULTS = {}
for result in PROCS:
    scheme = result.get()
    if scheme:
        epsilon = compute_energy.compute_errors(scheme,df)
        RESULTS[epsilon] = scheme

min_epsilon = min(RESULTS.keys())
best_scheme = RESULTS[min_epsilon]
print

B = compute_energy.sub_matrix(best_scheme, df)
print pretty_string(best_scheme), min_epsilon       
print
print B



