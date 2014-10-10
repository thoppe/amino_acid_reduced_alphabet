import numpy as np
from stirling import stirling_set
import compute_energy
from compute_energy import sub_matrix
import logging, argparse, time

desc = ''' TBW '''

parser = argparse.ArgumentParser(description=desc)
parser.add_argument('interaction_matrix', type=str,
                    help="Interaction matrix, must be one of MJ96 BT SJKG")
parser.add_argument('--bead_target', '-b', type=int, default=3,
                    help="Number of input beads")
cargs = vars(parser.parse_args())

# Start the logger
logging.basicConfig(level=logging.INFO)
np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

# Measure the start time
start_time = time.time()

df, L = compute_energy.load_interaction_matrix(cargs["interaction_matrix"])

#########################################################################

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
global_min = multiprocessing.Value('f', 20**2)


target_set = compute_energy.five_bead_schemes["Cieplak_2001"]
target_error = compute_energy.compute_errors(target_set, df, L,
                                             sub_matrix(target_set, df))

def branch_bound(scheme,**kwargs):

    # Always discard set if it has too many members
    if len(scheme)>cargs["bead_target"]: return True

    if scheme:
        B = compute_energy.sub_matrix(scheme, df)
        epsilon = compute_energy.compute_errors(scheme,df,L,B)
        if epsilon > global_min.value: 
            return True

    return False

def pretty_string(scheme):
    return ' '.join([''.join(sorted(x)) for x in scheme])

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
        epsilon = compute_energy.compute_errors(scheme,df,L,B)
        if epsilon < global_min.value and len(scheme)==cargs["bead_target"]:

            with global_min.get_lock():
                # Check again under the lock to be sure
                if epsilon < global_min.value:
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
        epsilon = compute_energy.compute_errors(scheme,df,L)
        RESULTS[epsilon] = scheme

min_epsilon = min(RESULTS.keys())
best_scheme = RESULTS[min_epsilon]
print

B = compute_energy.sub_matrix(best_scheme, df)
print pretty_string(best_scheme), min_epsilon       
print
print B

data = cargs.copy()
data["scheme"] = pretty_string(best_scheme)
data["epsilon"] = min_epsilon
data["bead_interaction_matrix"] = B.tolist()
data["computation_time"] = start_time - time.time()

import json
np.set_printoptions(formatter={'float': '{: 0.6f}'.format})

f_results = "results/{interaction_matrix}_{bead_target}.json"
with open(f_results.format(**cargs), 'w') as FOUT:
    FOUT.write(json.dumps(data,indent=4))
    





