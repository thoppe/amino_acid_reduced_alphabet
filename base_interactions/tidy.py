import numpy as np

np.set_printoptions(formatter={'float': '{: 0.2f}'.format})

f = "B.txt"
A = np.zeros((20,20))
row_n = 0
with open(f) as FIN:
    for line in FIN:
        line = line.strip()
        if line[0] != "#":
            missing = 20-len(line.split())
            terms = [0.0,]*missing + map(float,line.split())
            A[row_n] = terms
            row_n += 1

for i,j in zip(*np.triu_indices(20,k=1)):
    A[j,i] = A[i,j]

np.savetxt("tmp.txt",A,fmt="% .2f")
print A


