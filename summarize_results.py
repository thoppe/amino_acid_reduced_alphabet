import glob, json, collections
import numpy as np

import compute_energy
df, L = compute_energy.load_interaction_matrix("MJ96")
published_results = collections.defaultdict(dict)

for parts in [2,4,5,8]:
    for key,scheme in compute_energy.published_schemes[parts].items():    
        B = compute_energy.sub_matrix(scheme,df)
        epsilon = compute_energy.compute_errors(scheme, df, L, B)
        published_results[parts][key] = epsilon
   

BT_FILES = glob.glob("results/BT_*")
MJ96_FILES = glob.glob("results/MJ96_*")
FILES = glob.glob("results/*")

data_time = collections.defaultdict(dict)
data_epsilon = collections.defaultdict(dict)

for f in FILES:
    with open(f) as FIN:
        js = json.load(FIN)
    n    = js["bead_target"]
    name = js["interaction_matrix"]

    data_time[name][n] = abs(js["computation_time"])
    data_epsilon[name][n] = js["epsilon"]

import pylab as plt
import seaborn as sns
sns.set(style="white", palette="muted")
fig, axes = plt.subplots(2, 1, figsize=(8, 11))
print axes
#sns.despine(left=True)

# Plot the published results
for parts in published_results:
    N = np.ones(len(published_results[parts]))*parts
    Y = [published_results[parts][key] for key in published_results[parts]]

    if parts == 5: text = "Published results"
    else: text = ""
    axes[1].scatter(N,Y,color='k',s=95,marker='o',alpha=.75,
                    label=text)



for name in ["MJ96","SJKG"]:
    N    = sorted(data_time[name].keys())
    TIME = [data_time[name][n] for n in N]
    EPSILON = [data_epsilon[name][n] for n in N]

    print N, TIME
    axes[0].semilogy(N, TIME, label=name)
    axes[1].semilogy(N, EPSILON, label=name)


axes[0].set_xticks(np.arange(1, 12+1, 1.0))
axes[1].set_xticks(np.arange(1, 12+1, 1.0))
axes[0].set_xlim(1.0,12)
axes[1].set_xlim(1.0,12)
   
axes[0].set_xlabel("Number of partitions")
axes[0].set_ylabel(r"$\ln t$")
axes[0].set_title("Computation time")

axes[1].set_xlabel("Number of partitions")
axes[1].set_ylabel(r"RMSE to target")
axes[1].set_title("Accuracy of optimal solution")

axes[0].legend(loc="best",fontsize=18)
axes[1].legend(loc="best",fontsize=18)
axes[1].set_ylim(ymax=1.0)

text = "Exact minimum for reduced amino acid representations"
plt.suptitle(text, fontsize=18)

plt.savefig("figures/computation_summary.png")
plt.show()



plt.close()


