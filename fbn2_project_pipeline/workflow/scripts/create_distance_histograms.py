import numpy as np
from matplotlib import pyplot as plt
import tracklib as tl
import pandas as pd
import scipy.stats
from copy import deepcopy

import matplotlib
matplotlib.use("Agg")

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams.update({'font.size': 15})


data = tl.io.load.csv(snakemake.input.tagged_set_file, 
                    ['id', 't', 'x', 'y', 'z', 'x2', 'y2', 'z2','date','movie_index'],
                      skip_header=1)

data.makeSelection()

dists = np.concatenate([traj.relative().abs()[:].flatten() for traj in data])
np.save(snakemake.output.histogram_numpy,dists)

fig=plt.figure()
_ , _, h = plt.hist(dists, 
                    bins='auto', 
                    density=True, 
                    cumulative=True,
                    alpha=0.85, 
                    label=f" {snakemake.wildcards.sample_group}  (median = {int(np.nanmedian(dists)*1000)} nm)",
                    histtype='step',
                    lw=3,
                    )
plt.xlim([0,1.5])
plt.ylim([0,1])
plt.xlabel('Distance ($\mu m$)')
plt.ylabel('Cumulative Probability')
plt.legend()
plt.savefig(snakemake.output.histogram_cdf,bbox_inches="tight")
plt.close(fig)


fig=plt.figure()
_ , _, h = plt.hist(dists, 
                    bins=200, 
                    density=True, 
                    cumulative=False,
                    alpha=0.85, 
                    label=f" {snakemake.wildcards.sample_group}  (median = {int(np.nanmedian(dists)*1000)} nm)",
                    histtype='step',
                    lw=3,
                    )
plt.xlim([0,1.5])
plt.ylim([0,3.6])
plt.xlabel('Distance ($\mu m$)')
plt.ylabel('Probability density')
plt.legend()
plt.savefig(snakemake.output.histogram_pdf,bbox_inches="tight")
plt.close(fig)

