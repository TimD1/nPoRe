import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as lines
plt.rcParams.update({'font.size': 32})

fig, axs = plt.subplots(1,1, figsize=(15,7))

callvcfs = ['clair3', 'clair3-hap', 'clair3-npore-hap']
colors = ['purple', 'red', 'green']

truthvcfs = ['Truth VCF', 'Standardized Truth VCF']
markers = ['+', '.']

labels = [patches.Patch(color=c, label=l) for c,l in zip(colors, callvcfs)] + \
[lines.Line2D([0], [0], marker=m, color='k', linestyle="",
label=l, markerfacecolor='k', markersize=15) 
for m,l in zip(markers, truthvcfs)]
plt.legend(handles=labels, loc='center')

plt.tight_layout()
plt.savefig(f'img/happy_acc_legend.png')
plt.close()
