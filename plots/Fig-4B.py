# -*- coding: utf-8 -*-
# Time : 2021/11/19 10:56
# Author : 张志康

import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
from sklearn import preprocessing
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib as mpl
import matplotlib.cm as cm

import matplotlib
import matplotlib.colors
import seaborn as sns
import numpy as np
import matplotlib.patches as mpatches

def col(cc1,cc3,cc2):
    norm = matplotlib.colors.Normalize(-1,1)
    colors = [[norm(-1.0), cc1],
              #[norm(-0.85), cc3],
              [norm(0), "white"],
              [norm( 1.0), cc2]]
    
    return matplotlib.colors.LinearSegmentedColormap.from_list("", colors)

fpkm_cmap =  col("#003572",'#4B5E8F', "#FF2021")


#print(Mesenchymal)
###########################################
new = pd.read_csv('order.txt',sep='\t',index_col=0)



fpkm = pd.read_csv("data.txt",sep='\t',index_col=0).iloc[:-2,:]

lst = []
for i in new.index:
    if i in fpkm.columns:
        lst.append(i)

fpkm = fpkm[lst]
#print(PROM1_TNC_protein_p)

kegg = pd.read_csv('all_samples_ES.xls',sep='\t')
es1 = kegg[kegg.Term.str.contains('MAPK')].set_index('Term')
es1.columns = [n.split('_')[1] for n in es1.columns]
es1 = es1[lst]
es1 = es1.T
es1.columns = ['MAPK signature score']

new_red=[sns.color_palette("Blues", 200)[100:174]][0]

es2 = es1.sort_values(by='MAPK signature score')

c2 = {c:n for c,n in zip(es2.index,new_red)}


############################################
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['axes.unicode_minus'] = False

genes = list(fpkm.index)

fig = plt.figure(figsize=(12, 4))
ax = fig.add_axes([0.1,0.1,0.8,0.8])
ax.axis("off")
plt.xlim([0,0.5])
plt.ylim([0,1])
x= 0.075
y=0.025

ax1_ = []
ax2_ = []

for i in range(fpkm.shape[0]):
    #ax.plot([0.73, 0.73], [0.6 - i * x - i * y, 0.6 - (i + 1) * x - i * y], color='k')
    #ax.text(0.67, 0.605 - i * x - i * y - 0.05, ["mRNA","Protein"][i])
    ax1 = fig.add_axes([0.2, 0.3 - i * 0.08, 0.6, 0.09])
    ax1.text(1.01,0.35,genes[i])
    #ax1.axis("off")
    ax1_.append(ax1)

for i in range(fpkm.shape[0]):

    fpkm_line = fpkm.loc[genes[i]]
    order1 = fpkm_line.index
    fpkm_line_Z = pd.Series(preprocessing.scale(fpkm_line), index=fpkm_line.index)
    #fpkm_cmap = LinearSegmentedColormap.from_list('chaos', ["#003572", "white", "#FF2021"])
    fpkm_norm = mpl.colors.Normalize(vmin=min(fpkm_line_Z), vmax=max(fpkm_line_Z))
    fpkm_mapper = cm.ScalarMappable(norm=fpkm_norm, cmap=fpkm_cmap)
    #print(i)
    for j in range(len(order1)):
        
        ax1_[i].add_patch(matplotlib.patches.Rectangle(xy=(j / 74, 0), width=2 / 74, height=1,facecolor=fpkm_mapper.to_rgba(fpkm_line_Z[order1[j]]), edgecolor='#F1F1F1',linewidth=1))
    ax1_[i].axis("off")

im = ax.imshow(np.arange(1).reshape((1, 1)),cmap=fpkm_cmap)
im.remove()
c = plt.colorbar(im, cax = fig.add_axes([0.67, -0.28, 0.09, 0.06]),orientation="horizontal")
c.ax.set_xticklabels("")
c.ax.xaxis.set_tick_params(size=0)
c.outline.set_visible(False)
c.ax.set_xticklabels(['row min','','row max'])


ax2 = fig.add_axes([0.2, 0.4, 0.6, 0.05])
order = es1.index
esz = pd.Series(preprocessing.scale(es1['MAPK signature score']), index=es1.index)
es_norm = mpl.colors.Normalize(vmin=min(esz), vmax=max(esz))
#es_mapper = cm.ScalarMappable(norm=es_norm, cmap=es_cmap)
for j in range(len(order)):
    ax2.add_patch(matplotlib.patches.Rectangle(xy=(j / 74, 0), width=2 / 74, height=2,color=c2[order[j]]))
ax2.text(1.01,0.35,'MAPK signature score')
ax2.axis('off')

ax3 = fig.add_axes([0.2, 0.45, 0.6, 0.05])
cs = ['#ED0000']*19+ ['#42B540']*31+['#00468B']*24
for j in range(len(order)):
    ax3.add_patch(matplotlib.patches.Rectangle(xy=(j / 74, 0), width=2 / 74, height=2,color=cs[j]))
ax3.text(1.01,0.35,'Group')
ax3.axis('off')

patches = [mpatches.Patch(color=[ '#ED0000', '#42B540','#00468B'][i], label="{:s}".format(['Metabolic','Basal','Mesenchymal'][i]) ) for i in range(3)]
second = plt.legend(handles=patches, bbox_to_anchor = (0, -12.8), loc='upper left', ncol=1, title='Group', frameon=False)
ax10 = plt.gca().add_artist(second)
second._legend_box.align = "left"

patches = [mpatches.Patch(color=[new_red[0], new_red[50],new_red[-1]][i], label="{:s}".format(['0.0124','0.050','0.0787'][i])) for i in range(3)]
second = plt.legend(handles=patches, bbox_to_anchor = (0.38, -12.8), loc='upper left', ncol=1, title='MAPK signature score', frameon=False)
ax3 = plt.gca().add_artist(second)
second._legend_box.align = "left"

plt.savefig("heatmap-4B.pdf",bbox_inches='tight')
plt.savefig("heatmap-4B.png",dpi=300,bbox_inches='tight')

