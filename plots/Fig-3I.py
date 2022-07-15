# -*- coding: utf-8 -*-
# Time : 2021/9/22 13:52
# Author : 张志康

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from matplotlib.pyplot import MultipleLocator
import numpy as np

mese = pd.read_csv('Mesenchymal wikipathway_enrichment_result.csv')

####################################
basa = pd.read_csv("Basal wikipathway_enrichment_result.csv")

####################################
meta = pd.read_csv("Metabolic wikipathway_enrichment_result.csv")

####################################
pathway = pd.read_csv('list.txt',header=None,index_col=0)
mese = pd.merge(mese,pathway,left_on='Term_description',right_index=True,how='outer')
mese = mese[mese["Term_description"].isin(pathway.index)]
basa = pd.merge(basa,pathway,left_on='Term_description',right_index=True,how='outer')
basa = basa[basa["Term_description"].isin(pathway.index)]
meta = pd.merge(meta,pathway,left_on='Term_description',right_index=True,how='outer')
meta = meta[meta["Term_description"].isin(pathway.index)]
####################################################################################
pathway = list(pd.read_csv('list.txt',header=None,index_col=0).index)[::-1]
mese["Term_description"] = mese["Term_description"].astype("category")
mese["Term_description"].cat.reorder_categories(pathway, inplace=True)
mese.sort_values('Term_description', inplace=True)
######
basa["Term_description"] = basa["Term_description"].astype("category")
basa["Term_description"].cat.reorder_categories(pathway, inplace=True)
basa.sort_values('Term_description', inplace=True)
#######
meta["Term_description"] = meta["Term_description"].astype("category")
meta["Term_description"].cat.reorder_categories(pathway, inplace=True)
meta.sort_values('Term_description', inplace=True)

####################################
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['axes.unicode_minus'] = False
clist = ['#DEDFDF','#579EFA']
newcmp = LinearSegmentedColormap.from_list('chaos',clist)
fig = plt.figure(figsize=(10, 8))
#创建坐标轴，projection='rectilinear'笛卡尔坐标系
ax1 = fig.add_axes([0.45, 0.15, 0.25, 0.8])
sc =plt.scatter([0.1]*10,pathway,c=meta["Enrichment_score"],s=np.log2(-np.log10(meta["p-value"])+1.6)*130,cmap=newcmp)
plt.scatter([0.2]*10,pathway,c=basa["Enrichment_score"],s=np.log2(-np.log10(basa["p-value"])+1.6)*130,cmap=newcmp)
plt.scatter([0.3]*10,pathway,c=mese["Enrichment_score"],s=np.log2(-np.log10(mese["p-value"])+1.6)*130,cmap=newcmp)
plt.xlim([0.05,0.35])
plt.xticks([0.1,0.2,0.3],["Metabolic","Basal","Mesenchymal"],size=12)
plt.yticks(size=12)

######################################################################################
ax2 = fig.add_axes([0.75, 0.3, 0.035, 0.15])
ax2.axis('off')  # turn off all spines and ticks
unique_listhit = list(set(-np.log10(meta["p-value"]).dropna()) | set(-np.log10(basa["p-value"].dropna())) | set(-np.log10(mese["p-value"].dropna())))
siz = list(set(np.log2(-np.log10(meta["p-value"]).dropna()+1.6)) | set(np.log2(-np.log10(basa["p-value"].dropna())+1.6)) | set(np.log2(-np.log10(mese["p-value"].dropna())+1.6)))

# unique_listhit = [abs(i) for i in unique_listhit]
if len(unique_listhit) >= 3:
    x2 = [0] * 3
    # print(unique_listhit)
    num2 = np.percentile(unique_listhit,(100,50,10))
    num3 = np.percentile(siz,(100,50,10))
    # print(unique_listhit)
    print(num2)

else:
    x2 = [0] * len(unique_listhit)
    num2 = unique_listhit

s2 = [130*n for n in num3]
y2 = [4 * i for i in range(0, len(x2))]
ax2.set_ylim([-1.5, y2[-1] * 1.2])
ax2.scatter(x=x2, y=y2, s=s2, c='k', edgecolors='face')
for i, value in enumerate(num2):
       ax2.text(x=0.06, y=y2[i], s=round(value,2),
                va='center', ha='left', fontsize=10)

ax2.text(-0.1,  y2[-1] * 1.4, "-log$_{10}$(pvalue)", fontdict={'size': 12}, ha='left')
#########################################################################
cax = fig.add_axes([0.76, 0.6, 0.015, 0.20])
cb = fig.colorbar(sc, cax=cax)
cb.set_ticks(MultipleLocator(0.2))
ax1.text(0.37, 8, "Enrichment score",fontdict={'size': 12})
##########################################
plt.tight_layout()
plt.gcf().subplots_adjust()
plt.savefig("enrichment.pdf")
plt.savefig('enrichment.png',dpi=300)
plt.close()