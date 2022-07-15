import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns

import matplotlib.patches as mpatches
import matplotlib


e = list(pd.read_csv('ecm.txt',sep='\t',header=None,index_col=0).index)
i =list( pd.read_csv('immune.txt',sep='\t',header=None,index_col=0).index)
k = list(pd.read_csv('kinase.txt',sep='\t',header=None,index_col=0).index)

lists = k+i+e

def deal(ks,key):
    tmp = pd.DataFrame()
    for k in ks:
        dk = pd.read_csv(k,sep='\t')[['geneset','logFC','avgExp','t','pval','FDR']]
        tmp = pd.concat([tmp,dk])
    tmp['geneset'] =  tmp['geneset'].apply(lambda x:x.split('(')[0].capitalize().replace('4f','4F'))
    tmp = tmp[tmp['geneset'].isin(lists)]
    tmp['abs'] = abs(tmp['t'])
    tmp = tmp.sort_values(by='abs',ascending=False)
    tmp = tmp.drop_duplicates(subset='geneset',keep='first')
    tmp.set_index('geneset',inplace=True)
    ens = tmp[['t']]
    ens.columns = [key]
    fdr = tmp[['FDR']]
    fdr.columns = [key]
    return ens,fdr

    
ks = glob.glob('Kinase*/*/allexp_genesets_GSVA_score.xls')
    
iss = glob.glob('Immu*/*/allexp_genesets_GSVA_score.xls')

es = glob.glob('ECM*/*/allexp_genesets_GSVA_score.xls')

dks,fks = deal(ks,'Metabolic')

dis,fis = deal(iss,'Basal')

des,fes = deal(es,'Mesenchymal')


data = pd.concat([dks,dis,des],axis=1).T[lists].T
pvalue = pd.concat([fks,fis,fes],axis=1).T[lists].T
#pps = pd.concat([pks,pis,pes],axis=1).T[df0].T

anno = pd.DataFrame()
for index in data.index:
    for column in data.columns:
        cor = data.loc[index, column]
        p = pvalue.loc[index, column]
        new_p = sns.utils.sig_stars(p)
        anno.loc[index, column] = f"{new_p}" if new_p else ''
        

norm = matplotlib.colors.Normalize(0,1)
colors = [[norm(0), '#000080'],
         # [norm(0.25), '#00BFFF'],
          [norm(0.5), "#FFFFFF"],
          [norm( 0.75), '#FF8C00'],
          [norm( 1), '#FFFF00']]

cc = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)


c11 = {'ECM':'#00468B'}

c11.update({'Kinase':'#ED0000'})

c11.update({'Immune':'#42B540'})

plt.rcParams['font.sans-serif'] = ['Arial'] 
#matplotlib.rcParams['font.family'] = "sans-serif"



dc = pd.DataFrame(c11,index=[0]).T
data.index.name = ''        
height = max(6, int(data.shape[0] * 0.45))
width = max(8, int(data.shape[1] * .7))
fig, ax = plt.subplots(figsize=(width, height))
# 将每个module增加上对应的基因数量
ax = sns.heatmap(data, cmap=cc, center=0, annot=anno, square=True, 
                 fmt='s',cbar_kws={"shrink": 0.6,'use_gridspec': False, 'location': 'left','pad':0.2,'label':'t_value'})
#ax.plot([-1,-0.05],[0,9],color='k',lw=2)
ax.tick_params(right=True, labelright=True,left=False,labelleft=False,rotation=0)
ax.set_xticklabels(['Metabolic','Basal','Mesenchymal'],rotation=90)
ax1 = ax.inset_axes([-0.15,0.6,0.15,0.42])
ax1.bar(0.05,0.4,width=5,color='#ED0000')
ax1.axis('off')
ax2 = ax.inset_axes([-0.15,0.2,0.15,0.44])
ax2.bar(0.05,0.4,width=5,color='#42B540')
ax2.axis('off')
ax3 = ax.inset_axes([-0.15,0.,0.15,0.352])
ax3.bar(0.05,0.4,width=5,color='#00468B')
ax3.axis('off')

patches = [mpatches.Patch(color=['#ED0000', '#42B540', '#00468B'][i], label="{:s}".format(['Metabolic','Basal', 'Mesenchymal'][i]) ) for i in range(len(list(c11.values())))]
second = plt.legend(handles=patches, bbox_to_anchor = (-2, 1), loc='upper left', ncol=1, title='Subtype', frameon=False)
ax = plt.gca().add_artist(second)
plt.savefig('heatmap-fdr.png',dpi=300,bbox_inches='tight')
plt.savefig('heatmap-fdr.pdf',bbox_inches='tight')
plt.close()


