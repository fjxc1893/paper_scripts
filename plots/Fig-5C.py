import pandas as pd
from oebio.app import read
from oebio.plot import base
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
sns.set(color_codes=True)


def col(cc1,cc2):
    import matplotlib
    norm = matplotlib.colors.Normalize(-1,1)
    colors = [[norm(-1.0), cc1],
              [norm(0), "white"],
              [norm( 1.0), cc2]]
    
    return matplotlib.colors.LinearSegmentedColormap.from_list("", colors)

def deal(f,name,out,):
    df = read(f)
    df['geneset'] = df['geneset'].apply(lambda x: x.split(name)[-1].capitalize().replace('_',' '))
    df.to_csv(f'{out}.txt',sep='\t',index=False)
    

####################################
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['axes.unicode_minus'] = False
#clustc = sns.color_palette("seismic", 10000)[1000:9000]
#clustc = sns.color_palette("seismic", 10000)[3000:7000]
clustc = col('#0000FF','#FF8C00')

df = pd.read_csv('data.txt',sep='\t')
df1 = df[~df['name'].str.contains(remov,regex=True)]
df = df1.set_index('name')
new = pd.read_csv('order.txt',sep='\t',index_col=0)
lst = [i for i in new.index if i in df.columns]

color10_clu3 = ['#ED0000','#42B540',  '#00468B']
c10 = {c:n for c,n in zip(list(new.drop_duplicates(subset=['Group'])['Group']),color10_clu3)}

df1 = df[lst]
df1.index.name=''
cc = new.Group.map(c10).to_frame()

g = sns.clustermap(df1,col_colors=cc,z_score=0,figsize=(12, 6),cbar_pos=(0.065, .2, .02, .2),
                   cbar_kws={'label':''},cmap=clustc,center=0,col_cluster=False,
                   row_cluster=False,yticklabels=True,xticklabels=False)
g.ax_col_colors.yaxis.set_tick_params(size=0.5)
g.ax_col_colors.set_yticklabels(cc,size=9)
g.ax_col_colors.set_ylabel('')
##clo_colors以列表的形式传入就不会显示标签了：[col_colors]
patches = [mpatches.Patch(color=[ '#ED0000', '#42B540','#00468B'][i], label="{:s}".format(['Metabolic','Basal','Mesenchymal'][i]) ) for i in range(len(color10_clu3))]
second = plt.legend(handles=patches, bbox_to_anchor = (-0.6, 2.2), loc='upper left', ncol=1, title='Group', frameon=False)
ax10 = plt.gca().add_artist(second)
second._legend_box.align = "left"


plt.savefig('heatmap_peptides-5C.png',dpi=300,bbox_inches='tight')
plt.savefig('heatmap_peptides-5C.pdf',bbox_inches='tight')
plt.close()

