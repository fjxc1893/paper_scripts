import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def deal1(f,key):
    df = pd.read_csv(f)
    #df1 = df[~df['Classification_level1'].str.contains('Human Diseases')]
    df1 = df.set_index('Term_description')
    df1['p-value'] = -np.log10(df1['p-value'])
    if 'del' in key:
        df1['p-value'] = -df1['p-value']
    return df1

def deal2(f,key):
    df = pd.read_csv(f)
    df['len'] = df['Term_description'].apply(lambda x:len(x))
    #df = df[df.len<41]
    df1 = df.set_index('Term_description')
    df1['p-value'] = -np.log10(df1['p-value'])
    if 'del' in key:
        df1['p-value'] = -df1['p-value']
    return df1

def add_pos(f,key):
    df = pd.read_csv(f,sep='\t')
    df['pos'] = df['Gene Symbol'].apply(lambda x:dic[x])
    df1 = df.set_index('Gene Symbol')[[key,'pos']]
    df1['pos'] = df1['pos']/1000000000
    if 'del' in key:
        df1[key] = -df1[key]
    return df1


cs = ['#000000','#FFFFFF']*12

kar = pd.read_csv('karyotype.txt',sep='\t',header=None).drop(23)
kar[0] = kar[0].apply(lambda x:x.split('r')[-1])
kar['cumsum'] = kar[5].cumsum()
k1 = kar.copy()
k1['cumsum'] = k1['cumsum']/1000000000
k1 = k1.iloc[:-1,:]
kar['cumsum'] = [0]+list(kar['cumsum'])[:-1]
kar['t_pos'] = kar['cumsum']+kar[5]*0.1
kar = kar[[0,5,'t_pos']].set_index(0)/1000000000
kar = kar.iloc[:-1,:]



kar1 = kar[[5]].T

ax8 = kar1.plot(kind='barh',stacked=True,edgecolor='k',color=sns.color_palette(cs),width=0.02)
plt.legend('')

dic = pd.read_csv('gene_pos.xls',sep='\t',index_col=0).set_index('gene')['gene_pos'].to_dict()

ag = pd.read_csv('all_gene.xls',sep='\t')

kag = add_pos('kinase_amp/kinase_amp.xls','kinase_amp')
kdg = add_pos('kinase_del/kinase_del.xls','kinase_del')
iag = add_pos('immune_amp/immune_amp.xls','immune_amp')
idg = add_pos('immune_del/immune_del.xls','immune_del')
eag = add_pos('ECM_amp/ECM_amp.xls','ECM_amp')
edg = add_pos('ECM_del/ECM_del.xls','ECM_del')



kake = deal1('kinase_amp/KEGG_Enrichment/kinase_amp_KEGG_enrichment_result.csv','amp')
kake = kake[kake.index.str.contains('MAPK|PI3K-Akt|Central carbon metabolism in cancer')][['p-value']]
kake['sig'] = 'pos'
kdke = deal1('kinase_del/KEGG_Enrichment/kinase_del_KEGG_enrichment_result.csv','del')
kdke = kdke[kdke.index.str.contains('MAPK|PI3K-Akt|Central carbon metabolism in cancer')][['p-value']]
kdke['sig'] = 'neg'
kd = pd.concat([kdke,kake])

iake = deal2('immune_amp/GO_Enrichment/immune_amp_GO_enrichment_result.csv','amp')
iake = iake[iake.index.str.contains('^actin cytoskeleton organization$|eIF3m|yotube differentiation')][['p-value']]
iake['sig'] = 'pos'
idke1 = deal2('immune_del/GO_Enrichment/immune_del_GO_enrichment_result.csv','del')
idke1 = idke1[idke1.index.str.contains('yotube differentiation|ell motility')][['p-value']]
idke2 = deal1('immune_del/KEGG_Enrichment/immune_del_KEGG_enrichment_result.csv','del')
idke2 = idke2[idke2.index.str.contains('Regulation of actin cytoskeleton')][['p-value']]
idke = pd.concat([idke1,idke2]).sort_values(by='p-value')
idke['sig'] = 'neg'
ids = pd.concat([idke,iake])

eake = deal1('ECM_amp/KEGG_Enrichment/ECM_amp_KEGG_enrichment_result.csv','amp')
eake = eake[eake.index.str.contains('Focal adhesion|Tight junction|ECM-receptor interaction')][['p-value']]
eake['sig'] = 'pos'
edke = deal1('ECM_del/KEGG_Enrichment/ECM_del_KEGG_enrichment_result.csv','del')
edke = edke[edke.index.str.contains('Focal adhesion|Tight junction|ECM-receptor interaction')][['p-value']]
edke['sig'] = 'neg'
ed = pd.concat([edke,eake])

#plt.bar(x=kag['pos'],height=kag['kinase_amp'],width=0.015,color='r',edgecolor='')
#plt.bar(x=kdg['pos'],height=kdg['kinase_del'],width=0.015,color='b',edgecolor='')

plt.rcParams['font.sans-serif'] = ['Arial']
sns.set_style('white')  
fig,((ax1,ax2),(ax3,ax4),(ax5,ax6),(ax7,ax8)) = plt.subplots(4,2,figsize=(10,6))
plt.subplots_adjust(wspace=0.1, hspace=0)
#fig,(ax1,ax2) = plt.subplots(1,2)
ax1.bar(x=kag['pos'],height=kag['kinase_amp'],width=0.015,color='r',edgecolor='')
ax1.bar(x=kdg['pos'],height=kdg['kinase_del'],width=0.015,color='b',edgecolor='')
# =============================================================================
# ll =  list(ax1.get_ylim())
# for k in k1['cumsum']:
#     ax1.plot([k,k],ll,linestyle="-",alpha=0.6,linewidth=0.5,color='k')
# =============================================================================
ax1.set_xticklabels('')
ax1.set_yticklabels([-0.4,0.4,0,0.4])
ax1.set_ylabel('Metabolic')

#,sharex = ax1
kd['p-value'].plot(kind='barh',width=0.8,edgecolor='',color=kd.sig.map({'pos': 'r', 'neg': 'b'}), ax=ax2)
#ax2.set_xlim(-5,6)
ax2.set_xlim(-6,6)
ax2.set_ylabel('')
ax2.set_xticklabels('')
ax2.set_yticklabels('')
ax2.text(-5.95,4.8,'Central carbon metabolism in cancer',color='r',size=8)
ax2.text(-4.1,3.8,'MAPK signaling pathway',color='r',size=8)
ax2.text(-4.5,2.8,'PI3K-Akt signaling pathway',color='r',size=8)
ax2.text(0.1,1.8,'MAPK signaling pathway',color='b',size=8)
ax2.text(0.1,0.8,'Central carbon metabolism in cancer',color='b',size=8)
ax2.text(0.1,-0.2,'PI3K-Akt signaling pathway',color='b',size=8)

#ax2.yaxis.tick_right()
#ax2.yaxis.label.set_fontsize(20)

ax3.bar(x=iag['pos'],height=iag['immune_amp'],width=0.015,color='r',edgecolor='')
ax3.bar(x=idg['pos'],height=idg['immune_del'],width=0.015,color='b',edgecolor='')
ax3.set_xticklabels('')
ax3.set_yticklabels([-0.4,0.25,0,0.25])
ax3.set_ylabel('Basal')
# =============================================================================
# ll =  list(ax3.get_ylim())
# for k in k1['cumsum']:
#  ax3.plot([k,k],ll,linestyle="-",alpha=0.6,linewidth=0.5,color='k')
# =============================================================================

ids['p-value'].plot(kind='barh',width=0.8,edgecolor='',color=ids.sig.map({'pos': 'r', 'neg': 'b'}), ax=ax4)
#ax4.set_xlim(-5,6)
ax4.set_xlim(-6,6)
ax4.set_ylabel('')
ax4.set_xticklabels('')
ax4.set_yticklabels('')
ax4.text(-5.8,4.8,r'Eukaryotic translation initiation factor 3M',color='r',size=7)
ax4.text(-3.75,3.8,'Myotube differentiation',color='r',size=8)
ax4.text(-5.05,2.8,'Actin cytoskeleton organization',color='r',size=8)
ax4.text(0.1,1.8,'Cell motility',color='b',size=8)
ax4.text(0.1,0.8,'Regulation of actin cytoskeleton',color='b',size=8)
ax4.text(0.1,-0.2,'Myotube differentiation',color='b',size=8)

ax5.bar(x=eag['pos'],height=eag['ECM_amp'],width=0.015,color='r',edgecolor='')
ax5.bar(x=edg['pos'],height=edg['ECM_del'],width=0.015,color='b',edgecolor='')
# =============================================================================
# ll =  list(ax5.get_ylim())
# for k in k1['cumsum']:
#  ax5.plot([k,k],ll,linestyle="-",alpha=0.6,linewidth=0.5,color='k')
# =============================================================================
ax5.set_xticklabels('')
ax5.set_yticklabels([-0.4,0.25,0,0.25])
ax5.set_xlabel('')
ax5.set_ylabel('Mesenchymal')



ed['p-value'].plot(kind='barh',width=0.8,edgecolor='',color=ed.sig.map({'pos': 'r', 'neg': 'b'}),ax=ax6)
ax6.set_xlim(-6,6)
ax6.set_xlabel('-log$_{10}$(p-value)')
ax6.set_ylabel('')
ax6.set_xticklabels([6,4,2,0,2,4,6])
ax6.set_yticklabels('')
ax6.text(-4.15,4.8,'ECM-receptor interaction',color='r',size=8)
ax6.text(-2.25,3.8,'Tight junction',color='r',size=8)
ax6.text(-2.55,2.8,'Focal adhesion',color='r',size=8)
ax6.text(0.1,1.8,'ECM-receptor interaction',color='b',size=8)
ax6.text(0.1,0.8,'Tight junction',color='b',size=8)
ax6.text(0.1,-0.2,'Focal adhesion',color='b',size=8)
ax8.remove()

l1 = pd.concat([kar.iloc[:17,:],kar.iloc[18,:].to_frame().T,kar.iloc[20,:].to_frame().T,kar.iloc[22,:].to_frame().T])
l2 = pd.concat([kar.iloc[17,:].to_frame().T,kar.iloc[19,:].to_frame().T,kar.iloc[21,:].to_frame().T])
kar1.plot(kind='barh',stacked=True,edgecolor='k',linewidth=1,color=sns.color_palette(cs),width=0.04,ax=ax7,align='edge')


ax7.sharex(ax1)
ax7.legend('')
ax7.set_ylabel('')
ax7.set_xlabel('')
ax7.set_yticklabels('')
ax7.set_xticklabels('')
ax7.axis('off')
for lab,x,y in zip(l1.index,l1.t_pos,[0.05]*20):
    ax7.text(x,y,lab,fontsize=8)
for lab,x,y in zip(l2.index,l2.t_pos,[-0.05]*3):
    ax7.text(x,y,lab,fontsize=8)
#ax7.set_ylim(0,2)



plt.savefig('barplot.png',dpi=300,bbox_inches='tight')
plt.savefig('barplot.pdf',bbox_inches='tight')
plt.close()
