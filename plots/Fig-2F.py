import pandas as pd
from scipy.stats import zscore
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib as mpl
import matplotlib.cm as cm

##################读取处理数据#####################
BRAF_score = pd.read_excel("附件二、Scores.xlsx").dropna()
#BRAF_score由高到低排序
BRAF_score = BRAF_score.sort_values(by="BRAF_score",ascending=False)

protein_abundance = pd.read_csv("25个蛋白表达矩阵.xls", sep='\t')
del protein_abundance["Accession"]
protein_abundance.index = protein_abundance["Gene Name"]
protein_abundance.index.name = None
del protein_abundance["Gene Name"]
#按BRAD_score分数由高到低重排样本顺序
protein_abundance = protein_abundance[BRAF_score["Sample ID"]]
protein_abundance_zscore = pd.DataFrame(zscore(protein_abundance,axis=1),index=protein_abundance.index,columns=protein_abundance.columns)
#######################绘图########################
fig = plt.figure(figsize=(12, 6))
ax = fig.add_axes([0.1,0.1,0.8,0.8])

rows_num,cols_num = protein_abundance_zscore.shape
for i in range(rows_num):
    protein_cmap = LinearSegmentedColormap.from_list('chaos', ["#00BFFF","white","#B0080F"])
    protein_norm = mpl.colors.Normalize(vmin=min(protein_abundance_zscore.iloc[i,:]), vmax=max(protein_abundance_zscore.iloc[i,:]))
    protein_mapper = cm.ScalarMappable(norm=protein_norm, cmap=protein_cmap)
    for j in range(cols_num):
        ax.add_patch(plt.Rectangle(xy=(j/102, (27-i-3)/27), width=1/102, height=1/27,facecolor=protein_mapper.to_rgba(protein_abundance_zscore.iloc[i,j]), edgecolor='white',linewidth=0.1))

BRAF_score = BRAF_score["BRAF_score"].tolist()
BRAF_cmap = LinearSegmentedColormap.from_list('chaos', ["#268D21","yellow","#FF4500"])
BRAF_norm = mpl.colors.Normalize(vmin=min(BRAF_score),vmax=max(BRAF_score))
BRAF_mapper = cm.ScalarMappable(norm=BRAF_norm, cmap=BRAF_cmap)
for j in range(cols_num):
    ax.add_patch(plt.Rectangle(xy=(j/102, 25/27), width=1/102, height=1/27,facecolor=BRAF_mapper.to_rgba(BRAF_score[j]), edgecolor='white',linewidth=0.1))
    if protein_abundance_zscore.columns[j] == "FUSCC-41":
        ax.add_patch(plt.Rectangle(xy=(j / 102, 26 / 27), width=1 / 102, height=1 / 27,facecolor="grey", edgecolor='white', linewidth=0.1))
        plt.text(j / 102-0.05,1.01,"(BRAF T599dup)")
        plt.text(j / 102-0.03,1.04,"FUSCC-41")
    elif protein_abundance_zscore.columns[j] == "FUSCC-14":
        ax.add_patch(plt.Rectangle(xy=(j / 102, 26 / 27), width=1 / 102, height=1 / 27,facecolor="grey", edgecolor='white', linewidth=0.1))
        plt.text(j / 102-0.045,1.01,"(BRAF G469A)")
        plt.text(j / 102-0.03,1.04,"FUSCC-14")
    else:
        ax.add_patch(plt.Rectangle(xy=(j / 102, 26 / 27), width=1 / 102, height=1 / 27, facecolor="black",edgecolor='white', linewidth=0.1))

plt.yticks([1/27*i+1/54 for i in range(27)],list(protein_abundance_zscore.index)+["BRAF Score","BRAF mutation"])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
plt.xticks([])
ax.tick_params(bottom=False,top=False,left=False,right=False)


plt.savefig("BRAF.pdf")