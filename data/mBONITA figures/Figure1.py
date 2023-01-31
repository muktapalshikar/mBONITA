
# coding: utf-8

# # Rule inference on multi-omics networks - Data processing and Figure 1

# In[121]:


import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np
from upsetplot import from_contents, UpSet, plot
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
from math import ceil


# # Phosphoproteomics data analysis

# In[122]:


combined_clean = pd.read_excel(
    'raw_data/202007Quantitaive results_Cleaned.xlsx',
    sheet_name="Quan_Combined").fillna(0)
combined_clean.head()


# ## Log2 normalize data

# In[123]:


combined_clean.iloc[:, 4:51] = combined_clean.iloc[:, 4:51].apply(
    lambda x: np.log2(x + 1), raw=False)
combined_clean.iloc[0:4, 4:51]


# ## Make conditions matrix

# In[124]:


conditions = pd.read_excel('raw_data/202007Quantitaive results_Cleaned.xlsx',
                           sheet_name="condition_matrix").fillna(0)
conditions.head()


# ## Add gene names to expression matrix

# In[125]:


gene_mapping = pd.read_excel('raw_data/202007Quantitaive results_Cleaned.xlsx',
                             sheet_name="gene mapping").fillna(0)
gene_mapping.head()
combined_clean.loc[:, "Gene"] = [
    gene_mapping.To[gene_mapping.From == combined_clean.loc[
        i, "Protein AC"]].tolist()[0] if
    len(gene_mapping.To[gene_mapping.From == combined_clean.loc[i,
                                                                "Protein AC"]])
    > 0 else combined_clean.loc[i, "Protein AC"]
    for i in range(len(combined_clean))
]


# ## For rows corresponding to the same protein/gene - Retain only those rows with the highest abundance.

# In[126]:


combined_clean.select_dtypes(include=np.number).mean(axis=1)


# In[127]:


combined_clean.loc[:, "MeanAbundance"] = combined_clean.select_dtypes(
    include=np.number).mean(axis=1)
combined_clean.loc[:, "RowSums"] = combined_clean.select_dtypes(
    include=np.number).sum(axis=1)
#combined_clean = combined_clean.set_index("Gene", drop=False)
combined_clean = combined_clean.sort_values("MeanAbundance", ascending=False)
combined_clean = combined_clean.groupby("Gene", sort=False)
combined_clean = combined_clean.first()
combined_clean = combined_clean.reset_index()
combined_clean = combined_clean.drop(columns=["MeanAbundance", "RowSums"])
print(combined_clean.shape)


# # Process proteomics dataset

# In[128]:


proteomics = pd.read_csv("raw_data/ramos_genes.csv").fillna(0)
proteomics.loc[:, [
    "norm_19-A", "norm_19-B", "norm_19-C", "norm_1-A", "norm_1-B", "norm_1-C",
    "norm_1-1ug-A", "norm_1-1ug-B", "norm_1-1ug-C"
]] = proteomics.loc[:, [
    "norm_19-A", "norm_19-B", "norm_19-C", "norm_1-A", "norm_1-B", "norm_1-C",
    "norm_1-1ug-A", "norm_1-1ug-B", "norm_1-1ug-C"
]].apply(lambda x: np.log2(x + 1), raw=False)
proteomics.loc[:, [
    "norm_19-A", "norm_19-B", "norm_19-C", "norm_1-A", "norm_1-B", "norm_1-C",
    "norm_1-1ug-A", "norm_1-1ug-B", "norm_1-1ug-C"
]] = proteomics.loc[:, [
    "norm_19-A", "norm_19-B", "norm_19-C", "norm_1-A", "norm_1-B", "norm_1-C",
    "norm_1-1ug-A", "norm_1-1ug-B", "norm_1-1ug-C"
]].loc[~(proteomics.loc[:, [
    "norm_19-A", "norm_19-B", "norm_19-C", "norm_1-A", "norm_1-B", "norm_1-C",
    "norm_1-1ug-A", "norm_1-1ug-B", "norm_1-1ug-C"
]] == 0).all(axis=1)]
proteomics = proteomics.fillna(0)
prot_common_conditions = [
    "Gene", "norm_19-A", "norm_19-B", "norm_19-C", "norm_1-A", "norm_1-B",
    "norm_1-C", "norm_1-1ug-A", "norm_1-1ug-B", "norm_1-1ug-C"
]


# # Process transcriptomics dataset

# In[129]:


transcriptomics = pd.read_csv("raw_data/unfiltered_rpm_counts.txt",
                              sep="\t").fillna(0)
transcriptomics.loc[:, [
    "Ramos_19O2_NoCyclo_1", "Ramos_19O2_NoCyclo_2", "Ramos_19O2_NoCyclo_3",
    "Ramos_19O2_PlusCyclo_1", "Ramos_19O2_PlusCyclo_2",
    "Ramos_19O2_PlusCyclo_3", "Ramos_1O2_NoCyclo_1", "Ramos_1O2_NoCyclo_2",
    "Ramos_1O2_NoCyclo_3", "Ramos_1O2_PlusCyclo_1", "Ramos_1O2_PlusCyclo_2",
    "Ramos_1O2_PlusCyclo_3"
]] = transcriptomics.loc[:, [
    "Ramos_19O2_NoCyclo_1", "Ramos_19O2_NoCyclo_2", "Ramos_19O2_NoCyclo_3",
    "Ramos_19O2_PlusCyclo_1", "Ramos_19O2_PlusCyclo_2",
    "Ramos_19O2_PlusCyclo_3", "Ramos_1O2_NoCyclo_1", "Ramos_1O2_NoCyclo_2",
    "Ramos_1O2_NoCyclo_3", "Ramos_1O2_PlusCyclo_1", "Ramos_1O2_PlusCyclo_2",
    "Ramos_1O2_PlusCyclo_3"
]].apply(lambda x: np.log2(x + 1), raw=False)
transcriptomics.loc[:, [
    "Ramos_19O2_NoCyclo_1", "Ramos_19O2_NoCyclo_2", "Ramos_19O2_NoCyclo_3",
    "Ramos_19O2_PlusCyclo_1", "Ramos_19O2_PlusCyclo_2",
    "Ramos_19O2_PlusCyclo_3", "Ramos_1O2_NoCyclo_1", "Ramos_1O2_NoCyclo_2",
    "Ramos_1O2_NoCyclo_3", "Ramos_1O2_PlusCyclo_1", "Ramos_1O2_PlusCyclo_2",
    "Ramos_1O2_PlusCyclo_3"
]] = transcriptomics.loc[:, [
    "Ramos_19O2_NoCyclo_1", "Ramos_19O2_NoCyclo_2", "Ramos_19O2_NoCyclo_3",
    "Ramos_19O2_PlusCyclo_1", "Ramos_19O2_PlusCyclo_2",
    "Ramos_19O2_PlusCyclo_3", "Ramos_1O2_NoCyclo_1", "Ramos_1O2_NoCyclo_2",
    "Ramos_1O2_NoCyclo_3", "Ramos_1O2_PlusCyclo_1", "Ramos_1O2_PlusCyclo_2",
    "Ramos_1O2_PlusCyclo_3"
]].loc[~(transcriptomics.loc[:, [
    "Ramos_19O2_NoCyclo_1", "Ramos_19O2_NoCyclo_2", "Ramos_19O2_NoCyclo_3",
    "Ramos_19O2_PlusCyclo_1", "Ramos_19O2_PlusCyclo_2",
    "Ramos_19O2_PlusCyclo_3", "Ramos_1O2_NoCyclo_1", "Ramos_1O2_NoCyclo_2",
    "Ramos_1O2_NoCyclo_3", "Ramos_1O2_PlusCyclo_1", "Ramos_1O2_PlusCyclo_2",
    "Ramos_1O2_PlusCyclo_3"
]] == 0).all(axis=1)]
transcriptomics = transcriptomics.fillna(0)


# # Find common genes between all three datasets - all conditions, all measured genes included

# In[130]:


commonGenes = from_contents({
    'Proteomics': set(proteomics.index),
    'Transcriptomics': set(transcriptomics.index),
    'Phosphoproteomics': set(combined_clean.index)
})
commonGenes


# In[131]:


with plt.style.context('tableau-colorblind10'):
    plt.figure(figsize=(2, 5))
    plot(commonGenes,
         show_counts=True,
         min_subset_size=1,
         facecolor="midnightblue",
         orientation='horizontal')
    plt.savefig("figure1a_unfiltered.pdf")
    plt.savefig("figure1a_unfiltered.png", dpi=600)


# # Make datasets that just have the common genes and common conditions/samples

# ## Common conditions

# ### Phosphoproteomics

# In[132]:


group1 = conditions['TMT.group'] == 1
lowO2 = conditions.O2 == "Low"
noCXCL12 = conditions.CXCL12 == "No"
group2 = conditions['TMT.group'] == 2
highO2 = conditions.O2 == "High"
noCyA = conditions.CyA == "No"
phosph_common_conditions = conditions.loc[
    (group2 & highO2 & noCXCL12 & noCyA) | (group1 & lowO2 & noCXCL12),
    'Column labels - cleaned dataset'].tolist()[
        0:9]  #ignore last three samples - repeat
phosph_common_conditions.insert(0, "Gene")
print(phosph_common_conditions)


# ### Proteomics

# In[133]:


prot_common_conditions = [
    "Gene", "norm_19-A", "norm_19-B", "norm_19-C", "norm_1-A", "norm_1-B",
    "norm_1-C", "norm_1-1ug-A", "norm_1-1ug-B", "norm_1-1ug-C"
]


# ### Transcriptomics

# In[134]:


transcript_common_conditions = [
    'Gene', "Ramos_19O2_NoCyclo_1", "Ramos_19O2_NoCyclo_2",
    "Ramos_19O2_NoCyclo_3", "Ramos_1O2_NoCyclo_1", "Ramos_1O2_NoCyclo_2",
    "Ramos_1O2_NoCyclo_3", "Ramos_1O2_PlusCyclo_1", "Ramos_1O2_PlusCyclo_2",
    "Ramos_1O2_PlusCyclo_3"
]


# ## Create datasets

# In[135]:


# transcriptomics
transcript_common = transcriptomics[set(
    transcript_common_conditions).intersection(set(transcriptomics.columns))]
transcript_common.set_index("Gene", inplace=True)
#transcript_common.sort_index(inplace=True)
transcript_common = transcript_common[transcript_common.median(axis=1) > 0]

# phosphoproteomics
phospho_common = combined_clean[set(phosph_common_conditions).intersection(
    set(combined_clean))]
phospho_common.set_index("Gene", inplace=True)
#phospho_common.sort_index(inplace=True)
phospho_common = phospho_common[phospho_common.median(axis=1) > 0]

# proteomics
prot_common = proteomics[set(prot_common_conditions).intersection(
    set(proteomics.columns))]
prot_common.set_index("Gene", inplace=True)
#prot_common.sort_index(inplace=True)
prot_common = prot_common[prot_common.median(axis=1) > 0]


# In[136]:


phospho_common.to_csv("bonita_phosphoproteomics.csv")
prot_common.to_csv("bonita_proteomics.csv")
transcript_common.to_csv("bonita_transcriptomics.csv")
print(phospho_common.shape, transcript_common.shape, prot_common.shape)


# # Make figures showing overlap in genes, using filtered datasets

# In[137]:


commonGenes = from_contents({
    'Proteomics': set(prot_common.index),
    'Transcriptomics': set(transcript_common.index),
    'Phosphoproteomics': set(phospho_common.index)
})
commonGenes
with plt.style.context('tableau-colorblind10'):
    plt.figure(figsize=(2, 5))
    plot(commonGenes,
         show_counts=True,
         min_subset_size=1,
         facecolor="green",
         orientation='horizontal')
    plt.savefig("figure1a_filtered.pdf")
    plt.savefig("figure1a_filtered.png", dpi=600)


# In[138]:


commonGene = list(
    set(
        set(prot_common.index) & set(transcript_common.index)
        & set(phospho_common.index)))

transcript_common = transcript_common.loc[transcript_common.index.isin(
    commonGene)]
phospho_common = phospho_common.loc[phospho_common.index.isin(commonGene)]
prot_common = prot_common.loc[prot_common.index.isin(commonGene)]

print(transcript_common.shape, phospho_common.shape, prot_common.shape)


# # Calculate pairwise correlations between datasets

# In[139]:


temp = transcript_common.reset_index(drop=False).melt(id_vars='Gene')
temp['Condition'] = [i[:-2] for i in temp.variable]
for i in range(len(temp['Condition'])):
    if temp.loc[i, 'Condition'] == 'Ramos_19O2_NoCyclo':
        temp.loc[i, 'Condition'] = "19% O2, CyA-"
    elif temp.loc[i, 'Condition'] == 'Ramos_1O2_NoCyclo':
        temp.loc[i, 'Condition'] = "1% O2, CyA-"
    elif temp.loc[i, 'Condition'] == 'Ramos_1O2_PlusCyclo':
        temp.loc[i, 'Condition'] = "1% O2, CyA+"
a = pd.DataFrame(temp.groupby(['Condition', 'Gene']).median())  #.reset_index()
a


# In[140]:


temp2 = temp[['variable', 'Condition']].drop_duplicates().reset_index(drop=True)
temp2['O2'] = [i[:-6] for i in temp2.Condition]
temp2['CyA'] = [i[-4:] for i in temp2.Condition]
for i in set(temp2.Condition):
    temp2[i] = np.nan
    for j in range(len(temp2)):
        if (temp2.loc[j, 'Condition'] == i):
            temp2.loc[j, i] = 1
        else:
            temp2.loc[j, i] = 0
temp2[["1% O2, CyA+","1% O2, CyA-","19% O2, CyA-"]].to_csv("transcriptomics_conditions.csv")
temp2


# In[141]:


temp = phospho_common.reset_index(drop=False).melt(id_vars='Gene')
temp.columns = ['Gene', 'Condition', 'value']
for i in range(len(temp['Condition'])):
    oldcond = temp.loc[i, 'Condition']
    newcond = ''
    if conditions.loc[conditions['Column labels - cleaned dataset'] == oldcond,
                      "O2"].tolist()[0] == "Low":
        newcond = newcond + "1% O2, "
    else:
        newcond = newcond + "19% O2, "
    if conditions.loc[conditions['Column labels - cleaned dataset'] == oldcond,
                      "CyA"].tolist()[0] == "No":
        newcond = newcond + "CyA-"
    else:
        newcond = newcond + "CyA+"
    temp.loc[i, 'Condition'] = newcond
    temp.loc[i, 'variable'] = oldcond
b = pd.DataFrame(temp.groupby(['Condition', 'Gene']).median())  #.reset_index()
b


# In[142]:


b


# In[143]:


temp2 = temp[['variable', 'Condition']].drop_duplicates().reset_index(drop=True)
temp2['O2'] = [i[:-6] for i in temp2.Condition]
temp2['CyA'] = [i[-4:] for i in temp2.Condition]
for i in set(temp2.Condition):
    temp2[i] = np.nan
    for j in range(len(temp2)):
        if (temp2.loc[j, 'Condition'] == i):
            temp2.loc[j, i] = 1
        else:
            temp2.loc[j, i] = 0
temp2[["1% O2, CyA+","1% O2, CyA-","19% O2, CyA-"]].to_csv("phosphoproteomics_conditions.csv")
temp2


# In[144]:


temp = prot_common.reset_index(drop=False).melt(id_vars='Gene')
temp['Condition'] = [i[:-2] for i in temp.variable]

for i in range(len(temp['Condition'])):
    if temp.loc[i, 'Condition'] == 'norm_19':
        temp.loc[i, 'Condition'] = "19% O2, CyA-"
    elif temp.loc[i, 'Condition'] == 'norm_1':
        temp.loc[i, 'Condition'] = "1% O2, CyA-"
    elif temp.loc[i, 'Condition'] == 'norm_1-1ug':
        temp.loc[i, 'Condition'] = "1% O2, CyA+"

c = pd.DataFrame(temp.groupby(['Condition', 'Gene']).median())
c


# In[153]:


temp2 = temp[['variable', 'Condition']].drop_duplicates().reset_index(drop=True)
temp2['O2'] = [i[:-6] for i in temp2.Condition]
temp2['CyA'] = [i[-4:] for i in temp2.Condition]
for i in set(temp2.Condition):
    temp2[i] = np.nan
    for j in range(len(temp2)):
        if (temp2.loc[j, 'Condition'] == i):
            temp2.loc[j, i] = 1
        else:
            temp2.loc[j, i] = 0
temp2[["1% O2, CyA+","1% O2, CyA-","19% O2, CyA-"]].to_csv("proteomics_conditions.csv", index=False)
temp2


# In[146]:


medianDF = pd.DataFrame({
    "Transcriptomics": a.value,
    "Proteomics": c.value,
    "Phosphoproteomics": b.value
}).reset_index()
medianDF


# In[147]:


g = sns.pairplot(medianDF,
                 diag_kind='hist',
                 kind='reg',
                 plot_kws={'scatter_kws': {
                     'alpha': 0.25
                 }},
                 markers=["o", "s", "D"],
                 corner=False,
                 hue='Condition',
                 palette='colorblind')
g.savefig("figure1b.pdf")
g.savefig("figure1b.png", dpi=600)


# In[148]:


a = prot_common.T.reset_index(drop=True).corrwith(
    transcript_common.T.reset_index(drop=True), axis=0)
sns.histplot(a, color="blue")
b = phospho_common.T.reset_index(drop=True).corrwith(
    transcript_common.T.reset_index(drop=True), axis=0)
sns.histplot(b, color="red")
c = phospho_common.T.reset_index(drop=True).corrwith(
    prot_common.T.reset_index(drop=True), axis=0)
sns.histplot(c, color="brown")


# In[149]:


a = transcript_common.reset_index(drop=False).melt(
    id_vars="Gene").sort_values("Gene")
b = prot_common.reset_index(drop=False).melt(
    id_vars="Gene").sort_values("Gene")
c = phospho_common.reset_index(drop=False).melt(
    id_vars="Gene").sort_values("Gene")


# In[150]:


from scipy.stats import spearmanr, pearsonr
print(spearmanr(a.value, b.value), pearsonr(a.value, b.value))


# In[151]:


print(spearmanr(a.value, c.value), pearsonr(a.value, c.value))


# In[152]:


print(spearmanr(b.value, c.value), pearsonr(b.value, c.value))

"""
# # Make proteomics coexpression network

# In[ ]:


proteomics_coexp = proteomics[set(
    prot_common_conditions).intersection(set(proteomics.columns))]
proteomics_coexp.set_index("Gene", inplace=True)
proteomics_coexp = proteomics_coexp.T
proteomics_coexp = proteomics_coexp.corr(method="spearman")
proteomics_coexp.head
proteomics_coexp.to_csv("proteomics_corr.csv")


# In[ ]:


proteomics_coexp = pd.read_csv("proteomics_corr.csv", index_col=0)
proteomics_coexp['Gene'] = list(proteomics_coexp.index)
proteomics_coexp = proteomics_coexp.melt(id_vars=['Gene']).dropna()
proteomics_coexp.value = proteomics_coexp.value.abs()
proteomics_coexp = proteomics_coexp[proteomics_coexp.value >= 0.75]
proteomics_coexp = proteomics_coexp[proteomics_coexp.value != 1]
#proteomics_coexp.head
proteomics_net = nx.from_pandas_edgelist(proteomics_coexp,
                                         source="Gene",
                                         target="variable",
                                         edge_attr="value")
nx.set_node_attributes(proteomics_net,
                       nx.betweenness_centrality(proteomics_net),
                       "prot_betweenness")
#Gcc = sorted(nx.connected_components(proteomics_net), key=len, reverse=True)
#proteomics_net = proteomics_net.subgraph(Gcc[0])
nx.write_graphml_lxml(proteomics_net, "proteomics_net.graphml")
len(proteomics_net)
"""
proteomics_net = nx.read_graphml("proteomics_net.graphml")
len(proteomics_net)


# # Make transcriptomics coexpression network

# In[ ]:

"""
transcriptomics_coexp = transcriptomics[set(
    transcript_common_conditions).intersection(set(transcriptomics.columns))]
transcriptomics_coexp = transcriptomics_coexp.T
transcriptomics_coexp = transcriptomics_coexp.corr(method="spearman")
transcriptomics_coexp.head
transcriptomics_coexp.to_csv("transcriptomics_corr.csv")


# In[67]:


transcriptomics_coexp = pd.read_csv("transcriptomics_corr.csv", index_col=0)
transcriptomics_coexp['Gene'] = list(transcriptomics_coexp.index)
transcriptomics_coexp = transcriptomics_coexp.melt(id_vars=['Gene']).dropna()
transcriptomics_coexp.value = transcriptomics_coexp.value.abs()
transcriptomics_coexp = transcriptomics_coexp[
    transcriptomics_coexp.value >= 0.75]
transcriptomics_coexp = transcriptomics_coexp[transcriptomics_coexp.value != 1]


transcriptomics_net = nx.from_pandas_edgelist(transcriptomics_coexp,
                                              source="Gene",
                                              target="variable",
                                              edge_attr="value")
nx.set_node_attributes(transcriptomics_net,
                       nx.betweenness_centrality(transcriptomics_net),
                       "trans_betweenness")
#Gcc = sorted(nx.connected_components(transcriptomics_net), key=len, reverse=True)
#transcriptomics_net = transcriptomics_net.subgraph(Gcc[0])
nx.write_graphml_lxml(transcriptomics_net, "transcriptomics_net.graphml")
len(transcriptomics_net)
"""

transcriptomics_net = nx.read_graphml("transcriptomics_net.graphml")

# # Make phosphoproteomics coexpression network

# In[ ]:


phospho_coexp = combined_clean[set(
    phosph_common_conditions).intersection(set(combined_clean.columns))]
phospho_coexp = phospho_coexp.T
phospho_coexp = phospho_coexp.corr(method="spearman")
phospho_coexp.head
phospho_coexp.to_csv("phospho_corr.csv")


# In[ ]:


phospho_coexp = pd.read_csv("phospho_corr.csv", index_col=0)
phospho_coexp['Gene'] = list(phospho_coexp.index)
phospho_coexp = phospho_coexp.melt(id_vars=['Gene']).dropna()
phospho_coexp.value = phospho_coexp.value.abs()
phospho_coexp = phospho_coexp[phospho_coexp.value >= 0.75]
phospho_coexp = phospho_coexp[phospho_coexp.value != 1]

phospho_net = nx.from_pandas_edgelist(phospho_coexp,
                                      source="Gene",
                                      target="variable",
                                      edge_attr="value")
nx.set_node_attributes(phospho_net, nx.betweenness_centrality(phospho_net),
                       "phospho_betweenness")
#Gcc = sorted(nx.connected_components(phospho_net), key=len, reverse=True)
#phospho_net = phospho_net.subgraph(Gcc[0])
nx.write_graphml_lxml(phospho_net, "phospho_net.graphml")


# # Find edges overlapping between coexpression networks

# In[ ]:


def getOverlapGraph(netList=[], returnGiantComponent=True):
    overlapNet = nx.intersection_all(netList)
    if returnGiantComponent:
        Gcc = sorted(nx.connected_components(overlapNet), key=len, reverse=True)
        G0 = overlapNet.subgraph(Gcc[0])
        return G0
    else:
        return overlapNet


# In[ ]:


# transcriptomics - proteomics 
trans_prot = getOverlapGraph(netList=[transcriptomics_net, proteomics_net])
print(len(trans_prot))


# In[ ]:


# transcriptomics - phosphoproteomics 
trans_phospho = getOverlapGraph(netList=[transcriptomics_net, phospho_net])
print(len(trans_phospho))


# In[ ]:


# phosphoproteomics - proteomics
prot_phospho = getOverlapGraph(netList=[proteomics_net, phospho_net])
print(len(prot_phospho))


# In[ ]:


# make consensus network

consensus_net = getOverlapGraph(netList=[transcriptomics_net, proteomics_net, phospho_net], returnGiantComponent=False)

nx.set_node_attributes(consensus_net, nx.betweenness_centrality(consensus_net), "consensus_betweenness")

nx.write_graphml_lxml(consensus_net, "consensus_net.graphml")

consensus_net_giant = getOverlapGraph(netList=[transcriptomics_net, proteomics_net, phospho_net], returnGiantComponent=True)

nx.write_graphml_lxml(consensus_net_giant, "consensus_net_largest_connected_component.graphml")

print(len(consensus_net_giant))


# In[ ]:


# make upset plot showing intersections between coexpression network edges

commonNodes = from_contents({
    'Proteomics': set(proteomics_net.edges),
    'Transcriptomics': set(transcriptomics_net.edges),
    'Phosphoproteomics': set(phospho_net.edges)
})
with plt.style.context('tableau-colorblind10'):
    plt.figure(figsize=(2, 5))
    plot(commonNodes,
         show_counts=True,
         min_subset_size=1,
         facecolor="navy",
         orientation='horizontal')
    plt.savefig("figure2c.pdf")
    plt.savefig("figure2c.png", dpi=600)


# In[ ]:


def enrichmentNodes(networkNodes):
    consensusNodesEnrich = gp.enrichr(gene_list=list(networkNodes),
                                      gene_sets=['KEGG_2021_Human'],
                                      organism='Human',
                                      cutoff=0.05)
    consensusNodesEnrich.results.Genes = [
        temp.split(';') for temp in consensusNodesEnrich.results.Genes.tolist()
    ]
    enrichr_common_nodes = consensusNodesEnrich.results[
        consensusNodesEnrich.results['Adjusted P-value'] < 0.01]
    enrichr_common_nodes = enrichr_common_nodes.assign(
        log10_adjusted_p_value=[(-1) * np.log10(i)
                                for i in enrichr_common_nodes['Adjusted P-value']])
    return enrichr_common_nodes

def makeEnrichBubblePlot(enrichResult):
    with plt.style.context('tableau-colorblind10'):
        sns.set_context(
        "paper",
        rc={
        "font.size": 16,
        "axes.labelsize": 'medium',
        'ytick.labelsize': 'medium',
        'xtick.labelsize': 'medium',
        'axes.titlesize': 'medium',
        'legend.fontsize': 'medium',
        })
        sns.scatterplot(data=enrichResult,
                        y="Term",
                        x="log10_adjusted_p_value",
                        palette="Blues", s = 75)
        plt.xlabel("-log10 (adjusted p-value)")
        plt.ylabel("")
        axes = plt.gca()
        axes.yaxis.grid(color='grey',
                    linestyle=(0, (5, 10)),
                    linewidth=0.5)
        axes.xaxis.grid(color='grey',
                    linestyle=(0, (5, 10)),
                    linewidth=0.5)
        plt.xticks(range(0,ceil(max(enrichResult["log10_adjusted_p_value"])+1)))


# In[ ]:


# consensus network
consensusNodesEnrich = enrichmentNodes(consensus_net_giant.nodes)
plt.figure(figsize=(3,4))
makeEnrichBubblePlot(consensusNodesEnrich)


# In[ ]:


# proteomics only
len(set(transcriptomics_net.nodes))# - set(transcriptomics_net.nodes)) #.difference(set(phospho_net.nodes)))
#proteomics_net.subgraph(proteomics_net.nodes.in)
#proteomicsNodesEnrich = enrichmentNodes(G0.nodes)
#plt.figure(figsize=(5,10))
#makeEnrichBubblePlot(proteomicsNodesEnrich)


# In[ ]:


# transcriptomics only
overlapNet = nx.difference(transcriptomics_net, consensus_net.subgraph(transcriptomics_net.nodes))
Gcc = sorted(nx.connected_components(overlapNet), key=len, reverse=True)
G0 = overlapNet.subgraph(Gcc[0])
transcriptomicsNodesEnrich = enrichmentNodes(G0.nodes)
plt.figure(figsize=(5,10))
makeEnrichBubblePlot(transcriptomicsNodesEnrich)


# In[ ]:


# phosphoproteomics only
phosphoproteomicsNodesEnrich = enrichmentNodes(set(phospho_net.nodes) - set(consensus_net_giant.nodes))
plt.figure(figsize=(3,10))
makeEnrichBubblePlot(phosphoproteomicsNodesEnrich)


# In[ ]:


len(set(phospho_net.nodes) - set(consensus_net_giant.nodes))
len(set(proteomics_net.nodes) - set(consensus_net_giant.nodes))

