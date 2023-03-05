
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import statsmodels.stats.multitest
import os
import re
import requests
from ast import literal_eval

def unique(sequence):
    seen = set()
    return [x for x in sequence if not (x in seen or seen.add(x))]

def getPathwayName(hsaURL):
    """Use KEGG API to get the readable name using the KEGG code (eg. hsa00010 corresponds to the glycolysis pathway)"""
    fileReg = re.compile("NAME\s+(\w+.*)")
    pathwayFile = requests.get("http://rest.kegg.jp/get/" + hsaURL,
                               stream=True)
    for line in pathwayFile.iter_lines():
        line = line.decode("utf-8")
        result = fileReg.match(line)
        if result:
            return result.group(1)
    return hsaURL

def processBonitaPvalues(pvalues,
                         adjustPvalue=True,
                         method=[
                             "bonferroni", "sidak", "holm_sidak", "sidak",
                             "simes_hochberg", "hommel", "fdr_bh", "fdr_by",
                             "fdr_tsbh", "fdr_tsbky"
                         ],
                         alpha=0.05, showOnlyTopTen=True):
    pvalues = pd.read_csv(
        pvalues, index_col=0).reset_index(drop=True).set_index("pathway")
    cols = pvalues.columns.difference(["pathway","code","nodes"])
    if adjustPvalue:
        for i in cols:
            pvalues[i] = np.power(10, -1 * pvalues[i])
            pvalues[i] = statsmodels.stats.multitest.multipletests(
                pvalues[i], method=method, is_sorted=False,
                returnsorted=False)[1]
            pvalues[i] = [(-1) * np.log10(j) for j in pvalues[i]]

    melt_pvalues = pvalues.reset_index().melt(id_vars='pathway', value_vars=cols)
    melt_pvalues = melt_pvalues[melt_pvalues.value > ((-1) * np.log10(alpha))]
    melt_pvalues = melt_pvalues.sort_values('value', ascending=False)
    melt_pvalues.columns = ("Pathway Name", "Condition", "value")
    if showOnlyTopTen:
        #pvalues = pvalues.loc[pvalues[cols].sum(axis=1).sort_values(ascending=False)[0:10].index]
        print(len(unique(melt_pvalues['Pathway Name'])[0:min(10, len(melt_pvalues))]))
        melt_pvalues = melt_pvalues[melt_pvalues['Pathway Name'].isin(unique(melt_pvalues['Pathway Name'])[0:min(10, len(melt_pvalues))])]
        print(unique(melt_pvalues['Pathway Name'])[0:min(10, len(melt_pvalues))])
        return melt_pvalues
    return melt_pvalues

method='bonferroni'
alpha=1

proteomics_pvalues = processBonitaPvalues(
    pvalues="pvalues_proteomics_20220601.csv",
    adjustPvalue=True,
    method=method,
    alpha=alpha)

sns.set_theme(context='paper', style='ticks')
g = sns.scatterplot(data=proteomics_pvalues,
                    x='value',
                    y="Pathway Name",
                    hue="Condition",
                    palette="Greens",
                    s=75)
plt.legend(
    loc='upper center',  #'center left',  #
    bbox_to_anchor=(0.15, -0.15),  # (1.1, 0.6),  #
    fancybox=True,
    shadow=False,
    ncol=1,
    borderaxespad=0,
    facecolor="white")
plt.ylabel("")
plt.xlabel("-log10 (adjusted p-value)")
axes = plt.gca()
axes.yaxis.grid(color='grey', linestyle=(0, (5, 10)), linewidth=0.5)
g.figure.set_figheight(8)
g.figure.set_figwidth(5)
g.figure.tight_layout()
g.figure.savefig("proteomics_bonita.pdf")  #, height=10, width=4)
g.figure.savefig("proteomics_bonita.png", dpi=300)  #, height=10, width=4)
plt.clf()
transcriptomics_pvalues = processBonitaPvalues(
    pvalues="pvalues_transcriptomics_20220601.csv",
    adjustPvalue=True,
    method=method,
    alpha=alpha)
sns.set_theme(context='paper', style='ticks')
g = sns.scatterplot(data=transcriptomics_pvalues,
                    x='value',
                    y="Pathway Name",
                    hue="Condition",
                    palette="Greens",
                    s=75)
plt.legend(
    loc='upper center',  #'center left',  #
    bbox_to_anchor=(0.15, -0.15),  # (1.1, 0.6),  #
    fancybox=True,
    shadow=False,
    ncol=1,
    borderaxespad=0,
    facecolor="white")
plt.ylabel("")
plt.xlabel("-log10 (p-value)")
axes = plt.gca()
axes.yaxis.grid(color='grey', linestyle=(0, (5, 10)), linewidth=0.5)
g.figure.set_figheight(8)
g.figure.set_figwidth(5)
g.figure.tight_layout()
g.figure.savefig("transcriptomics_bonita.pdf")  #, height=10, width=4)
g.figure.savefig("transcriptomics_bonita.png", dpi=300)  #, height=10, width=4)
plt.clf()
phosphoproteomics_pvalues = processBonitaPvalues(
    pvalues="pvalues_phosphoproteomics_20220601.csv",
    adjustPvalue=True,
    method=method,
    alpha=alpha)
sns.set_theme(context='paper', style='ticks')
g = sns.scatterplot(data=phosphoproteomics_pvalues,
                    x='value',
                    y="Pathway Name",
                    hue="Condition",
                    palette="Greens",
                    s=75)
plt.legend(
    loc='upper center',  #'center left',  #
    bbox_to_anchor=(0.15, -0.15),  # (1.1, 0.6),  #
    fancybox=True,
    shadow=False,
    ncol=1,
    borderaxespad=0,
    facecolor="white")
plt.ylabel("")
plt.xlabel("-log10 (p-value)")
axes = plt.gca()
axes.yaxis.grid(color='grey', linestyle=(0, (5, 10)), linewidth=0.5)
g.figure.set_figheight(8)
g.figure.set_figwidth(5)
g.figure.tight_layout()
g.figure.savefig("phosphoproteomics_bonita.pdf")
g.figure.savefig("phosphoproteomics_bonita.png", dpi=300)
plt.clf()
concatenated_pvalues = processBonitaPvalues(
    pvalues="pvalues_concatenated_20220601.csv",
    adjustPvalue=True,
    method=method,
    alpha=alpha)
sns.set_theme(context='paper', style='ticks')
g = sns.scatterplot(data=concatenated_pvalues,
                    x='value',
                    y="Pathway Name",
                    hue="Condition",
                    palette="Greens",
                    s=75)
plt.legend(
    loc='upper center',  #'center left',  #
    bbox_to_anchor=(0.15, -0.15),  # (1.1, 0.6),  #
    fancybox=True,
    shadow=False,
    ncol=1,
    borderaxespad=0,
    facecolor="white")
plt.ylabel("")
plt.xlabel("-log10 (p-value)")
axes = plt.gca()
axes.yaxis.grid(color='grey', linestyle=(0, (5, 10)), linewidth=0.5)
g.figure.set_figheight(8)
g.figure.set_figwidth(5)
g.figure.tight_layout()
g.figure.savefig("concatenated_bonita.pdf")
g.figure.savefig("concatenated_bonita.png", dpi=300)
plt.clf()

transcriptomics_pvalues["Dataset"] = "Transcriptomics"
phosphoproteomics_pvalues["Dataset"] = "Phosphoproteomics"
proteomics_pvalues["Dataset"] = "Proteomics"
all_pvalues = pd.concat(
    [transcriptomics_pvalues, proteomics_pvalues, phosphoproteomics_pvalues])
print(all_pvalues.shape)
all_pvalues.columns = ["Pathway Name", "Condition", "-log10 (adjusted p-value)", "Dataset"]
sns.set_theme(context='talk', style='ticks')
g = sns.relplot(
    data=all_pvalues,
    x='-log10 (adjusted p-value)',
    y="Pathway Name",
    hue="Condition",
    col ="Dataset",
    col_order = ["Proteomics", "Phosphoproteomics", "Transcriptomics"], 
    style = "Condition",
    palette="colorblind",
    s=200,
    facet_kws={
        'sharex': True,
        'sharey': False,
        'legend_out': True
    })
leg = g._legend
leg._loc = 7
# leg.set_bbox_to_anchor([1, 0.1])
g.figure.set_figheight(12)
g.figure.set_figwidth(30)
g.figure.tight_layout()

for a in g.axes_dict:
    g.axes_dict[a].yaxis.grid(color='black', linestyle=(0, (5, 10)), linewidth=0.5)
    g.axes_dict[a].xaxis.grid(color='black', linestyle=(0, (5, 10)), linewidth=0.5)

g.figure.savefig("all_pvalues_bonita.pdf")
g.figure.savefig("all_pvalues_bonita.png", dpi=300)