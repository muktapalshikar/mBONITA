df = pd.read_csv("local1Data.csv", index_col=0)
df['indegree'] = [
    len(set([item for sublist in literal_eval(i) for item in sublist]))
    for i in df.andNodeList
]
df = df.replace('concatenated', 'Concatenated')
df.head()


# %%
# Get optimized networks
def maxERS(inDegree):
    if inDegree == 0:
        maxERS = 1
    else:
        inDegree = min(inDegree, 3)
        if inDegree == 2:
            maxERS = 15
        else:
            maxERS = 127
    return maxERS


df.loc[:, 'maxERS'] = [maxERS(i) for i in df.indegree]
df.loc[:, 'hasReducedERS'] = df.loc[:, 'ruleLengths'] < df.loc[:, 'maxERS']


# %%
commonNodesRules = df[df.nodeList.isin(
    nt2.index.get_level_values('Display Name').tolist())]
netid_to_name = {
    key: value
    for (key, value) in
    zip(set(nt2.index.get_level_values('Network')),
        [getPathwayName(f) for f in set(nt2.index.get_level_values('Network'))])
}
nt2['networkName'] = [netid_to_name[i] for i in nt2.index.get_level_values('Network')]
commonNodesRules = commonNodesRules[commonNodesRules.networkName.isin(
    nt2.networkName)].reset_index(drop=True)
commonNodesRules.head()

# %%
cnr_transcriptomics = commonNodesRules[
    commonNodesRules['methodName'].str.contains('Transcriptomics')]
print(cnr_transcriptomics.shape)
cnr_phosphoproteomics = commonNodesRules[
    commonNodesRules['methodName'].str.contains('Phosphoproteomics')]
print(cnr_phosphoproteomics.shape)
cnr_proteomics = commonNodesRules[commonNodesRules['methodName'].str.match(
    'Proteomics'
)]  #note difference between contains and match - DO NOT CHANGE THIS
print(cnr_proteomics.shape)

# %%
tallyDF_final = commonNodesRules
tallyDF_final["percentOverlap"] = [np.nan for i in range(0, tallyDF_final.shape[0])]
print(tallyDF_final.shape)

# %%
def percentOverlap2(equivs, exactCompare=False):
    """Compare n lists to get a percentage overlap compared to the "universe" (i.e. the union of all n lists)"""
    overlap = set()
    universe = set()
    for l in equivs:
        rule = set([tuple(i) for i in l])
        universe = universe | rule
        if len(overlap) == 0:
            overlap = rule
        else:
            overlap = rule & overlap
    if len(universe) == 0:
        return(np.nan)
    else:
        if exactCompare:
            percentOverlap = float(len(overlap)) / len(universe) *100
        else:
            if len(overlap) > 0:
                return(float(1))
            else:
                return(float(0))
    return(percentOverlap)

# %%
def getNodeMinimalRule(ruleSet, andNodeList):
    numOrRules = [sum(ruleSet[j]) for j in range(len(ruleSet))]
    deciderVar = max(numOrRules)  # maximal or rules
    maxOrRules = ruleSet[
        numOrRules.index(deciderVar)
    ]  # rules with maximum or terms
    maxUpstreamNodes = 0
    minimalRule = []
    for orRule in [maxOrRules]:
        if sum(orRule) > 0:
            numUpstreamNodes = [
                andNodeList[orTerm]
                for orTerm in range(len(orRule))
                if orRule[orTerm] == 1
            ]
        else:
            minimalRule = orRule
            continue
        numUpstreamNodes = [len(element) for element in numUpstreamNodes]
        numUpstreamNodes = sum(numUpstreamNodes)
        if numUpstreamNodes > maxUpstreamNodes:
            maxUpstreamNodes = numUpstreamNodes
            minimalRule = orRule
        else:
            maxUpstreamNodes = maxUpstreamNodes
            minimalRule = minimalRule
    return minimalRule

# %%
for network in list(set(commonNodesRules.networkName))[0:15]:
    print("############")
    print(network.upper())
    print("######\n")
    net_trans = cnr_transcriptomics[
        cnr_transcriptomics['networkName'].str.contains(network, regex=False)]
    net_phosph = cnr_phosphoproteomics[
        cnr_phosphoproteomics['networkName'].str.contains(network,
                                                          regex=False)]
    net_prot = cnr_proteomics[cnr_proteomics['networkName'].str.contains(
        network, regex=False)]
    allNodes = pd.concat(
        [net_trans.nodeList, net_phosph.nodeList, net_prot.nodeList])
    for node in set(allNodes):

        tempTrans = pd.DataFrame(net_trans[net_trans["nodeList"] == node])
        tempProt = pd.DataFrame(net_prot[net_prot["nodeList"] == node])
        tempPhosph = pd.DataFrame(net_phosph[net_phosph["nodeList"] == node])
        #print(tempTrans.shape, tempPhosph.shape, tempProt.shape)
        """
        percentOverlap = percentOverlap2([
            tempTrans.plainRules.tolist(),
            tempProt.plainRules.tolist(),
            tempPhosph.plainRules.tolist()
        ])
        if percentOverlap != 1:
            print(node, network)
            print([
                set(tempTrans.plainRules.tolist()),
                set(tempProt.plainRules.tolist()),
                set(tempPhosph.plainRules.tolist())
            ], percentOverlap)
        """
        minimalTrans = getNodeMinimalRule(
            literal_eval(tempTrans.equivs.iloc[0]),
            literal_eval(tempTrans.andNodeList.iloc[0]))
        minimalTrans = literal_eval(tempTrans.plainRules.iloc[0])[literal_eval(
            tempTrans.equivs.iloc[0]).index(minimalTrans)]
        minimalProt = getNodeMinimalRule(
            literal_eval(tempProt.equivs.iloc[0]),
            literal_eval(tempProt.andNodeList.iloc[0]))
        minimalProt = literal_eval(tempProt.plainRules.iloc[0])[literal_eval(
            tempProt.equivs.iloc[0]).index(minimalProt)]
        minimalPhosph = getNodeMinimalRule(
            literal_eval(tempPhosph.equivs.iloc[0]),
            literal_eval(tempPhosph.andNodeList.iloc[0]))
        minimalPhosph = literal_eval(tempPhosph.plainRules.iloc[0])[literal_eval(
            tempPhosph.equivs.iloc[0]).index(minimalPhosph)]
        print(node)
        print("Trans: ", minimalTrans, "\nProt:", minimalProt, "\nPhosph:",
              minimalPhosph)
        print("######\n")

# %%
minimalTrans = "RB1*=(CDKN1A) or (CDK2 and CDKN1A) or (CDK2)"
minimalProt = "RB1*=(GRB2) or (TP53 and CDK2 and GRB2) or (TP53 and CDK2) or (TP53 and GRB2) or (TP53) or (CDK2 and GRB2) or (CDK2)"
minimalPhosph = "RB1*=(not PTK2B) or (not TP53 and  not MAP3K7 and  not PTK2B) or (not TP53 and  not MAP3K7) or (not TP53) or (not MAP3K7 and  not PTK2B) or (not MAP3K7)"

# simple intersection
def makeRuleList(rule):
    line=rule
    line=line.strip()
    matcher=re.search("(.*)\*=(.*)", line)
    node=matcher.group(1).strip()
    rule=matcher.group(2).strip()
    preds=re.sub("and|not|or|\(|\)", "", rule)
    preds=preds.strip()
    preds=re.sub("\s+", ",", preds)
    preds=tuple(preds.split(","))
    rule=rule.split(" or ")
    return rule, preds
    
transRule, transPreds = makeRuleList(minimalTrans)
protRule, protPreds = makeRuleList(minimalProt)
phosphRule, phosphPreds = makeRuleList(minimalPhosph)

print(set(transRule),set(protRule),set(phosphRule))
print(set(transRule) & set(protRule) & set(phosphRule))

# %%
from utils import Get_expanded_network
Get_expanded_network([literal_eval(t)[0] for t in net_trans.plainRules.tolist()], equal_sign="*=").nodes


# %%
for network in list(set(commonNodesRules.networkName)):
    net_trans = cnr_transcriptomics[cnr_transcriptomics['networkName'].str.contains(network, regex=False)]
    net_phosph = cnr_phosphoproteomics[cnr_phosphoproteomics['networkName'].str.contains(network, regex=False)]
    net_prot = cnr_proteomics[cnr_proteomics['networkName'].str.contains(network, regex=False)]
    print(net_trans.shape)
    for node in set(df.nodeList):
        tempTrans = pd.DataFrame(net_trans[net_trans["nodeList"] == node])
        tempProt = pd.DataFrame(net_prot[net_prot["nodeList"] == node])
        tempPhosph = pd.DataFrame(net_phosph[net_phosph["nodeList"] == node])
        if (net_trans.shape[0] >= 0 and net_prot.shape[0] >= 0 and net_phosph.shape[0]>= 0):
            percentOverlap = percentOverlap2([tempTrans.plainRules.tolist(), tempProt.plainRules.tolist(), tempPhosph.plainRules.tolist()])
            temp2 = [t1 & t2 for t1, t2 in zip(list(tallyDF_final["networkName"]== network), list(tallyDF_final["nodeList"] == node))]
            tallyDF_final.loc[temp2, "percentOverlap"] = percentOverlap

# %%
tallyDF_final.head()
tallyDF_final.to_csv("tallyDF.csv")
"""Aggregate by network and show average percent overlap"""
temp1 = pd.DataFrame(tallyDF_final.groupby(["networkName"]).mean())
temp1.head()


