from moBonita_kegg_parser import *
#find_pathways_kegg()

from glob import glob
import networkx as nx
import re

pathways = glob("*_before.graphml")
gmt = {}

def getPathwayName(hsaURL):
    fileReg = re.compile("NAME\s+(\w+.*)")
    pathwayFile = requests.get("http://rest.kegg.jp/get/" + hsaURL, stream=True)
    for line in pathwayFile.iter_lines():
        line = line.decode("utf-8")
        result = fileReg.match(line)
        if result:
            return result.group(1)
    print(hsaURL)
    return hsaURL

regex1 = re.compile(r'.*cpd:.*')
regex2 = re.compile(r'.*path:.*')
regex3 = re.compile(r'E\d\.\.*') 
regex4 = re.compile(r'gl:.*')
regex5 = re.compile(r'.*\_.*') 
regex6 = re.compile(r'.*\_.*')
regex7 = re.compile(r'K\d\.\.*')
regex8 = re.compile(r'.*dr:.*')
f = open("kegg_networks.gmt", "w")
for p in pathways:
    graph = nx.read_graphml(p)
    filtered = [i for i in graph.nodes() if not regex1.match(i)]
    filtered = [i for i in filtered if not regex2.match(i)]
    filtered = [i for i in filtered if not regex3.match(i)]
    filtered = [i for i in filtered if not regex4.match(i)]
    filtered = [i for i in filtered if not regex6.match(i)]
    filtered = [i for i in filtered if not regex7.match(i)]   
    filtered = [i for i in filtered if not regex8.match(i)]   
    gmt[p[0:8]] = []
    for i in filtered:
        if regex5.match(i):
            i = i.split("_")
            prefix = i[0][:-1]
            if i[len(i)-1].isdigit() and i[0][-1].isdigit():
                temp = [prefix+str(j) for j in range(int(i[0][-1]), int(i[len(i)-1])+1)]
                gmt[p[0:8]].extend(temp)
            else:
                temp = [prefix+i[j] for j in range(1, len(i))]
                temp.append(prefix+i[0][-1])
                gmt[p[0:8]].extend(temp)
        else:
            gmt[p[0:8]].append(i)
    f.write("\t".join([getPathwayName(str(p[0:8])), str(p[0:8]),"\t".join(gmt[p[0:8]])]))
    f.write("\n")
    print("\t".join([getPathwayName(str(p[0:8])), str(p[0:8]),"\t".join(gmt[p[0:8]])]))
f.close()
