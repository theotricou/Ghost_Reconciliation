#!/usr/bin/env python3
# by Theo Tricou

from ete3 import Tree as tr
import sys
import glob
import re

tool = sys.argv[1]

ale = tr("aletree", format = 1)
ext = tr("SAMPLE_1/SampledSpeciesTree.nwk", format = 1)

map = dict()

map[ale.name] = [ext.name, [], []]

for i in ext:
    ext_node = i
    ale_node = ale.search_nodes(name=i.name)[0]
    while ale_node.name not in map.keys() :
        map[ale_node.name] = [ext_node.name, [], []]
        ext_node = ext_node.up
        ale_node = ale_node.up

import glob


uTs = glob.glob("*" + tool +".uTs")
for i in uTs:
    with open(i) as file:
        for lines in file:
            line = re.split("\t| ", lines.strip())
            map[line[0]][1].append(float(line[2]))
            map[line[1]][2].append(float(line[2]))

map_out=open("Nodes_mapping_" + tool + ".txt","w+")
map_out.write(" ".join(["Ale_ID", "Node", "From-sum", "To_sum"])+"\n")
for i in sorted(map.keys()):
    map_out.write(i + " " + map[i][0] + " " + str(round(sum(map[i][1]), 6)) + " " + str(round(sum(map[i][2]), 6)) +"\n")

map_out.close()


# GNU Ghost
