#!/usr/bin/env python3
# by Theo Tricou



from ete3 import Tree as tr


ale = tr("aletree", format = 1)
ext = tr("SAMPLE_1/SampledSpeciesTree.nwk", format = 1)

map = dict()

map[ale.name] = [ext.name, 0, 0]

for i in ext:
    ext_node = i.copy()
    ale_node = ale.search_nodes(name=i.name)[0]
    while ext_node.name not in map.keys() :
        map[ale_node.name] = [ext_node.name, [], []]
        ext_node = ext_node.up
        ale_node = ale_node.up


#
# import glob
#
# uTs = glob.glob("*uTs")
# for i in uTs:
#     with open(i) as file:
#         for lines in file:
#             if not lines.startswith("#"):
#                 line = lines.strip().split(" ")
#                 map[line[1]][1] += line[3]
#                 map[line[2]][1] += line[3]
#



map_out=open("Nodes_mapping_ale.txt","w+")
map_out.write(" ".join(["Ale_ID", "Node", "From-sum", "To_sum"])+"\n")
for i in sorted(map.keys()):
    # map_out.write(i + " " + map[i][0] + " " + map[i][1] + " " + map[i][2] +"\n")
    map_out.write(i + " " + map[i][0] + "\n")

map_out.close()


# GNU Ghost
