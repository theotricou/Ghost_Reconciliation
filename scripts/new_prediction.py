#!/usr/bin/env python3
# by Theo Tricou

from ete3 import Tree as tr

def get_ghost_prediction(com_node):
    # com_node = com.search_nodes(name = node.name)[0]
    while com_node.name not in back_bone_nodes:
        com_node = com_node.up
    ext_node = bak_bone_tree.search_nodes(name = com_node.name)[0]
    while len(ext_node.get_children()) not in [0,2]:
        ext_node = ext_node.get_children()[0]
    return(ext_node.name)

def get_ghost_banch_length(node):
    if node.name in nodes_contemporary:
        return(node.dist)
    else:
        dist_to_root = node.get_distance(com)
        distup_to_root = node.up.get_distance(com)
        if dist_to_root > extroot_to_comroot and distup_to_root > extroot_to_comroot:
            ext_bl = round(node.dist, 6)
        elif dist_to_root > extroot_to_comroot and distup_to_root < extroot_to_comroot:
            ext_bl = round((dist_to_root - distup_to_root) - (extroot_to_comroot - distup_to_root), 6)
        else:
            ext_bl = 0
        return(ext_bl)

com = tr("T/CompleteTree.nwk", format = 1)
ext = tr("SAMPLE_1/SampledSpeciesTree.nwk", format = 1)
ext.dist = 0
back_bone_nodes = []
back_bone_nodes.append(com.name)

for i in ext:
    node = com.search_nodes(name = i.name)[0]
    while node.name != com.name:
        back_bone_nodes.append(node.name)
        node = node.up

bak_bone_tree = com.copy()
bak_bone_tree.prune(list(set(back_bone_nodes)))
# print(bak_bone_tree.get_ascii(show_internal=True))

anc = com.search_nodes(name=ext.name)[0]
extroot_to_comroot = anc.get_distance(com)

nodes_contemporary = [node.name for node in anc.iter_descendants()]

mapping = dict()
stat_by_node = dict()

for i in com.traverse():
    prediction = get_ghost_prediction(i)
    if not i.is_root():
        i.add_features(ghost_prediction = prediction, ghost_dist = get_ghost_banch_length(i))
    else:
        i.add_features(ghost_prediction = prediction, ghost_dist = 0)
    if i.name in list(set(back_bone_nodes)):
        mapping[i.name] = [prediction, prediction]
    else:
        mapping[i.name] = [prediction, "Ghost"]
    if prediction in stat_by_node:
        stat_by_node[prediction][1] += 1
        stat_by_node[prediction][2] += i.ghost_dist
    else:
        stat_by_node[prediction] = [ext.search_nodes(name=prediction)[0].dist, 1, i.ghost_dist, 0, 0]


list_tr_mapped = dict()

with open("Transfers_sim.txt") as file:
    for lines in file:
        line = lines.strip().split(" ")
        node_from = mapping[line[2]][0]
        node_to   = mapping[line[3]][0]
        if line[0] not in list_tr_mapped.keys():
            list_tr_mapped[line[0]] = [[node_from, node_to]]
        elif [node_from, node_to] not in list_tr_mapped[line[0]]:
            list_tr_mapped[line[0]].append([node_from, node_to])

for i in list_tr_mapped:
    for j in list_tr_mapped[i]:
        if j[0] != j[1]:
            if j[0] != ext.name:
                if ext.search_nodes(name = j[0])[0].get_sisters()[0].name != j[1]:
                    stat_by_node[j[0]][3] += 1
                    stat_by_node[j[1]][4] += 1
            else:
                stat_by_node[j[0]][3] += 1
                stat_by_node[j[1]][4] += 1


map_out=open("Nodes_mapping_sim.txt","w+")
map_out.write(" ".join(["Node", "Ghost_donor", "Ghost_recipient"])+"\n")
for i in sorted(mapping.keys()):
    map_out.write(" ".join([i, mapping[i][0], mapping[i][1]])+"\n")

map_out.close()

stat_out=open("Stats_simulation.txt","w+")
stat_out.write(" ".join(["Node", "branch_length", "N_ghost_node", "L_ghost_branch", "N_transfers_donor", "N_transfers_recip"])+"\n")
for i in sorted(stat_by_node.keys()):
    stat_out.write(i + " " + " ".join(str(round(v, 6)) for v in stat_by_node[i])+"\n")

stat_out.close()


tree = tr("T/CompleteTree.nwk", format = 1)

for i in tree.traverse():
    i.add_features(tr_to=dict(), tr_from=dict())



# creates a dictionnary with all events
tr_events = dict()
for lines in reversed(list(open("Transfers_sim.txt"))):
    line = lines.strip().split(" ")
    if line[0] in tr_events.keys() and line[2] in tr_events[line[0]].keys() and line[3] in tr_events[line[0]][line[2]].keys():
        tr_events[line[0]][line[2]][line[3]].append(line[1])
    elif line[0] not in tr_events.keys():
        recipient = dict()
        donor = dict()
        recipient[line[3]] = [line[1]]
        donor[line[2]] = recipient
        tr_events[line[0]] = donor
    elif line[2] not in tr_events[line[0]].keys():
        recipient = dict()
        recipient[line[3]] = [line[1]]
        tr_events[line[0]][line[2]] = recipient
    elif line[3] not in tr_events[line[0]][line[2]].keys():
        tr_events[line[0]][line[2]][line[3]] = [line[1]]

# remove events taht are received by ghost with no subsequent transfers from
# the ghost or on of its iter_descendants

gene = 0
for lines in reversed(list(open("Transfers_sim.txt"))):
    line = lines.strip().split(" ")
    if line[3] not in back_bone_nodes and gene != line[0]:
        tr_events[line[0]][line[2]][line[3]].remove(line[1])
        if len(tr_events[line[0]][line[2]][line[3]]) == 0:
            del tr_events[line[0]][line[2]][line[3]]
        if len(tr_events[line[0]][line[2]]) == 0:
            del tr_events[line[0]][line[2]]
        if len(tr_events[line[0]]) == 0:
            del tr_events[line[0]]
        print(line)
    else:
        gene = line[0]
        if line[3] not in back_bone_nodes:
            recip_node = tree.search_nodes(name = line[3])[0]
            for donor_node in tr_events[line[0]]:
                if line[3] in tr_events[line[0]][donor_node].keys():
                    for time tr_events[line[0]][donor_node][line[3]]:
                        if time > line[1]:
                            for node_recip2 in tr_events[line[0]][line[3]]:


import pandas as pd

dd = pd.read_csv("Transfers_sim.txt", sep=" ", header = None, names = ["gene","time","from","to"])
dd["from"]  = dd["from"].astype("string")
dd["to"]  = dd["to"].astype("string")
dd["to_ghost"] = True
dd.loc[dd.to.isin(back_bone_nodes), "to_ghost"] = False
dd["from_prediction"] = None
dd["to_prediction"] = None


# remove first transfers to ghost with no impact
for i in sorted(list(dd.gene.unique())):
    sub = dd[dd["gene"] == i].sort_values(['time'], ascending=[False])
    for index, row in sub.iterrows():
        is_True = row["to_ghost"]
        if is_True == True:
            dd = dd.drop(index)
        else:
            break
# donor prediction
for i in sorted(list(dd.gene.unique())):
    sub = dd[dd["gene"] == i].sort_values(['time'], ascending=[False])
    for index, row in sub.iterrows():
        if not row["to_ghost"]:
            fr_node = tree.search_nodes(name = row["from"])[0]




fr_node = tree.search_nodes(name = "10055")[0]
node = fr_node
time = sub[0:1]["time"].item()
true_donor = "not_found"
while not len(sub[sub["to"] == node.name]) == 0 or node.name not in back_bone_nodes or true_donor == "not_found":
    print(node.name)
    print(time)
    previous_node = node
    if len(sub[sub["to"] == node.name].query('time < @time')) != 0:
        closest_transfer = sub[sub["to"] == node.name].query('time < @time')[:1]
        time = closest_transfer["time"].item()
        node = tree.search_nodes(name = closest_transfer["from"].item())[0]
    elif node.name not in back_bone_nodes:
        node = node.up
    elif node.name in back_bone_nodes:
        child = node.get_children()[1]
        if child == previous_node:
            child = node.get_children()[0]
        # sub.loc[index, "from_prediction"] = mapping[node.name][0]
        # true_donor = mapping[node.name][0]
        all_descendants = [child.name] + [node.name for node in child.get_descendants()]
        bone_all_descendants = [x for x in all_descendants if x in back_bone_nodes]
        for i in bone_all_descendants:
            if len(sub[sub["from"] == i]) >= 1:
                true_donor = i


sub = dd[dd["gene"] == 9].sort_values(['time'], ascending=[False])
for index, row in sub.iterrows():
    if not row["to_ghost"]:
        fr_node = tree.search_nodes(name = row["from"])[0]
        node = fr_node
        time = row["time"]
        true_donor = "not_found"
        while not len(sub[sub["to"] == node.name]) == 0 or node.name not in back_bone_nodes or true_donor == "not_found":
            print(node.name)
            print(time)
            previous_node = node
            if len(sub[sub["to"] == node.name].query('time < @time')) != 0:
                closest_transfer = sub[sub["to"] == node.name].query('time < @time')[:1]
                time = closest_transfer["time"].item()
                node = tree.search_nodes(name = closest_transfer["from"].item())[0]
            elif node.name not in back_bone_nodes:
                node = node.up
            elif node.name in back_bone_nodes:
                child = node.get_children()[1]
                if child == previous_node:
                    child = node.get_children()[0]
                # sub.loc[index, "from_prediction"] = mapping[node.name][0]
                # true_donor = mapping[node.name][0]
                all_descendants = [child.name] + [node.name for node in child.get_descendants()]
                bone_all_descendants = [x for x in all_descendants if x in back_bone_nodes]




print(time)
print(node.name)
print(previous_node)





for index, row in dd.iterrows():
    if row["to_ghost"]:
        to_node = tree.search_nodes(name = row["to"])[0]
        list_descendants = [node.name for node in to_node.get_descendants()]
        list_descendants.append(to_node.name)
        for i in to_node.traverse("preorder"):
            node.name

to_node = tree.search_nodes(name = "10006")[0]
list_descendants = [node.name for node in to_node.get_descendants()]
list_descendants.append(to_node.name)
for node in to_node.traverse("preorder"):
    if



for i in sorted(list(dd.gene.unique())):
    sub = dd[dd["gene"] == i].sort_values(['time'], ascending=[False])
    for index, row in sub.iterrows():



dd[dd["gene"] == 9]


tr_events = dict
with reversed(list(open("Transfers_sim.txt"))) as file:
    for lines in file:
        line = lines.strip().split(" ")
        to_node = tree.search_nodes(name = line[2])[0]
        from_node = tree.search_nodes(name = line[3])[0]

        to_node.tr_to[line[1]] = from_node.name
        from_node.tr_from[line[1]] = to_node.name





        node_from = mapping[line[1]][0]
        node_to   = mapping[line[2]][0]
        if line[0] not in list_tr_mapped.keys():
            list_tr_mapped[line[0]] = [[node_from, node_to]]
        elif [node_from, node_to] not in list_tr_mapped[line[0]]:
            list_tr_mapped[line[0]].append([node_from, node_to])


for i in tree.traverse():






# GNU Ghost
