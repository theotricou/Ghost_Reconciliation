#!/usr/bin/env python3
# by Theo Tricou



from ete3 import Tree as tr
import glob, os



def all_mapping(dir):


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
    ext = tr(os.path.join(dir, "SampledSpeciesTree.nwk"), format = 1)

    ext.dist = 0
    back_bone_nodes = []
    back_bone_nodes.append(com.name)

    for i in ext:
        node = com.search_nodes(name = i.name)[0]
        while node.name != com.name:
            back_bone_nodes.append(node.name)
            node.name
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
            if i.name not in bak_bone_tree:
                stat_by_node[prediction][1] += 1
                stat_by_node[prediction][2] += i.ghost_dist
        else:
            stat_by_node[prediction] = [ext.search_nodes(name=prediction)[0].dist, 0, 0, 0, 0]


    list_tr_mapped = dict()

    with open("Transfers_sim.txt") as file:
        for lines in file:
            line = lines.strip().split(" ")
            node_from = mapping[line[2]][0]
            node_to   = mapping[line[3]][1]
            if line[0] not in list_tr_mapped.keys():
                list_tr_mapped[line[0]] = [[node_from, node_to]]
            elif [node_from, node_to] not in list_tr_mapped[line[0]]:
                list_tr_mapped[line[0]].append([node_from, node_to])


    for i in list_tr_mapped:
        for j in list_tr_mapped[i]:
            if j[1] != 'Ghost' and j[0] != j[1]:
                if j[0] != ext.name:
                    if ext.search_nodes(name = j[0])[0].get_sisters()[0].name != j[1]:
                        stat_by_node[j[0]][3] += 1
                        stat_by_node[j[1]][4] += 1
                else:
                    stat_by_node[j[0]][3] += 1
                    stat_by_node[j[1]][4] += 1




    sampling = dir.split("_")[1]

    map_out=open("Nodes_mapping_" + sampling + "_sim.txt","w+")
    map_out.write(" ".join(["Node", "Ghost_donor", "Ghost_recipient"])+"\n")
    for i in sorted(mapping.keys()):
        map_out.write(" ".join([i, mapping[i][0], mapping[i][1]])+"\n")

    map_out.close()

    stat_out=open("Stats_simulation_" + sampling + ".txt","w+")
    stat_out.write(" ".join(["Node", "branch_length", "N_ghost_node", "L_ghost_branch", "N_transfers_donor", "N_transfers_recip"])+"\n")
    for i in sorted(stat_by_node.keys()):
        stat_out.write(i + " " + " ".join(str(round(v, 6)) for v in stat_by_node[i])+"\n")

    stat_out.close()


dir = glob.glob("Sample_*")

for i in dir:
    all_mapping(i)



# GNU Ghost
