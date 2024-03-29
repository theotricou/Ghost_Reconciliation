#!/usr/bin/env python3
# by Theo Tricou

from ete3 import Tree as tr
import os
import pandas as pd
import glob

def build_detach(tc, ts):
    # create a detach gene tree
    sample_spe = ts.get_leaf_names()
    to_keep = []
    for i in sample_spe:
    	node = tc.search_nodes(name = i)[0]
    	while node:
    		to_keep.append(node.name)
    		node = node.up
    to_keep = list(set(to_keep))
    tc_detach = tc.copy()
    for i in tc_detach.iter_descendants("postorder"):
    	if i.name not in to_keep:
    	    node = tc_detach.search_nodes(name = i.name)[0]
    	    removed_node = node.detach()
    return(tc_detach)

def are_they_sisters(to_name, from_name, tree):
    # print(from_name)
    to_node = tree.search_nodes(name = to_name)[0]
    from_node = tree.search_nodes(name = from_name)[0]
    if to_node.is_root():
        return(False)
    else:
        if to_node.get_sisters()[0] == from_node:
            return(False)
        else:
            return(True)

def get_prediction_R(dd, detach_tree, tr_prediction, samp):
    # print(dd["to"] + "_" + dd["to_id"])
    try:
        R_node = detach_tree.search_nodes(name = dd["to"] + "_" + dd["to_id"])[0]
        while len(R_node.get_children()) not in [0,2]:
            R_node = R_node.get_children()[0]
        list_recip = list(tr_prediction["to_prediction"] + "_" + tr_prediction["to_prediction_id"])
        nodes_recip = [i for i in R_node.traverse() if i.name in list_recip]
        leaves_in_nodes_recip = []
        for i in nodes_recip:
            for l in i:
                leaves_in_nodes_recip.append(l.name)
        node_no_tr = list(set(R_node.get_leaf_names()) - set(leaves_in_nodes_recip))
        if len(node_no_tr) == 1:
            R_node = detach_tree.search_nodes(name = node_no_tr[0])[0]
        elif len(node_no_tr) == 0:
            return(["TO_REMOVE", "TO_REMOVE"])
        else:
            R_node = detach_tree.get_common_ancestor(node_no_tr)
        # node_samp = [i.name for i in samp.traverse()]
        # for i in R_node.traverse():
        #     if i.name.split("_")[0] in node_samp and i.name not in list_recip:
        #         R_node = i
        #         break
        if R_node.name.split("_")[0] == dd["to"] and R_node.name.split("_")[1] > dd["to_id"]:
            R_node = detach_tree.search_nodes(name = dd["to"] + "_" + dd["to_id"])[0]
        return(R_node.name.split("_"))
    except IndexError:
        return(["TO_REMOVE", "TO_REMOVE"])



def make_list(node, speciation, tr_prediction):
    name = node.name.split("_")[0]
    child1 = speciation.loc[speciation.c1 == name]
    child2 = speciation.loc[speciation.c2 == name]
    spe_node_up = pd.concat([child1, child2])
    if len(spe_node_up) == 0:
        spe_node_up = speciation.loc[speciation.up == name]
    spe_list_recip = list(tr_prediction[tr_prediction["TIME"] >= spe_node_up.iloc[0]['TIME']]["to_prediction"] + "_" + tr_prediction[tr_prediction["TIME"] >= spe_node_up.iloc[0]['TIME']]["to_prediction_id"])
    return(spe_list_recip)


def get_prediction_D(dd, detach_tree, tr_prediction, samp, speciation):
    if dd["sister"] == False:
        return(["TO_REMOVE", "TO_REMOVE"])
    else:
        # dd=tr_prediction.loc[38]
        R_node = detach_tree.search_nodes(name = dd["to_prediction"] + "_" + dd["to_prediction_id"])[0]
        D_node = detach_tree.search_nodes(name = R_node.name)[0]
        list_recip = list(tr_prediction[tr_prediction["TIME"] >= dd["TIME"]]["to_prediction"] + "_" + tr_prediction[tr_prediction["TIME"] >= dd["TIME"]]["to_prediction_id"])
        list_recip_all = list(tr_prediction["to_prediction"] + "_" + tr_prediction["to_prediction_id"])
        Root_nodes = []
        root_node = detach_tree.get_common_ancestor(detach_tree.get_leaf_names())
        while root_node:
            Root_nodes.append(root_node.name)
            root_node = root_node.up
        if D_node.name in Root_nodes:
            return(["TO_REMOVE", "TO_REMOVE"])
        while len(D_node.get_sisters()) != 1:
            D_node = D_node.up
        D_node = D_node.get_sisters()[0]
        if D_node.up.name not in Root_nodes:
            leaf_in_D_node = D_node.get_leaf_names()
            nodes_recip = [i for i in D_node.traverse() if i.name in make_list(D_node, speciation, tr_prediction)]
            leaves_in_nodes_recip = []
            for i in nodes_recip:
                for l in i:
                    leaves_in_nodes_recip.append(l.name)
            count = 0
            while set(leaf_in_D_node) == set(leaves_in_nodes_recip) and D_node.name not in Root_nodes:
                count += 1
                if D_node.up.name in Root_nodes:
                    D_node = D_node.get_sisters()[0]
                    while len(D_node.get_children()) not in [0,2]:
                        D_node = D_node.get_children()[0]
                else:
                    D_node = D_node.up
                try: leaf_in_D_node = D_node.get_leaf_names()
                except:
                    return(["TO_REMOVE", "TO_REMOVE"])
                # print(D_node.name)
                nodes_recip = [i for i in D_node.traverse() if i.name in make_list(D_node, speciation, tr_prediction)]
                leaves_in_nodes_recip = []
                for i in nodes_recip:
                    for l in i:
                        leaves_in_nodes_recip.append(l.name)
        while len(D_node.get_children()) not in [0,2]:
            D_node = D_node.get_children()[0]
        leaf_in_D_node = D_node.get_leaf_names()
        nodes_recip = [i for i in D_node.traverse() if i.name in make_list(D_node, speciation, tr_prediction)]
        leaves_in_nodes_recip = []
        for i in nodes_recip:
            for l in i:
                leaves_in_nodes_recip.append(l.name)
        node_no_tr = list(set(leaf_in_D_node) - set(leaves_in_nodes_recip))
        if len(node_no_tr) == 1:
            D_node = detach_tree.search_nodes(name = node_no_tr[0])[0]
        elif len(node_no_tr) > 1:
            D_node = D_node.get_common_ancestor(list(set(leaf_in_D_node) - set(leaves_in_nodes_recip)))
        if D_node.name.split("_")[0] not in [i.name for i in samp.traverse()]:
            D_node = com.search_nodes(name = D_node.name.split("_")[0])[0]
            prediction = com.search_nodes(name = get_ghost_prediction(D_node, back_bone_nodes, bak_bone_tree))[0]
            return([prediction.name, "X"])
        else:
            return(D_node.name.split("_"))


def complete_prediction(gene, dir, back_bone_nodes, samp):
    num = gene.split("/")[1].split("_")[0]
    tr_data = pd.read_csv("G/Gene_families/" + num + "_events.tsv", sep="\t")
    speciation = tr_data.loc[tr_data["EVENT"] == "S"].copy()
    speciation[["up", "up_id", "c1", "c1_id", "c2", "c2_id"]] = speciation["NODES"].apply(lambda x: pd.Series(x.split(";")))
    tr_data = tr_data[tr_data["EVENT"] == "T"]
    if tr_data.empty:
        return
    tr_data[["from", "from_id", "from_d", "from_d_id", "to", "to_id"]] = tr_data["NODES"].apply(lambda x: pd.Series(x.split(";")))
    # tr_data = tr_data.drop(columns = "NODES")
    tr_prediction = tr_data.loc[tr_data.to.isin(back_bone_nodes)].copy()
    if tr_prediction.empty:
        return
    # read genes trees
    t_gene_p = tr('G/Gene_trees/' + num + '_completetree.nwk',format = 1)
    t_gene_s = tr(gene ,format = 1)
    detach_tree = build_detach(t_gene_p, t_gene_s)
    # print(detach_tree.get_ascii(attributes=["name"], show_internal=True))
    # remove duplicate recipient
    tr_prediction = tr_prediction.sort_values('TIME', ascending=False).drop_duplicates(['to']).sort_values('TIME')
    # remove transfers beween sister species in the complete tree, if any:
    # tr_prediction["sister"] = tr_prediction.apply(lambda x: are_they_sisters(x["to"], x["from"], com), axis = 1)
    # tr_prediction = tr_prediction[tr_prediction.sister]
    # predict recipient for the remaining transfers
    tr_prediction[["to_prediction", "to_prediction_id"]] = ["None", "None"]
    for index, row in tr_prediction.sort_values('TIME', ascending=False).iterrows():
        tr_prediction.loc[index, ["to_prediction", "to_prediction_id"]] = get_prediction_R(row, detach_tree, tr_prediction, samp)
    # tr_prediction[["to_prediction", "to_prediction_id"]] = tr_prediction.apply(lambda x: pd.Series(get_prediction_R(x, detach_tree, tr_prediction, samp)), axis = 1)
    DandR = ["None", "None"]
    for index, row in tr_prediction.iterrows():
        if set(tr_prediction.loc[index, ["from", "to"]]) == set(DandR):
            tr_prediction.loc[index, ["to_prediction"]] = "TO_REMOVE"
        DandR = list(tr_prediction.loc[index, ["from", "to"]])
    tr_prediction = tr_prediction[tr_prediction.to_prediction != "TO_REMOVE"]
    if tr_prediction.empty:
        return
    # predict donor
    tr_prediction[["from_prediction", "from_prediction_id"]] = ["None", "None"]
    tr_prediction["sister"] = tr_prediction.apply(lambda x: are_they_sisters(x["to"], x["from"], com), axis = 1)
    tr_prediction[["from_prediction", "from_prediction_id"]] = tr_prediction.apply(lambda x: pd.Series(get_prediction_D(x, detach_tree, tr_prediction, samp, speciation)), axis = 1)
    tr_prediction = tr_prediction[tr_prediction.from_prediction != "TO_REMOVE"]
    if tr_prediction.empty:
        return
    # once again we remove transfers beween sister species, in sampled tree this time
    tr_prediction["sister"] = tr_prediction.apply(lambda x: are_they_sisters(x["to"], x["from"], com), axis = 1)
    tr_prediction = tr_prediction[tr_prediction.sister]
    if tr_prediction.empty:
        return
    tr_prediction["sister"] = tr_prediction.apply(lambda x: are_they_sisters(x["to_prediction"], x["from_prediction"], samp), axis = 1)
    tr_prediction = tr_prediction[tr_prediction.sister]
    if tr_prediction.empty:
        return
    tr_prediction["Gene"] = gene.split("/")[1].split("_")[0]
    tr_prediction = tr_prediction.drop(columns=["EVENT", "NODES", "sister"])
    name_out = gene.split(".")[0]
    tr_prediction.to_csv(name_out + "_prediction", sep=' ', index = False)
    return(tr_prediction)

def is_to_strict_descendant(name_from ,name_to, tree):
    n_from = tree.search_nodes(name = name_from)[0]
    n_to = tree.search_nodes(name = name_to)[0]
    if n_to in n_from.get_descendants():
        return(True)
    else:
        return(False)

def is_to_time_descendant(name_from ,name_to, tree):
    n_from = tree.search_nodes(name = name_from)[0]
    n_to = tree.search_nodes(name = name_to)[0].up
    if n_from.get_distance(samp) >= n_to.get_distance(samp):
        return(False)
    elif n_from.get_distance(samp) <= n_to.get_distance(samp):
        return(True)
    else:
        return("ERROR")

def get_ghost_prediction(com_node, back_bone_nodes, bak_bone_tree):
    # com_node = com.search_nodes(name = node.name)[0]
    while com_node.name not in back_bone_nodes:
        com_node = com_node.up
    ext_node = bak_bone_tree.search_nodes(name = com_node.name)[0]
    while len(ext_node.get_children()) not in [0,2]:
        ext_node = ext_node.get_children()[0]
    return(ext_node.name)

def get_ghost_banch_length(node, nodes_contemporary, com):
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


# loop for sample start here (and gene ?)
com = tr("T/CompleteTree.nwk",format = 1)

dirs = glob.glob("Sampl*")
for dir in dirs:
    samp_path = os.path.join(dir, "SampledSpeciesTree.nwk")
    samp = tr(samp_path,format = 1)
    back_bone_nodes = []
    for i in samp:
        node = com.search_nodes(name = i.name)[0]
        while node:
            back_bone_nodes.append(node.name)
            node = node.up
    back_bone_nodes = list(set(back_bone_nodes))
    bak_bone_tree = com.copy()
    bak_bone_tree.prune(list(set(back_bone_nodes)))
    genes = glob.glob(dir + "/*_sampledtree.nwk")
    all_prediction = []
    for gene in genes:
        try:
            all_prediction.append(complete_prediction(gene, dir, back_bone_nodes, samp))
        except:
            print(os.getcwd() + "/" + gene)
    try:
        full_prediction = pd.concat(all_prediction, ignore_index=True)
        full_prediction["to_direct_descendant"] = full_prediction.apply(lambda x: is_to_strict_descendant(x["from_prediction"], x["to_prediction"], samp), axis = 1)
        full_prediction["to_descendant"] = full_prediction.apply(lambda x: is_to_time_descendant(x["from_prediction"], x["to_prediction"], samp), axis = 1)
    except:
        full_prediction = pd.DataFrame()
    samp.dist = 0
    res_df = pd.DataFrame()
    res_df["Node"] = [i.name for i in samp.traverse()]
    res_df["br_length"] = res_df.apply(lambda x: samp.search_nodes(name = x["Node"])[0].dist, axis = 1)
    res_df["dist_to_root"] = res_df.apply(lambda x: round(samp.search_nodes(name = x["Node"])[0].get_distance(samp),6), axis = 1)
    if full_prediction.empty:
        res_df["N_transfers_donor"] = 0
        res_df["N_transfers_recip"] = 0
        res_df["to_direct_descendant"] = 0
        res_df["to_descendant"] = 0
    else:
        res_df["N_transfers_donor"] = res_df.apply(lambda x: full_prediction.loc[full_prediction["from_prediction"] == x["Node"]].shape[0], axis = 1)
        res_df["N_transfers_recip"] = res_df.apply(lambda x: full_prediction.loc[full_prediction["to_prediction"] == x["Node"]].shape[0], axis = 1)
        res_df["to_direct_descendant"] = res_df.apply(lambda x: full_prediction.loc[(full_prediction.from_prediction == x["Node"]) & (full_prediction.to_direct_descendant)].shape[0], axis = 1)
        res_df["to_descendant"] = res_df.apply(lambda x: full_prediction.loc[(full_prediction.from_prediction == x["Node"]) & (full_prediction.to_descendant)].shape[0], axis = 1)
    anc = com.search_nodes(name=samp.name)[0]
    extroot_to_comroot = anc.get_distance(com)
    nodes_contemporary = [node.name for node in anc.iter_descendants()]
    stat_by_node = dict()
    for i in com.traverse():
        prediction = get_ghost_prediction(i, back_bone_nodes, bak_bone_tree)
        if not i.is_root():
            i.add_features(ghost_prediction = prediction, ghost_dist = get_ghost_banch_length(i, nodes_contemporary, com))
        else:
            i.add_features(ghost_prediction = prediction, ghost_dist = 0)
        if prediction in stat_by_node:
            if i.name not in bak_bone_tree:
                stat_by_node[prediction][0] += 1
                stat_by_node[prediction][1] += i.ghost_dist
        else:
            stat_by_node[prediction] = [0, 0]
    res_df["N_ghost_node"] = res_df.apply(lambda x: stat_by_node[x["Node"]][0], axis = 1)
    res_df["L_ghost_branch"] = res_df.apply(lambda x: round(stat_by_node[x["Node"]][1],6), axis = 1)
    res_df.to_csv("Stats_simulation_" + dir + "_V2.txt", sep=' ', index = False)
    full_prediction.to_csv("Trans" + dir + "_V2.txt", sep=' ', index = False)




# GNU Ghost
