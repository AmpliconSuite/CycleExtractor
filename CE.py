from gurobipy import Model, GRB
import re
import networkx as nx
import random
from typing import Generator, List, Tuple, Optional
import numpy as np
from typing import cast
from collections import Counter
import time
from pathlib import Path
import sys
import pickle
from datetime import datetime
import os
import argparse
from pathlib import Path
##############################
def parse_graph_file(graph_file_path):
    segments = {}
    Path_Constraints = []
    seg_id = 1

    with open(graph_file_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            # --- Parse sequence (segment) lines ---
            if line.startswith("sequence"):
                parts = line.split()
                chrom_start = parts[1]  # e.g. chr12:75430274-
                chrom_end = parts[2]    # e.g. chr12:75531009+
                chrom = chrom_start.split(":")[0]
                start = int(chrom_start.split(":")[1].rstrip("+-"))
                end   = int(chrom_end.split(":")[1].rstrip("+-"))

                segments[seg_id] = {
                    "id": seg_id,
                    "chrom": chrom,
                    "start": start,
                    "end": end
                }
                seg_id += 1

            # --- Parse path constraint lines ---
            elif line.startswith("path_constraint"):
                # Example:
                # path_constraint e6-:1,c5+:1,e5-:1,d1+:1,e3-:1,c2+:1,e2-:1,d6+:1,e8-:1 66
                parts = re.split(r'\s+', line.strip(), maxsplit=2)
                if len(parts) < 3:
                    continue
                constraint_str = parts[1].strip()
                support_value = parts[2].strip()

                # Convert support to int if possible
                try:
                    support = int(float(support_value))
                except ValueError:
                    support = None

                # All "satisfied" set to True by default
                Path_Constraints.append({
                    "index": len(Path_Constraints) + 1,
                    "constraint_string": constraint_str,
                    "support": support,
                    "satisfied": True
                })

    return segments, Path_Constraints
############# Multiplicity ###################
class DiscordantEdge:
    def __init__(self, read_count: int):
        self.lr_count = read_count

class EdgeMultiplicityModel:
    def __init__(self, edge_dict: dict[int, dict]):
        self.discordant_edges = [
            DiscordantEdge(info["ReadCount"]) for info in edge_dict.values()
        ]

    @staticmethod
    def check_valid_discordant_rc_partition(
        rc_list: List[int], partition: List[int], max_multiplicity: int = 5
    ) -> Optional[Tuple[int, float]]:
        """Verify if a given partition of discordant edges meets the CoRAL criteria."""
        if partition[0] == partition[1]:
            return (partition[0], 0.0)
        rc_list_p = rc_list[partition[0] : partition[1] + 1]
        if rc_list_p[-1] < rc_list_p[0] * 2.0:
            return (partition[1], 0.0)
        base_ri = 0
        while base_ri < len(rc_list_p) and rc_list_p[base_ri] < rc_list_p[0] * 2.0:
            base_ri += 1
        base_avg_rc = cast(float, np.average(rc_list_p[:base_ri]))
        if rc_list_p[-1] / base_avg_rc >= max_multiplicity + 0.5:
            return None
        score = -10.0
        best_ri = base_ri
        sum_deviation = 1.0
        for base_ri_ in range(base_ri, 0, -1):
            base_avg_rc = cast(float, np.average(rc_list_p[:base_ri_]))
            base_size = len(rc_list_p[:base_ri_])
            sizes = {}
            li = base_ri_
            multiplicity = 2
            if rc_list_p[base_ri_] / base_avg_rc < multiplicity - 0.5:
                continue
            while rc_list_p[base_ri_] / base_avg_rc >= multiplicity + 0.5:
                multiplicity += 1
            sum_gap = np.log2(rc_list_p[base_ri_]) - np.log2(rc_list_p[base_ri_ - 1])
            for i in range(base_ri_, len(rc_list_p)):
                if rc_list_p[i] / base_avg_rc >= multiplicity + 0.5:
                    sum_gap += np.log2(rc_list_p[i]) - np.log2(rc_list_p[i - 1])
                    sizes[multiplicity] = [li, i - 1]
                    li = i
                    while rc_list_p[i] / base_avg_rc >= multiplicity + 0.5:
                        multiplicity += 1
            sizes[multiplicity] = [li, len(rc_list_p) - 1]
            if multiplicity > max_multiplicity:
                continue
            size_flag = True
            for m in range(2, multiplicity + 1):
                if m in sizes and sizes[m][1] - sizes[m][0] >= base_size:
                    size_flag = False
                    break
            if not size_flag:
                continue
            sum_deviation_ = sum(
                [
                    np.abs(
                        m
                        - np.average(
                            rc_list_p[sizes[m][0] : sizes[m][1] + 1] / base_avg_rc  # type: ignore[operator]
                        )
                    )
                    for m in range(2, multiplicity + 1)
                    if m in sizes
                ],
                0,
            )
            if sum_gap - sum_deviation_ > score:
                score = sum_gap - sum_deviation_
                sum_deviation = sum_deviation_
                best_ri = base_ri_
        if sum_deviation < 1.0:
            return (best_ri + partition[0] - 1, score)
        return None

    @staticmethod
    def enumerate_partitions(
        k: int, start: int, end: int
    ) -> Generator[List[List[int]], None, None]:
        """Generate all partitions of the interval [start, end] into k parts."""
        if k == 0:
            yield [[start, end]]
        else:
            for i in range(1, end - start - k + 2):
                for res in EdgeMultiplicityModel.enumerate_partitions(k - 1, start + i, end):
                    yield [[start, start + i - 1], *res]

    def infer_discordant_edge_multiplicities(
        self, max_multiplicity: int = 5
    ) -> List[int]:
        """Estimate the upper bound of multiplicities for each discordant edge."""
        rc_list = [de.lr_count for de in self.discordant_edges]  # Read counts
        if len(rc_list) == 0:
            return []
        rc_indices = np.argsort(rc_list)
        rc_list = sorted(rc_list)
        if np.log2(rc_list[-1]) - np.log2(rc_list[0]) < 1.0:
            return [1 for i in rc_indices]
        num_clusters = 1
        valid_clustering = False
        best_score_all = -10.0
        best_partitions = []
        distinct_all = []
    
        while not valid_clustering:
            valid_clustering = False
            for partitions in EdgeMultiplicityModel.enumerate_partitions(
                num_clusters - 1, 0, len(rc_list) - 1
            ):
                valid_partition = True
                score_all = 0.0
                distinct = []
                for pi in range(len(partitions)):
                    partition = partitions[pi]
                    scored_partition = EdgeMultiplicityModel.check_valid_discordant_rc_partition(
                        rc_list,
                        partition,
                        max_multiplicity,
                    )
                    if scored_partition is None:
                        valid_partition = False
                        break
                    base_ri, score = scored_partition
                    score_all += score
                    distinct.append([partitions[pi][0], base_ri])
                    if pi > 0:
                        score_all += np.log2(rc_list[partitions[pi][0]]) - np.log2(
                            rc_list[partitions[pi - 1][1]]
                        )
                if valid_partition:
                    valid_clustering = True
                    if score_all > best_score_all:
                        best_score_all = score_all
                        best_partitions = partitions
                        distinct_all = distinct
            if not valid_clustering:
                num_clusters += 1
        multiplicities_sorted = []
        for pi in range(len(best_partitions)):
            partition = best_partitions[pi]
            base_ = distinct_all[pi]
            for i in range(base_[0], base_[1] + 1):
                multiplicities_sorted.append(1)
            base_ri = base_[1] + 1
            if base_ri > partition[1]:
                continue
            base_avg_rc = np.average(rc_list[base_[0] : base_[1] + 1])
            multiplicity = 2
            while rc_list[base_ri] / base_avg_rc >= multiplicity + 0.5:
                multiplicity += 1
            for i in range(base_ri, partition[1] + 1):
                while rc_list[i] / base_avg_rc >= multiplicity + 0.5:
                    multiplicity += 1
                multiplicities_sorted.append(multiplicity)
        return [
            multiplicities_sorted[list(rc_indices).index(i)]
            for i in range(len(rc_list))
        ]
#################
def get_incident_edges(G, node, edge_type_filter=None):
    """
    Returns a list of (neighbor, edge_key, edge_data) for edges incident to `node`.
    Optionally filter by edge_type: 'sequence', 'concordant', 'discordant', or None for all.

    G: networkx.MultiGraph
    node: node label (str)
    edge_type_filter: str or None

    Returns:
        List of tuples: (neighbor_node, edge_key, edge_data_dict)
    """
    incident_edges = []
    for neighbor in G.neighbors(node):
        edges_between = G.get_edge_data(node, neighbor)
        if edges_between is None:
            continue
        for edge_key, edge_data in edges_between.items():
            if edge_type_filter is None or edge_data.get('type') == edge_type_filter:
                incident_edges.append((neighbor, edge_key, edge_data))
    return incident_edges 
#### Create the graph file with integer node for foldbacks
#### There iare no "s" and "t" here
def Create_graph(graph_file_path):
    G = nx.MultiGraph()
    i = 1  # Reset foldback ID counter

    path_constraints = []
    
    seq_edges_pc = []
    conc_edges_pc = []
    disc_edges_pc = []

    seq_edges_dict = dict()
    conc_edges_dict = dict()
    disc_edges_dict = dict()
    
    seq_counter = 0
    conc_counter = 0
    disc_counter = 0
    with open(graph_file_path, "r") as file:
        section = None
        
        for line in file:
            line = line.strip()
            
            # Identify sections
            if line.startswith("SequenceEdge:"):
                section = "sequence"
                continue
            elif line.startswith("BreakpointEdge:"):
                section = "breakpoint"
                continue
            elif line.startswith("PathConstraint:"):
                section = "path_constraint"
                continue
            
            
            
            # Skip empty lines
            if not line:
                continue
            
            # Extract nodes and sequence edges
            if section == "sequence" and line.startswith("sequence"):
                #parts = line.split("\t")
                parts = line.strip().split()
                _, start, end, copy_count, avg_coverage, size, num_reads = parts
                
                copy_count = float(copy_count)
                size = int(size)
                  
                G.add_nodes_from([start,end])
                G.add_edge(start, end, type='sequence', capacity=copy_count, length=size)
    
                seq_edges_dict[seq_counter] = {
                    "StartPosition": start,
                    "EndPosition": end,
                    "PredictedCN": copy_count,
                    "AverageCoverage": float(avg_coverage),
                    "Size": size,
                    "NumberOfLongReads": int(num_reads)
                }
                seq_counter += 1
            # Extract concordant and discordant edges
            if section == "breakpoint" and line.startswith("concordant"):
                parts = line.split("\t")
                start_end = parts[1]  # Start and End positions
                
                # Extract start and end nodes of breakpoint edges
                start, end = start_end.split("->")
                
            
                # if ref folder
                # copy_count = parts[3] if len(parts) > 3 else "0.0"  # Default to "0.0" if missing
                # if reconstructions folder
                copy_count = parts[2] if len(parts) > 3 else "0.0"  # Default to "0.0" if missing
               
                try:
                    copy_count = float(copy_count)  # Convert to float
                except ValueError:
                    copy_count = 0.0  # If conversion fails, set to 0.0
                
                
                G.add_nodes_from([start,end])
                G.add_edge(start, end, type='concordant', capacity=copy_count)
                
                conc_edges_dict[conc_counter] = {
                    "StartPosition": start,
                    "EndPosition": end,
                    "PredictedCN": copy_count
                }
                conc_counter += 1
                
            if section == "breakpoint" and line.startswith("discordant"):
                parts = line.split("\t")
                start_end = parts[1]  # Start and End positions
                
                # Extract start and end nodes of breakpoint edges
                start, end = start_end.split("->")
                
                # if ref folder
                # copy_count = parts[3] if len(parts) > 3 else "0.0"  # Default to "0.0" if missing
                # if reconstructions folder
                copy_count = parts[2] if len(parts) > 3 else "0.0"  # Default to "0.0" if missing
                read_count = parts[3] if len(parts) > 3 else "0.0"  # Default to "0.0" if missing
                try:
                    copy_count = float(copy_count)  # Convert to float
                except ValueError:
                    copy_count = 0.0  # If conversion fails, set to 0.0
                if start != end: # I want to avoid loops for now
                   G.add_nodes_from([start,end])
                   G.add_edge(start, end, type='discordant', capacity=copy_count, Read_Count=read_count)
                   
                   disc_edges_dict[disc_counter] = {
                        "StartPosition": start,
                        "EndPosition": end,
                        "PredictedCN": copy_count,
                        "ReadCount": read_count
                   }
                   disc_counter += 1
                else:
                   G.add_nodes_from([start, str(i), str(i + 1)])
                   G.add_edge(start, str(i), type='discordant', capacity=copy_count, Read_Count=read_count)
                   disc_edges_dict[disc_counter] = {
                        "StartPosition": start,
                        "EndPosition": str(i),
                        "PredictedCN": copy_count,
                        "ReadCount": read_count
                    }
                   disc_counter += 1
                   
                   G.add_edge(str(i),str(i+1), type='sequence', capacity=copy_count, length=1)
                   G.add_edge(str(i+1),start, type='discordant', capacity=copy_count, Read_Count=read_count)
                   disc_edges_dict[disc_counter] = {
                        "StartPosition": str(i + 1),
                        "EndPosition": start,
                        "PredictedCN": copy_count,
                        "ReadCount": read_count
                    }
                   disc_counter += 1
                   
                   i=i+2
                   
            if section == "path_constraint" and line.startswith("path_constraint"):
                    
                    parts = line.split("\t")
                    if len(parts) != 3:
                        continue  # skip malformed lines
        
                    _, path_str, support = parts
                    support = int(support)
                    # Extract the edge tokens, e.g., ['e4-:1', 'c3+:1', ...]
                    tokens = path_str.replace("path_constraint", "").strip()
                    edge_tokens = tokens.split(',')
                    
                    seq_set = []
                    conc_set = []
                    disc_set = []
            
        
                    for token in edge_tokens:
                        if ":" not in token:
                            continue
                        edge_token, _ = token.split(":")
                        edge_type = edge_token[0]
                        edge_num = int(edge_token[1:-1])
                    
                        if edge_type == "e" and edge_num-1 in seq_edges_dict:
                            seq_set.append(seq_edges_dict[edge_num-1])
                        elif edge_type == "c" and edge_num-1 in conc_edges_dict:
                            conc_set.append(conc_edges_dict[edge_num-1])
                        elif edge_type == "d" and edge_num-1 in disc_edges_dict:
                            disc_set.append(disc_edges_dict[edge_num-1])
                            
                    seq_edges_pc.append(seq_set)
                    conc_edges_pc.append(conc_set)
                    disc_edges_pc.append(disc_set)
                    
                    path_constraints.append({
                            "paths": edge_tokens,
                            "support": support
                        })
                        
                    
                
        disc_edges_readcount= {
            k: {"ReadCount": float(v["ReadCount"])}
            for k, v in disc_edges_dict.items()
        }
        
        # Pass it to the model class
        Multiplicity_bound_model = EdgeMultiplicityModel(disc_edges_readcount)

        # Call the function to get multiplicities
        K = Multiplicity_bound_model.infer_discordant_edge_multiplicities(max_multiplicity=5)
        
        for edge_key, k in zip(disc_edges_dict.keys(), K):
            disc_edges_dict[edge_key]["K"] = k
            
            
        for edge_info in disc_edges_dict.values():
            u = edge_info["StartPosition"]
            v = edge_info["EndPosition"]
            k_val = edge_info["K"]
        
            for key in G[u][v]:
                attr = G[u][v][key]
                if attr.get("type") == "discordant":
                    G[u][v][key]["K"] = k_val
                    
        nodes=G.nodes
        
        concordant_edges = {tuple(sorted((u, v))) for u, v, attr in G.edges(data=True) if attr.get('type') == 'concordant'}
        capacity_concordant_edges = {tuple(sorted((u,v))): attr['capacity'] for u, v, attr in G.edges(data=True) if attr.get('type') == 'concordant'}
        
        discordant_edges = {tuple(sorted((u, v))) for u, v, attr in G.edges(data=True) if attr.get('type') == 'discordant'}
        capacity_discordant_edges = {tuple(sorted((u,v))): attr['capacity'] for u, v, attr in G.edges(data=True) if attr.get('type') == 'discordant'}
        read_count_discordant_edges = {tuple(sorted((u,v))): attr['Read_Count'] for u, v, attr in G.edges(data=True) if attr.get('type') == 'discordant'}
        K_discordant_edges = {tuple(sorted((u, v))): attr["K"]for u, v, key, attr in G.edges(keys=True, data=True)if attr.get("type") == "discordant" and "K" in attr}
        
        sequence_edges = {tuple(sorted((u,v))) for u, v, attr in G.edges(data=True) if attr.get('type') == 'sequence'}
        capacity_sequence_edges = {tuple(sorted((u, v))): attr['capacity'] for u, v, attr in G.edges(data=True) if attr.get('type') == 'sequence'}
        length_sequence_edges = {tuple(sorted((u, v))): attr['length'] for u, v, attr in G.edges(data=True) if attr.get('type') == 'sequence'}
        
        
        #### define p_ij^k which indictae if edge ij is in P_k or not
        #### it is ones the start and the end of the edges nd it shows the number of the pc that edge belongs to -1
        ##### ("chr3:76247852, chr3:98628764,0)=1" 0 means pc 1
        

        # These will be dictionaries of (i, j, k) → 1
        p_concordant_edges = {}
        p_discordant_edges = {}
        p_sequence_edges = {}
        
        for k, path_constraint in enumerate(path_constraints):
            edge_tokens = path_constraint["paths"]
        
            for token in edge_tokens:
                if ":" not in token:
                    continue
                edge_token, _ = token.split(":")
        
                match = re.match(r"([ecd])(\d+)", edge_token)
                if not match:
                    continue
        
                edge_type, edge_num = match.groups()
                edge_num = int(edge_num)
        
                if edge_type == "e" and (edge_num - 1) in seq_edges_dict:
                    edge_data = seq_edges_dict[edge_num - 1]
                    i, j = sorted((edge_data["StartPosition"], edge_data["EndPosition"]))
                    p_sequence_edges[(i, j, k)] = 1
        
                elif edge_type == "c" and (edge_num - 1) in conc_edges_dict:
                    edge_data = conc_edges_dict[edge_num - 1]
                    i, j = sorted((edge_data["StartPosition"], edge_data["EndPosition"]))
                    p_concordant_edges[(i, j, k)] = 1
        
                elif edge_type == "d" and (edge_num - 1) in disc_edges_dict:
                    edge_data = disc_edges_dict[edge_num - 1]
                    i, j = sorted((edge_data["StartPosition"], edge_data["EndPosition"]))
                    p_discordant_edges[(i, j, k)] = 1
        
    return (
            G,
            nodes,
            concordant_edges,
            capacity_concordant_edges,
            discordant_edges,
            capacity_discordant_edges,
            read_count_discordant_edges,
            K_discordant_edges,
            sequence_edges,
            capacity_sequence_edges,
            length_sequence_edges,
            p_sequence_edges,
            p_concordant_edges,
            p_discordant_edges,
            path_constraints
        )

############ creat graph with nodes "s" and "t" for one connected cycle extraction
def creat_s_t_graph_for_connected_cycles(G,capacity_discordant_edges,capacity_sequence_edges,K_discordant_edges,length_sequence_edges):
    G.add_node('s')
    G.add_node('t')
    
    new_capacity_discordant_edges = dict(capacity_discordant_edges)
    new_capacity_sequence_edges = dict(capacity_sequence_edges)
    
    for node in G.nodes():
        if node in {'s', 't'} or (isinstance(node, str) and node.isdigit()):
            continue
        
        # Compute capacity from incident edges
        dc_cns = [data.get('capacity', 0) for _, _, data in get_incident_edges(G, node, edge_type_filter='discordant')]
        seq_cns = [data.get('capacity', 0) for _, _, data in get_incident_edges(G, node, edge_type_filter='sequence')]
        capacity = max(dc_cns + seq_cns) if (dc_cns + seq_cns) else 1
        
        # Add edges (undirected)
        G.add_edge('s', node, type='discordant', capacity=capacity, Read_Count=1000, K=1)
        G.add_edge('t', node, type='sequence', capacity=capacity, Read_Count=1000, length=1)
        
        # Update capacities using sorted tuple
        new_capacity_discordant_edges[tuple(sorted(('s', node)))] = capacity
        new_capacity_sequence_edges[tuple(sorted(('t', node)))] = capacity
    
    discordant_edges = {tuple(sorted((u, v))) for u, v, attr in G.edges(data=True) if attr.get('type') == 'discordant'}
    sequence_edges = {tuple(sorted((u, v))) for u, v, attr in G.edges(data=True) if attr.get('type') == 'sequence'}
    
    K_discordant_edges = {tuple(sorted((u, v))): attr["K"] for u, v, key, attr in G.edges(keys=True, data=True) if attr.get("type") == "discordant" and "K" in attr}
    length_sequence_edges = {tuple(sorted((u, v))): attr.get("length", 1) for u, v, key, attr in G.edges(keys=True, data=True) if attr.get("type") == "sequence"}
    
    return G, G.nodes, discordant_edges, new_capacity_discordant_edges, sequence_edges, new_capacity_sequence_edges, K_discordant_edges, length_sequence_edges

########## Calculate Length weighted copy number of the graph
def compute_Length_weighted_copy_number_graph(capacity_sequence_edges,length_sequence_edges):
    valid_edges = [
        tuple(sorted((i, j)))
        for (i, j) in sequence_edges
        if i != "t" and j != "t" and i != "s" and j != "s" and not i.isdigit() and not j.isdigit()
    ]
    weight = sum(
        capacity_sequence_edges[edge] * length_sequence_edges[edge]
        for edge in valid_edges
    )
    return weight
########## Calculate Length weighted copy number of the graph. Here we do not consider segments with CN<5. (CN of the graph)
def compute_LWCN_graph_excluding_some_edges(capacity_sequence_edges,length_sequence_edges):
    valid_edges = [
        tuple(sorted((i, j)))
        for (i, j) in sequence_edges
        if i != "t" and j != "t" 
           and i != "s" and j != "s" and not i.isdigit() and not j.isdigit()
           and capacity_sequence_edges[tuple(sorted((i, j)))]>5
        
    ]
    weight = sum(
        capacity_sequence_edges[edge] * length_sequence_edges[edge]
        for edge in valid_edges
    )
    return weight
######### Build the model
def build_cycle_model(
                G,
                nodes,
                concordant_edges,
                capacity_concordant_edges,
                discordant_edges,
                capacity_discordant_edges,
                read_count_discordant_edges,
                K_discordant_edges,
                sequence_edges,
                capacity_sequence_edges,
                length_sequence_edges,
                p_sequence_edges,
                p_concordant_edges,
                p_discordant_edges,
                path_constraints,
                enforce_connectivity,
                gamma):
    model = Model("ILP_CoRAL_cycle")
    # Decision variables   
    X_concordant_edge = model.addVars(concordant_edges, vtype=GRB.BINARY, name="X_concordant_edge")
    X_discordant_edge = model.addVars(discordant_edges, vtype=GRB.BINARY, name="X_discordant_edge")
    X_sequence_edge = model.addVars(sequence_edges, vtype=GRB.BINARY, name="X_sequence_edge")   
    f_concordant_edge = model.addVars(concordant_edges, vtype=GRB.CONTINUOUS, lb=0, name="f_concordant_edge")
    f_discordant_edge = model.addVars(discordant_edges, vtype=GRB.CONTINUOUS, lb=0, name="f_discordant_edge")
    f_sequence_edge = model.addVars(sequence_edges, vtype=GRB.CONTINUOUS, lb=0, name="f_sequence_edge")    
    F = model.addVar(name="f_min", lb=0, vtype=GRB.CONTINUOUS)  
    P = model.addVars(len(path_constraints), vtype=GRB.BINARY, name="P_path_constraint")
    y_discordant = {}
    for i, j in discordant_edges:
        edge_key = tuple(sorted((i, j)))
        for m in range(1, K_discordant_edges[edge_key] + 1):
            y_discordant[edge_key + (m,)] = model.addVar(name=f"y_discordant_{edge_key[0]}_{edge_key[1]}_{m}", vtype=GRB.BINARY)
    z_discordant = {}
    for i, j in discordant_edges:
        edge_key = tuple(sorted((i, j)))
        for m in range(1, K_discordant_edges[edge_key] + 1):
            z_discordant[edge_key + (m,)] = model.addVar(lb=0, name=f"z_discordant_{edge_key[0]}_{edge_key[1]}_{m}", vtype=GRB.CONTINUOUS)

    
    ### Objective function
    if len(path_constraints)>0:
        multiplier= gamma* sum(capacity_sequence_edges[tuple(sorted((i, j)))] * length_sequence_edges[tuple(sorted((i, j)))] for (i, j) in sequence_edges if i != 't' and j != 't')/ len(path_constraints)
    else:
        multiplier=0
    
   ### Objective: Maximize total flow on sequence edges weighted by length and number of path constraints satisfied
    model.setObjective(
        sum(f_sequence_edge[tuple(sorted((i, j)))] * length_sequence_edges[tuple(sorted((i, j)))] for (i, j) in sequence_edges if i != 't' and j != 't')+ multiplier*sum(P[k] for k in range(len(path_constraints))),
        GRB.MAXIMIZE
    )
    # Constraint 1: Copy number of an edge ≤ capacity × selection variable (for undirected edges)
    
    for (i, j) in concordant_edges:
        edge = tuple(sorted((i, j)))
        model.addConstr(
            f_concordant_edge[edge] <= capacity_concordant_edges[edge] * X_concordant_edge[edge],
            name=f"ConcordantEdgeCapacity_{edge[0]}_{edge[1]}")    
    for (i, j) in discordant_edges: 
        edge = tuple(sorted((i, j)))
        model.addConstr(
            f_discordant_edge[edge] <= capacity_discordant_edges[edge] * X_discordant_edge[edge],
            name=f"DiscordantEdgeCapacity_{edge[0]}_{edge[1]}")    
    for (i, j) in sequence_edges: 
        edge = tuple(sorted((i, j)))
        model.addConstr(
            f_sequence_edge[edge] <= capacity_sequence_edges[edge] * X_sequence_edge[edge],
            name=f"SequenceEdgeCapacity_{edge[0]}_{edge[1]}")
    ### Constraint 3: copy number sequence edge incident = copy number discordant/concordant edge incident 
    ##########################################################################    
    for k in nodes:
        if k in {"s", "t"}:
            continue
    
        # Define undirected incident edges
        incident_edges = {
          tuple(sorted((i, k)))
          for i in nodes
          if tuple(sorted((i, k))) in concordant_edges | discordant_edges | sequence_edges and i != k}
        # Flow into k (sequence edges)
        flow_in_sequence = sum(
            f_sequence_edge[edge]
            for edge in incident_edges
            if edge in sequence_edges and k in edge)
    
        # Flow out of k (concordant + discordant)
        flow_out_concordant = sum(
            f_concordant_edge[edge]
            for edge in incident_edges
            if edge in concordant_edges and k in edge)
        flow_out_discordant = sum(
            f_discordant_edge[edge]
            for edge in incident_edges
            if edge in discordant_edges and k in edge)
    
        model.addConstr(
            flow_in_sequence == flow_out_concordant + flow_out_discordant,
            name=f"FlowConservation2_{k}")
    ### Constraint 7 #Big-M constant (ensure it's large enough)    
    M = 1000000    
    for edge in concordant_edges:
        i, j = edge
        edge_key = tuple(sorted((i, j)))
        model.addConstr(
            F <= f_concordant_edge[edge_key] + M * (1 - X_concordant_edge[edge_key]),
            name=f"ConcordantEdgeCapacity_{i}_{j}")    
    for edge in discordant_edges:
        i, j = edge
        edge_key = tuple(sorted((i, j)))
        model.addConstr(
            F <= f_discordant_edge[edge_key] + M * (1 - X_discordant_edge[edge_key]),
            name=f"DiscordantEdgeCapacity_{i}_{j}")
    for i, j in discordant_edges:
        edge_key = tuple(sorted((i, j)))
        model.addConstr(
            sum(y_discordant[edge_key + (m,)] for m in range(1, K_discordant_edges[edge_key] + 1)) == X_discordant_edge[edge_key],
            name=f"one_multiple_discordant_{edge_key[0]}_{edge_key[1]}")
    for i, j in discordant_edges:
        edge_key = tuple(sorted((i, j)))
        for m in range(1, K_discordant_edges[edge_key] + 1):
            z = z_discordant[edge_key + (m,)]
            y = y_discordant[edge_key + (m,)]
            model.addConstr(z <= F, name=f"z_leq_fmin_discordant_{edge_key[0]}_{edge_key[1]}_{m}")
            model.addConstr(z <= M * y, name=f"z_leq_My_discordant_{edge_key[0]}_{edge_key[1]}_{m}")
            model.addConstr(z >= F - M * (1 - y), name=f"z_geq_fmin_minus_M_discordant_{edge_key[0]}_{edge_key[1]}_{m}")
            model.addConstr(z >= 0, name=f"z_geq_0_discordant_{edge_key[0]}_{edge_key[1]}_{m}")
    for i, j in discordant_edges:
        edge_key = tuple(sorted((i, j)))
        model.addConstr(
            f_discordant_edge[edge_key] == sum(m * z_discordant[edge_key + (m,)] for m in range(1, K_discordant_edges[edge_key] + 1)),
            name=f"flow_def_discordant_{edge_key[0]}_{edge_key[1]}")
    #### Constarints 9 Path constraints
    for (i, j, k) in p_sequence_edges:
        edge_key = tuple(sorted((i, j)))
        model.addConstr(
            X_sequence_edge[edge_key] >= P[k] * p_sequence_edges[(i, j, k)])
    for (i, j, k) in p_concordant_edges:
        edge_key = tuple(sorted((i, j)))
        model.addConstr(
            X_concordant_edge[edge_key] >= P[k] * p_concordant_edges[(i, j, k)])
    for (i, j, k) in p_discordant_edges:
        edge_key = tuple(sorted((i, j)))
        model.addConstr(
            X_discordant_edge[edge_key] >= P[k] * p_discordant_edges[(i, j, k)])
    
    if enforce_connectivity:
        #####################################################    
        ##### Ensure Connectivity Variables
        #####################################################
        all_edges = concordant_edges | discordant_edges | sequence_edges
        # Generate all directed versions of directed variables on undirected edges
        directed_sequence_edges = {(u, v) for u, v in sequence_edges} | {(v, u) for u, v in sequence_edges}
        directed_concordant_edges = {(u, v) for u, v in concordant_edges} | {(v, u) for u, v in concordant_edges}
        directed_discordant_edges = {(u, v) for u, v in discordant_edges} | {(v, u) for u, v in discordant_edges}
        c_seq  = model.addVars(directed_sequence_edges,   vtype=GRB.BINARY, name="c_seq")
        c_conc = model.addVars(directed_concordant_edges, vtype=GRB.BINARY, name="c_conc")
        c_disc = model.addVars(directed_discordant_edges, vtype=GRB.BINARY, name="c_disc")
        d = model.addVars(G.nodes, vtype=GRB.INTEGER, lb=0, ub=len(G.nodes), name="d")
        
        ##### Connectivity constraints for s and t  Constraint  "s" and "t" edges and the copy number of their incident edges
        flow_out_s_discordant = sum(
                f_discordant_edge[tuple(sorted(("s", k)))]
                for k in nodes
                if tuple(sorted(("s", k))) in discordant_edges
            )
            
            ### Total flow into sink 't' on discordant edges
        flow_in_t_sequence = sum(
                f_sequence_edge[tuple(sorted((k, "t")))]
                for k in nodes
                if tuple(sorted((k, "t"))) in sequence_edges
            )
            
            ### Enforce the flow match
        model.addConstr(
                flow_out_s_discordant == flow_in_t_sequence,
                name="TotalFlowFromS_Equals_FlowIntoT"
            )
            
        model.addConstr(
                sum(
                    X_discordant_edge[tuple(sorted(("s", k)))]
                    for k in nodes
                    if tuple(sorted(("s", k))) in discordant_edges
                ) == 1,
                name="OneEdgeOutOfS"
            )
    
        model.addConstr(
                sum(
                    X_sequence_edge[tuple(sorted((k, "t")))]
                    for k in nodes
                    if tuple(sorted((k, "t"))) in sequence_edges
                ) == 1,
                name="OneEdgeIntoT"
            )
        #if node in {'s', 't'} or (isinstance(node, str) and node.isdigit()):
        for k in nodes:
            if k in {"s", "t"} or k.isdigit():
               continue
            model.addConstr(
                X_discordant_edge[tuple(sorted(("s", k)))] == X_sequence_edge[tuple(sorted((k, "t")))],
                name=f"s_t_same_node_{k}"
            )
        
        ######################################################
        ######### Connectivity constraint 0 #################
        ###### Just about 's' which is the start node #######
        ##### d(s)=1, y_plus ('s', only one node)=1 that is where X_...=1
        ##### y_minu('s', every other node)=0, y_plus(every other node,'s')=0
        ###################################################### 
       
        for u, v in directed_discordant_edges:
            if u != 's':
                continue
            edge = tuple(sorted((u, v)))  # Get the undirected version
        
            # Only add constraint if this edge exists in X variables
            
            if edge in X_discordant_edge:
                x_var = X_discordant_edge[edge]
            else:
                continue  # Skip if edge not in any X set
        
            model.addConstr(
                c_disc[u, v] ==x_var,
                name=f"y_dir_limit_{u}_{v}"
            )  
            model.addConstr(
                c_disc[v,u] ==0,
                name=f"y_dir_limit_{u}_{v}"
            ) 
            
            
        model.addConstr(d['s'] == 1, name='d_start_node_one')
        ######################################################
        ######### Connectivity constraint 1 #################
        ######################################################
        num_nodes = len(G.nodes)
    
        for u in G.nodes:
            incident_x_vars = []
            for v in G.neighbors(u):
                edge = tuple(sorted((u, v)))
                if edge in sequence_edges:
                    incident_x_vars.append(X_sequence_edge[edge])
                if edge in concordant_edges:
                    incident_x_vars.append(X_concordant_edge[edge])
                if edge in discordant_edges:
                    incident_x_vars.append(X_discordant_edge[edge])
            model.addConstr(
                d[u] <= num_nodes * sum(incident_x_vars),
                name=f"only_assign_d_if_used_{u}"
            )
            
        ###################################################
        ###### Connectivity constraint 2 ##################
        ###################################################
                            
        for u, v in sequence_edges:
            edge = tuple(sorted((u, v)))
            model.addConstr(c_seq[u, v] + c_seq[v, u] <= X_sequence_edge[(u, v)], name=f"dir_seq_use_{u}_{v}")
        
        for u, v in concordant_edges:
            edge = tuple(sorted((u, v)))
            model.addConstr(c_conc[u, v] + c_conc[v, u] <= X_concordant_edge[(u, v)], name=f"dir_conc_use_{u}_{v}")
        
        for u, v in discordant_edges:
            edge = tuple(sorted((u, v)))
            model.addConstr(c_disc[u, v] + c_disc[v, u] <= X_discordant_edge[(u, v)], name=f"dir_disc_use_{u}_{v}")
        
        ###################################################
        ######### Connectivity constraint 3 ###############
        ###################################################
        all_nodes = list(G.nodes)
        for u in all_nodes:
            if u == 's':
                continue  # skip source
            incoming_c = []
    
            for v in all_nodes:
                if (v, u) in c_seq:
                    incoming_c.append(c_seq[v, u])
                if (v, u) in c_conc:
                    incoming_c.append(c_conc[v, u])
                if (v, u) in c_disc:
                    incoming_c.append(c_disc[v, u])
        
            incoming_sum = sum(incoming_c)
        
            # Compute total x variables that involve u (undirected edge participation)
            x_incident = []
            for v in G.neighbors(u):
                edge = tuple(sorted((u, v)))  # undirected
                if edge in X_sequence_edge:
                    x_incident.append(X_sequence_edge[edge])
                if edge in X_concordant_edge:
                    x_incident.append(X_concordant_edge[edge])
                if edge in X_discordant_edge:
                    x_incident.append(X_discordant_edge[edge])
        
            x_total = sum(x_incident)
            print('x_total for', u, 'is', x_total)
        
            # Constraint: if node u is in the solution, it must have at least one incoming edge
            model.addConstr(
                incoming_sum >= (1/num_nodes)*x_total,
                name=f"incoming_constraint_{u}"
            )
        ###################################################
        ######### Connectivity constraint 4 ###############
        ###################################################
                        
        for u, v in directed_sequence_edges:
            if (u, v) in c_seq:
                model.addConstr(
                    d[v] - d[u] + num_nodes * (1 - c_seq[u, v]) >= 1,
                    name=f"d_order_sequence_{u}_{v}"
                )
        for u, v in directed_concordant_edges:
            if (u, v) in c_conc:
                model.addConstr(
                    d[v] - d[u] + num_nodes * (1 - c_conc[u, v]) >= 1,
                    name=f"d_order_concordant_{u}_{v}"
                )
        for u, v in directed_discordant_edges:
            if (u, v) in c_disc:
                model.addConstr(
                    d[v] - d[u] + num_nodes * (1 - c_disc[u, v]) >= 1,
                    name=f"d_order_discordant_{u}_{v}"
                )
        ###################################################
        ######### End of constraints ###############
        ###################################################
    return model, {
        "X_concordant_edge": X_concordant_edge,
        "X_discordant_edge": X_discordant_edge,
        "X_sequence_edge": X_sequence_edge,
        "f_concordant_edge": f_concordant_edge,
        "f_discordant_edge": f_discordant_edge,
        "f_sequence_edge": f_sequence_edge,
        "F": F,
        "P": P,
        "y_discordant": y_discordant,
        "z_discordant": z_discordant,
    }
########## Some small functions to use in the Heuristic later


def use_edge(counter, edge):
    counter[edge] -= 1
    if counter[edge] == 0:
        del counter[edge]
###### test if the path of walk is alternating betwenn sequence and concordant/discordant edge
def is_edge_type_alternating(walk):
    
    # Check that edge labels alternate between:
    # 'sequence edge' and ('discordant edge' or 'concordant edge')
    # Nodes are ignored.

    edge_sequence = [label for (label, _) in walk if "edge" in label]
    if not edge_sequence:
        return True  # Nothing to validate
    
    expected = "sequence edge" if edge_sequence[0] == "sequence edge" else "non-sequence"
    
    for label in edge_sequence:
        if expected == "sequence edge":
            if label != "sequence edge":
                return False
            expected = "non-sequence"
        else:
            if label not in {"discordant edge", "concordant edge"}:
                return False
            expected = "sequence edge"
    
    return True
######### Before merging, rotate the walk if it is necessary
def find_and_rotate_walk(walk1, walk2):
    
    # Try merging walk2 into walk1 at shared nodes,
    # testing both entry directions, and preserve edge-type alternation.

    nodes1 = {entry[1] for entry in walk1 if entry[0] == "node"}
    shared_indices = [(i, entry[1]) for i, entry in enumerate(walk2)
                      if entry[0] == "node" and entry[1] in nodes1]

    if not shared_indices:
        print("⚠️ Warning: No shared node found. Skipping merge.")
        return walk1

    for shared_index, shared_node in shared_indices:
       # for offset in (0, -2): # in daghighan roo shared node shoroonemikard vase hamin bazi vaghta javab ghalat midad. masalan 2 ta node posht sare ham midad
        for offset in (0,):    
            walk2_copy = walk2[:]

            # Handle redundant closure
            if walk2_copy[0][0] == "node" and walk2_copy[-1][0] == "node" and walk2_copy[0][1] == walk2_copy[-1][1]:
                walk2_copy = walk2_copy[:-1]

            rotated_start = shared_index + offset
            if not (0 <= rotated_start < len(walk2_copy)):
                continue

            rotated = walk2_copy[rotated_start:] + walk2_copy[:rotated_start]

            # Remove repeated start node if present
            if rotated and rotated[0] == ("node", shared_node):
                rotated = rotated[1:]
            rotated.append(("node", shared_node))  # Always close the loop

            # Insert rotated walk2 at shared_node in walk1
            insert_index = next(i for i, (lbl, node) in enumerate(walk1)
                                if lbl == "node" and node == shared_node)

            merged_candidate = walk1[:insert_index + 1] + rotated + walk1[insert_index + 1:]
            # print("First part of walk 1: ", walk1[:insert_index + 1])
            # print("walk 2 rotated: ", rotated)
            # print("Second part of walk 1: ", walk1[insert_index + 1:])
            if is_edge_type_alternating(merged_candidate):
                return merged_candidate
            # --- Try reversed rotated walk ---
            # Keep closure at start/end, reverse middle
            middle = rotated[:-1]
            reversed_middle = []
            # reverse the order of edges + nodes (node after edge becomes node before edge)
            for i in range(len(middle)-1, 0, -2):
                reversed_middle.append(middle[i])      # node
                reversed_middle.append(middle[i-1])    # edge
            # print("reversed middle walk2 beofre append: ", reversed_middle)
            reversed_middle.append(middle[0]) 
            # print("reversed middle walk2: ", reversed_middle)
            
            #reversed_rotated = [("node", shared_node)] + reversed_middle + [("node", shared_node)]
            reversed_rotated =  reversed_middle + [("node", shared_node)]

            # print("reversed walk 2 rotated: ", reversed_rotated)
            merged_candidate_reversed = walk1[:insert_index + 1] + reversed_rotated + walk1[insert_index + 1:]
            if is_edge_type_alternating(merged_candidate_reversed):
                return merged_candidate_reversed

    print("⚠️ Warning: No valid merge preserving alternation. Skipping merge.")
    return walk1
########## Heuristic to create cycles
def create_closed_walks(sequence_edges, f_sequence_edge, discordant_edges, f_discordant_edge, 
                     concordant_edges, f_concordant_edge, max_attempts,segments, F):
    
    node_to_segment_id = {}
    for seg_id, info in segments.items():
        node_start = f"{info['chrom']}:{info['start']}-"
        node_end = f"{info['chrom']}:{info['end']}+"
        node_to_segment_id[node_start] = seg_id
        node_to_segment_id[node_end] = seg_id
    def edge_segment_rank(edge):
        return min(
            node_to_segment_id.get(edge[0], float('inf')),
            node_to_segment_id.get(edge[1], float('inf'))
        )
    # Filter and normalize undirected sequence edges
    solution_sequence_edges = []
    for (i, j) in sequence_edges:
        if i == 't' or j == 't':
           continue  # This happens only if we enforce connectivity, I find this if check eaiser than the enforce connectivity flag 
        if f_sequence_edge[i, j] > 0:
            multiplicity = round(f_sequence_edge[i, j] / F)
            solution_sequence_edges.extend([tuple(sorted((i, j)))] * multiplicity)
    # print(" Sequence edges in the solution considering multiplicity: ")
    # for edge in solution_sequence_edges:
    #     print(edge)
        
    solution_discordant_edges = []
    for (i, j) in discordant_edges:
        if i == 's' or j == 's':
           continue  # This happens only if we enforce connectivity, I find this if check eaiser than the enforce connectivity flag 
        if f_discordant_edge[i, j] > 0:
            multiplicity = round(f_discordant_edge[i, j] / F)
            solution_discordant_edges.extend([tuple(sorted((i, j)))] * multiplicity)
    # print(" Discordant edges in the solution considering multiplicity: ")
    # for edge in solution_discordant_edges:
    #     print(edge)
        
    solution_concordant_edges = []
    for (i, j) in concordant_edges:
        if f_concordant_edge[i, j] > 0:
            multiplicity = round(f_concordant_edge[i, j] / F)
            solution_concordant_edges.extend([tuple(sorted((i, j)))] * multiplicity)
    # print(" Concordant edges in the solution considering multiplicity: ")
    # for edge in solution_concordant_edges:
    #     print(edge)       
    
    attempt = 0
    error_in_closed_walks = 1  # default to failure
    
    while attempt < max_attempts and error_in_closed_walks == 1:
        attempt += 1
        # print(f"Attempt {attempt} to construct valid closed walks...")
    
        # Make fresh copies of the edge sets
        remaining_sequence_edges   = Counter(solution_sequence_edges)
        remaining_discordant_edges = Counter(solution_discordant_edges)
        remaining_concordant_edges = Counter(solution_concordant_edges)
    
        all_closed_walks = []
        #first_walk = True
    
        while remaining_sequence_edges:
            # if first_walk:
            #     sorted_candidates = sorted(remaining_sequence_edges, key=edge_segment_rank)
    
            #     for candidate in sorted_candidates:
            #         n1, n2 = candidate
            #         bad = False
            #         for d_edge in remaining_discordant_edges:
            #             if (n1 in d_edge and any(x.isdigit() and not x.startswith("chr") for x in d_edge if x != n1)) or \
            #                (n2 in d_edge and any(x.isdigit() and not x.startswith("chr") for x in d_edge if x != n2)):
            #                 bad = True
            #                 break
            #         if not bad:
            #             start_edge = candidate
            #             break
            #     else:
            #         start_edge = sorted_candidates[0]
    
            #     first_walk = False
            # else:
            start_edge = random.choice(list(remaining_sequence_edges.elements()))
    
            ordered_closed_walk = [
                ("node", start_edge[0]),
                ("sequence edge", start_edge),
                ("node", start_edge[1])
            ]
            use_edge(remaining_sequence_edges, start_edge)
    
            current_node = start_edge[1]
            edge_type = "non-sequence"
    
            while True:
                if edge_type == "non-sequence":
                    # match = None
                    # for edge in sorted(remaining_discordant_edges, key=edge_segment_rank):
                    #     if current_node in edge:
                    #         match = ("discordant edge", edge)
                    #         use_edge(remaining_discordant_edges, edge)
                    #         break
                    # if match is None:
                    #     for edge in sorted(remaining_concordant_edges, key=edge_segment_rank):
                    #         if current_node in edge:
                    #             match = ("concordant edge", edge)
                    #             use_edge(remaining_concordant_edges, edge)
                    #             break
                    # if match is None:
                    #     break
                
                    match = None
                    # Combine discordant and concordant edges into one list
                    candidate_edges = [
                        ("discordant edge", edge) for edge in remaining_discordant_edges.elements()
                    ] + [
                        ("concordant edge", edge) for edge in remaining_concordant_edges.elements()
                    ]
                    
                    
                    #### Pick randomly 
                    candidate_edges = [x for x in candidate_edges if current_node in x[1]]
                    
                    if candidate_edges:
                        match = random.choice(candidate_edges)
                        edge_label, edge = match
                        if edge_label == "discordant edge":
                            use_edge(remaining_discordant_edges, edge)
                        else:
                            use_edge(remaining_concordant_edges, edge)
                    else:
                        break 
                    
                    if match is None:
                        break
    
                    edge_label, edge = match
                    next_node = edge[1] if edge[0] == current_node else edge[0]
                    ordered_closed_walk.append((edge_label, edge))
                    ordered_closed_walk.append(("node", next_node))
                    current_node = next_node
                    edge_type = "sequence"
    
                elif edge_type == "sequence":
                    match = None
                    
                    candidate_edges = [edge for edge in remaining_sequence_edges.elements() if current_node in edge]
                    
                    if candidate_edges:
                        match = random.choice(candidate_edges)  # pick one randomly
                        use_edge(remaining_sequence_edges, match)
                    else:
                        break
                     
                    if match is None:
                        break
    
                    next_node = match[1] if match[0] == current_node else match[0]
                    ordered_closed_walk.append(("sequence edge", match))
                    ordered_closed_walk.append(("node", next_node))
                    current_node = next_node
                    edge_type = "non-sequence"
    
                if current_node == ordered_closed_walk[0][1]:
                    break
    
            all_closed_walks.append(ordered_closed_walk)
            
        # --- Evaluate correctness after each attempt ---
        all_walks_closed = all(
            walk[0][1] == walk[-1][1]
            for walk in all_closed_walks
        )
    
        all_edges_used = (
            not remaining_sequence_edges and
            not remaining_concordant_edges and
            not remaining_discordant_edges
        )
        
        # if not all_edges_used:
        #     print("Some edges are not used in the final closed walks:" )
        #     print("These are the remaining sequence edges:", remaining_sequence_edges )
        #     print("These are the remaining concordant edges:", remaining_concordant_edges )
        #     print("These are the remaining discordant edges:", remaining_discordant_edges )
            
        # if not all_walks_closed:
        #     print(f"There is a walk with CN= {F} that is not closed here:" , all_closed_walks)

    
        error_in_closed_walks = 0 if all_walks_closed and all_edges_used else 1
    return all_closed_walks, error_in_closed_walks

######### merge the closed walks
def merge_all_closed_walks(walks):
    if not walks:
        return [], []   # Always return a tuple
    
    if len(walks) == 1:
        return walks[0], []   # Always return a tuple

    merged_walk = walks[0]
    not_merged_walks = []

    for next_walk in walks[1:]:
        old_merged = merged_walk
        merged_candidate = find_and_rotate_walk(merged_walk, next_walk)

        # If merge failed, find_and_rotate_walk returns old walk unchanged
        if merged_candidate == old_merged:
            not_merged_walks.append(next_walk)
        else:
            merged_walk = merged_candidate  # accept new merged walk

    return merged_walk, not_merged_walks

def test_closed_walk(walk):
    
    # Validates a closed walk based on four rules:
    # 1. Alternating edge types: sequence, then non-sequence, etc.
    # 2. No two consecutive nodes.
    # 3. No two consecutive edges.
    # 4. Walk must start and end with the same node.

    # Returns: (is_valid)

    if not walk:
        return False, "❌ Walk is empty."

    if walk[0][0] != "node" or walk[-1][0] != "node":
        return False, "❌ Walk must start and end with a node."

    if walk[0][1] != walk[-1][1]:
        return False, f"❌ Start and end node differ: {walk[0][1]} ≠ {walk[-1][1]}"

    prev_type = None
    prev_edge_type = None

    for i, entry in enumerate(walk):
        entry_type, label = entry

        # Rule 2 & 3: No two nodes or two edges in a row
        if prev_type == entry_type:
            if entry_type == "node":
                
                return False, f"❌ Two consecutive nodes at index {i-1} and {i}. The walk is {walk}"
            else:
                return False, f"❌ Two consecutive edges at index {i-1} and {i}."

        if entry_type == "edge":
            # Rule 1: Enforce alternating edge type
            edge_type = label[0]  # 'e', 'c', or 'd'
            if prev_edge_type is not None:
                if (edge_type == "e" and prev_edge_type == "e") or \
                   (edge_type in ("c", "d") and prev_edge_type in ("c", "d")):
                    return False, f"❌ Edge alternation failed at index {i}: {prev_edge_type} → {edge_type}"
            prev_edge_type = edge_type

        prev_type = entry_type

    return True, "✅ Walk is valid."
###### test if all closed walks are merged 
###### by testing if all edges of ILP are used in the merged final walk
def test_merged_closed_walk(walk, sequence_edges, f_sequence_edge, discordant_edges, f_discordant_edge, 
                         concordant_edges, f_concordant_edge, F,Path):
    
    # Filter and normalize undirected sequence edges
    number_of_sequence_edges_in_solution = 0
    for (i, j) in sequence_edges:
        if Path == False:
            if i == "t" or j == "t": # when we enforce connectivity we should not count seq edge connected to t for a closed walk, 
                continue
        if f_sequence_edge[i, j] > 0:
            multiplicity = round(f_sequence_edge[i, j] / F)
            number_of_sequence_edges_in_solution += multiplicity    
            
    number_of_discordant_edges_in_solution = 0
    for (i, j) in discordant_edges:
        if Path == False:
            if i == "s" or j == "s": # when we enforce connectivity we should not count disc edge connected to s for a closed walk but we should consider it for a path which is a s-t-walk
                continue
        if f_discordant_edge[i, j] > 0:
            multiplicity = round(f_discordant_edge[i, j] / F)
            number_of_discordant_edges_in_solution += multiplicity   
    
    number_of_concordant_edges_in_solution = 0
    for (i, j) in concordant_edges:
        # if i == "t" or j == "t":
        #     continue
        if f_concordant_edge[i, j] > 0:
            multiplicity = round(f_concordant_edge[i, j] / F)
            number_of_concordant_edges_in_solution += multiplicity    
    number_of_edges_in_ILP_Solution = number_of_sequence_edges_in_solution+number_of_concordant_edges_in_solution+number_of_discordant_edges_in_solution
    number_of_edges_in_merged_walk = (len(walk)-1)/2
    if number_of_edges_in_ILP_Solution == number_of_edges_in_merged_walk: # this means that all the edges of ILP solution are used in the merged walk
       return True
    else:
       return False

########## In case od two disconnected cycles, we might be able to increase the CN of one of them
def Increase_CN(walk, F, capacity_sequence_edges, capacity_concordant_edges, capacity_discordant_edges):
    capacities_of_all_edge_types_in_walk = []

    for edge_type, edge in walk:
        if edge_type == 'sequence edge':
            cap = capacity_sequence_edges.get(edge, 0)
        elif edge_type == 'concordant edge':
            cap = capacity_concordant_edges.get(edge, 0)
        elif edge_type == 'discordant edge':
            cap = capacity_discordant_edges.get(edge, 0)
        else:
            continue  # ignore nodes etc.
        
        if cap > 0:  # only consider positive capacities
            capacities_of_all_edge_types_in_walk.append(cap)

    if not capacities_of_all_edge_types_in_walk:
        return 0.0  # no positive capacities

    min_capacity = min(capacities_of_all_edge_types_in_walk)
    print("min capacity is:", min_capacity)
    print("F is:", F)

    delta_F = max(0.0, min_capacity - F)
    return delta_F
########## parsing function
def parse_node_colon_format(node_str):
    
     # Parse a node name like 'chr15:98778369-' into ('chr15', 98778369)
    
    match = re.match(r"(chr[\w\d]+):(\d+)[+-]", node_str)
    if match:
        chrom = match.group(1)
        pos = int(match.group(2))
        return chrom, int(pos)
    return None, None
############# convert closed walk to segment string
def convert_walk_to_segment_string(walk, segments, copy_count):
    segment_path = []

    for i in range(len(walk)):
        entry = walk[i]
        if entry[0] != "sequence edge":
            continue

        node1, node2 = entry[1]

        if i == 0 or i + 1 >= len(walk):
            continue
        if walk[i - 1][0] != "node" or walk[i + 1][0] != "node":
            continue

        prev_node = walk[i - 1][1]
        next_node = walk[i + 1][1]

        chrom1, pos1 = parse_node_colon_format(node1)
        chrom2, pos2 = parse_node_colon_format(node2)

        if chrom1 != chrom2:
            continue

        #for seg_id, (seg_chrom, seg_start, seg_end) in segments.items():
        for seg_id, seg_info in segments.items():
            seg_chrom = seg_info["chrom"]
            seg_start = seg_info["start"]
            seg_end   = seg_info["end"]
            if seg_chrom != chrom1:
                continue

            # Match either forward or reverse segment definition
            if {seg_start, seg_end} == {pos1, pos2}:
                # Determine direction based on strand
                from_strand = prev_node[-1]
                to_strand = next_node[-1]

                if from_strand == '+' and to_strand == '-':
                    segment_path.append(f"{seg_id}-")
                elif from_strand == '-' and to_strand == '+':
                    segment_path.append(f"{seg_id}+")
                break

    return f"Cycle={1};Copy_count={copy_count};Segments={','.join(segment_path)}"
#### check if any of the path constraints are subsequents of walks
def add_e_prefix(walk_raw):
    #Convert ['2+','5-']  → ['e2+','e5-']  (keeps any existing prefix).
    return [f"e{seg}" if seg and seg[0] not in "ecd" else seg for seg in walk_raw if seg]
#################
def reverse_complement(edge_list):
    if not edge_list:
        print("⚠️ Warning: reverse_complement called with empty edge_list.")
        return []

    flip = {"+": "-", "-": "+"}
    result = []
    for e in reversed(edge_list):
        if len(e) < 2 or e[-1] not in flip:
            print(f" Malformed edge label: '{e}' (expected format like 'e12+' or 'c5-')")
            continue  # or raise an error if you prefer to fail fast
        result.append(f"{e[:-1]}{flip[e[-1]]}")
    return result
####################
def is_circular_subsequence(sub, seq):
    ###Check if sub is a contiguous slice of seq (linear or wrapping once).
    n, m = len(seq), len(sub)
    if m > n:
        return False

    # Check all linear slices
    for i in range(n - m + 1):
        if seq[i:i + m] == sub:
            return True

    # Now check wraparound (sub spans the end and the start)
    for i in range(1, m):
        # Last i elements + first (m-i) elements
        if seq[-i:] + seq[:m - i] == sub:
            return True
    return False
#### ---------- main matcher ------------------------------------------------------
def match_e_constraints(walk_raw, path_constraints):
    walk = add_e_prefix(walk_raw)
    walk_rc = reverse_complement(walk)

    matched = []
    for pc in path_constraints:
        if not (isinstance(pc, dict)
                and isinstance(pc.get("constraint_string"), str)
                and "index" in pc):
            continue

        e_edges = re.findall(r"e\d+[+-]", pc["constraint_string"])
        if not e_edges:
            continue
        e_edges_rc = reverse_complement(e_edges)

        if (is_circular_subsequence(e_edges,    walk) or
            is_circular_subsequence(e_edges_rc, walk) or
            is_circular_subsequence(e_edges,    walk_rc) or
            is_circular_subsequence(e_edges_rc, walk_rc)):
            matched.append(pc["index"])

    return matched
######## These function are parsing
def extract_chr(pos):
    #Extract 'chrN' from a string like 'chr3:20660863+', or use full string/int otherwise.
    if isinstance(pos, str):
        match = re.match(r'^(chr[\w]+):', pos)
        return match.group(1) if match else pos
    return str(pos)

def remove_chr_prefix(pos):
    #Remove 'chrN:' prefix if present, else return as string.
    if isinstance(pos, str):
        return re.sub(r'^chr[\w]+:', '', pos)
    return str(pos)
######## Write the ILP output file
def write_ILP_output(output_path,model,sequence_edges,f_sequence_edge,X_sequence_edge,concordant_edges,f_concordant_edge,
                     X_concordant_edge,discordant_edges,f_discordant_edge, X_discordant_edge,P,path_constraints):
    eps = 1e-5
    with open(output_path, "w") as f:
        if model.status == GRB.OPTIMAL:
            f.write(f"Optimization took {end_time - start_time:.2f} seconds.\n")
            # Header with the new "X" column
            f.write(f"{'Edge':<15}{'Amplicon':<15}{'Start position':<20}{'End position':<20}{'X':<12}{'Flow':<10}{'Length':<12}{'K':<8}\n")
            f.write("=" * 120 + "\n")
    
            for (i, j) in sequence_edges:
                if f_sequence_edge[i, j] > eps:
                 # if i != "t" and j != "t":
                    chrom = extract_chr(i)
                    start = remove_chr_prefix(i)
                    end = remove_chr_prefix(j)
                    # Get X_sequence_edge value
                    X_value = X_sequence_edge[i, j]
                    f.write(f"{'Sequence':<15}{chrom:<15}{start:<20}{end:<20}{X_sequence_edge[i, j]:<12.1f}{f_sequence_edge[i, j]:<10.1f}{length_sequence_edges[i, j]:<12.1f}\n")
            f.write("=" * 120 + "\n")
            for (i, j) in concordant_edges:
                if f_concordant_edge[i, j] > eps:
                    chrom = extract_chr(i)
                    start = remove_chr_prefix(i)
                    end = remove_chr_prefix(j)
                    # Get X_concordant_edge value
                    X_value = X_concordant_edge[i, j]
                    f.write(f"{'Concordant':<15}{chrom:<15}{start:<20}{end:<20}{X_concordant_edge[i, j]:<12.1f}{f_concordant_edge[i, j]:<10.1f}\n")
            # f.write("=" * 120 + "\n")
            # for (i, j) in discordant_edges:
            #     if f_discordant_edge[i, j].x > 0:
            #       if i != "s" and j != "s":
            #         chrom = extract_chr(i)
            #         start = remove_chr_prefix(i)
            #         end = remove_chr_prefix(j)
            #         # Get X_discordant_edge value
            #         X_value = X_discordant_edge[i, j].x
            #         f.write(f"{'Discordant':<15}{chrom:<15}{start:<20}{end:<20}{X_discordant_edge[i, j].x:<12.1f}{f_discordant_edge[i, j].x:<10.1f}\n")
            f.write("=" * 120 + "\n")
            for (i, j) in discordant_edges:
                if f_discordant_edge[i, j] > eps:
                #  if i != "s" and j != "s":
                    chrom = extract_chr(i)
                    start = remove_chr_prefix(i)
                    end = remove_chr_prefix(j)
                    # Get X_discordant_edge value
                    X_value = X_discordant_edge[i, j]
                    
                    edge_key = tuple(sorted((i, j)))
                    K_val = G[edge_key[0]][edge_key[1]].get("K", "NA")
                    
                    f.write(f"{'Discordant':<15}{chrom:<15}{start:<20}{end:<20}{X_discordant_edge[i, j]:<12.1f}{f_discordant_edge[i, j]:<10.1f}{'-':<12}{str(K_val):<8}\n")

            
            
            
            f.write("\nSatisfied path constraints :\n")
            f.write("=" * 50 + "\n")
            for k in range(len(path_constraints)):
                if P[k] > 0.5:  # Check if path constraint k is satisfied
                    support = path_constraints[k]["support"]
                    path_str = ", ".join(path_constraints[k]["paths"])
                    f.write(f"path constraint {k+1} : {path_str}, Support: {support}\n")
                    
        
        else:
            f.write("No optimal solution found.\n")

####### Write the ILP cycle file
def write_ILP_cycle_file(segments, new_output_file,final_segment_str, Iteration,
                         Final_Closed_Walk_segments, Final_closed_walk_test, Final_closed_walk_completion_test, error_in_closed_walks,
                         Path_Constraints, matched_index_set,max_attempts):
    
    parts = final_segment_str.split(';', 1)  # split at first semicolon
    final_segment_str = f"Cycle={Iteration};{parts[1]}"
    # --------------- CASE 1: FIRST ITERATION ------------------
    if Iteration == 1:

        with open(new_output_file, "w") as out_f:
            out_f.write("List of cycle segments\n")

            for seg_id in sorted(segments.keys()):
                seg = segments[seg_id]
                out_f.write(
                    f"Segment\t{seg['id']}\t{seg['chrom']}\t{seg['start']}\t{seg['end']}\n"
                )

            # Write header
            out_f.write("List of longest subpath constraints\n")

            # Write constraints with satisfaction status
            for pc in Path_Constraints:
                status = "Satisfied" if pc["index"] in matched_index_set else "Unsatisfied"
                out_f.write(
                    f"Path constraint\t{pc['index']}\t{pc['constraint_string']}"
                    f"\tSupport={pc['support']}\t{status}\n"
                )

            # Write section marker
            out_f.write("List of extracted cycles/paths\n")

            # Write first cycle
            out_f.write(final_segment_str)
            path_str = ",".join(str(i) for i in sorted(matched_index_set))
            out_f.write(f";Path_constraints_satisfied={path_str}\n")

    # --------------- CASE 2: SUBSEQUENT ITERATIONS ------------
    else:
        with open(new_output_file, "a") as out_f:

            # Append new cycle only
            out_f.write(final_segment_str)
            path_str = ",".join(str(i) for i in sorted(matched_index_set))
            out_f.write(f";Path_constraints_satisfied={path_str}\n")
            
####### Write the meesage about the solutions 
def write_Messages_file(Solution_Type,output_file, Final_closed_walk_test, 
                        Final_closed_walk_completion_test, error_in_closed_walks,max_attempts,Iteration):
    if Iteration==1:
        with open(output_file, "w") as out_f:
            # test the validation of final closed walk
            if not Final_closed_walk_test:
                out_f.write(f"Final closed walk is invalid: {test_message}\n\n")
                
            # test if the final closed walk is using all the edges of the ILP solution    
            if not Final_closed_walk_completion_test and Solution_Type=="Cycle":
                out_f.write(f"Cycle {Iteration} is a connected component from a disconnected solution. \n\n")
            if not Final_closed_walk_completion_test and Solution_Type=="Path":
                out_f.write(f"Path {Iteration} is a connected component from a disconnected solution. \n\n")

            
            if error_in_closed_walks == 1:
                out_f.write(
                f"There is an error in the cycle printed below because after {max_attempts} attempts "
                f"the heuristic failed to find valid closed walks. \n\nYou should try the heuristic more times.\n\n"
            )    
    else:
        with open(output_file, "a") as out_f:
            
            # test the validation of final closed walk
            if not Final_closed_walk_test:
                out_f.write(f"Final closed walk is invalid: {test_message}\n\n")
                
            # test if the final closed walk is using all the edges of the ILP solution    
            if not Final_closed_walk_completion_test and Solution_Type=="Cycle":
                out_f.write(f"Cycle {Iteration} is a connected component from a disconnected solution. \n\n")
            if not Final_closed_walk_completion_test and Solution_Type=="Path":
                out_f.write(f"Path {Iteration} is a connected component from a disconnected solution.\n\n")
    
            
            if error_in_closed_walks == 1:
                out_f.write(
                f"There is an error in the cycle printed below because after {max_attempts} attempts "
                f"the heuristic failed to find valid closed walks. \n\nYou should try the heuristic more times.\n\n"
            ) 
###### compute Length weighted copy number solution
def compute_Length_weighted_copy_number_solution(sequence_edges,f_sequence_edge,length_sequence_edges, Final_closed_walk):
    # valid_edges = [
    #     tuple(sorted((i, j)))
    #     for (i, j) in sequence_edges
    #     if i != "t" and j != "t" and i != "s" and j != "s" and not i.isdigit() and not j.isdigit()
    # ]
    # # Remove duplicates in case (i, j) and (j, i) both appear
    # valid_edges = list(set(valid_edges))
    
    # weight_of_solution = sum(
    #     f_sequence_edge[edge] * length_sequence_edges[edge]
    #     for edge in valid_edges
    # )
    # 1️⃣ Extract only the sequence edges from Final_closed_walk
    walk_sequence_edges = [
        edge for (etype, edge) in Final_closed_walk
        if etype == "sequence edge"
    ]

    # 2️⃣ Normalize the edges (sorted tuples) to match keys in f_sequence_edge / length_sequence_edges
    valid_edges = [
        tuple(sorted(edge))
        for edge in walk_sequence_edges
        if edge[0] not in ("s", "t") and edge[1] not in ("s", "t")
           and not edge[0].isdigit() and not edge[1].isdigit()
    ]

    # 3️⃣ Remove duplicates in case of reverse direction
    valid_edges = list(set(valid_edges))

    # 4️⃣ Sum up weight
    weight_of_solution = sum(
        f_sequence_edge[e] * length_sequence_edges[e]
        for e in valid_edges
        if e in f_sequence_edge and e in length_sequence_edges
    )
    
    return weight_of_solution                
###### compute Length weighted copy number solution, here we do not consider seq edges with CN<5. (THE CN of the graph)
def compute_LWCN_solution_excludeing_some_segments(capacity_sequence_edges,sequence_edges,f_sequence_edge,length_sequence_edges, Final_closed_walk):
    # valid_edges = [
    #     tuple(sorted((i, j)))
    #     for (i, j) in sequence_edges
    #     if i != "t" and j != "t" and i != "s" and j != "s" and not i.isdigit() and not j.isdigit()
    # ]
    # # Remove duplicates in case (i, j) and (j, i) both appear
    # valid_edges = list(set(valid_edges))
    
    # weight_of_solution = sum(
    #     f_sequence_edge[edge] * length_sequence_edges[edge]
    #     for edge in valid_edges
    # )
    # 1️⃣ Extract only the sequence edges from Final_closed_walk
    walk_sequence_edges = [
        edge for (etype, edge) in Final_closed_walk
        if etype == "sequence edge"
    ]

    # 2️⃣ Normalize the edges (sorted tuples) to match keys in f_sequence_edge / length_sequence_edges
    valid_edges = [
        tuple(sorted(edge))
        for edge in walk_sequence_edges
        if edge[0] not in ("s", "t") and edge[1] not in ("s", "t")
           and not edge[0].isdigit() and not edge[1].isdigit()
           and capacity_sequence_edges[edge]>5
    ]

    # 3️⃣ Remove duplicates in case of reverse direction
    valid_edges = list(set(valid_edges))

    # 4️⃣ Sum up weight
    weight_of_solution = sum(
        f_sequence_edge[e] * length_sequence_edges[e]
        for e in valid_edges
        if e in f_sequence_edge and e in length_sequence_edges
    )
    
    return weight_of_solution
############# write the LWCN of cycles in LWCN output file
def write_LWCN_cycle_in_file(metrics_output_path
                             ,weight_cycle_ILP,weight_graph,weight_ratio_cycle_ILP, length_cycle_ILP
                             ,weight_cycle_ILP_excluding_some_segments,weight_graph_excluding_some_segments,weight_ratio_cycle_ILP_excluding_some_segments, length_cycle_ILP_excluding_some_segments
                             ,F,Iteration,Num_paths):
    #with open(metrics_output_path, "a") as f:
        # Only write Coral cycle info for the first iteration
        if Iteration == 1:
            with open(metrics_output_path, "w") as f:
                i=0
                # for i, cycle in enumerate(cycle_weights, start=1):
                #     f.write(f"Reconstructed Cycle {i}:\n")
                #     f.write(f"{'Copy count':<35}: {cycle['copy_count']:.4f}\n")
                #     f.write(f"{'Total segment length':<35}: {cycle['total_length']}\n")
                #     f.write(f"{'Length weighted copy number':<35}: {cycle['weight']:.2f}\n")
                #     f.write(f"{'length weighted copy number of (original) graph':<40}: {weight_graph:.2f}\n")
                #     f.write(f"{'Weight ratio of reconstructed cycle':<35}: {cycle['weight_ratio']:.4f}\n")
                    
                #     f.write(f"{'length of reconstructed cycle excluding some segments':<40}: {cycle['Length_excluding_segments']:.2f}\n")
                #     f.write(f"{'length weighted copy number of reconstructed cycle excluding some segments':<40}: {cycle['weight_excluding_segments']:.2f}\n")
                #     f.write(f"{'length weighted copy number of (original) graph excluding some segments':<40}: {weight_graph_excluding_some_segments:.2f}\n")
                #     f.write(f"{'weight ratio of reconstructed cycle excluding some segments':<40}: {cycle['weight_ratio_excluding_segments']:.4f}\n")

                #     f.write("\n")
                # for j, path in enumerate(path_weights, start=1):
                #     f.write(f"Reconstructed path {j+i}:\n")
                #     f.write(f"{'Copy count':<35}: {path['copy_count']:.4f}\n")
                #     f.write(f"{'Total segment length':<35}: {path['total_length']}\n")
                #     f.write(f"{'Length weighted copy number':<35}: {path['weight']:.2f}\n")
                #     f.write(f"{'length weighted copy number of (original) graph':<40}: {weight_graph:.2f}\n")
                #     f.write(f"{'Weight ratio of reconstructed path':<35}: {path['weight_ratio']:.4f}\n")
                    
                #     f.write(f"{'length of reconstructed path excluding some segments':<40}: {path['Length_excluding_segments']:.2f}\n")
                #     f.write(f"{'length weighted copy number of reconstructed path excluding some segments':<40}: {path['weight_excluding_segments']:.2f}\n")
                #     f.write(f"{'length weighted copy number of (original) graph excluding some segments':<40}: {weight_graph_excluding_some_segments:.2f}\n")
                #     f.write(f"{'weight ratio of reconstructed path excluding some segments':<40}: {path['weight_ratio_excluding_segments']:.4f}\n")

                    
                #     f.write("\n")
                    
                    
                f.write("\n")
                f.write(f"{'ILP Cycle ':}{Iteration}: \n")
                f.write(f"{'Copy count':<40}: {F:.2f}\n")
                f.write(f"{'length of ILP cycle':<40}: {length_cycle_ILP:.2f}\n")
                f.write(f"{'length weighted copy number of ILP cycle':<40}: {weight_cycle_ILP:.2f}\n")
                f.write(f"{'length weighted copy number of (original) graph':<40}: {weight_graph:.2f}\n")
                f.write(f"{'weight ratio of ILP cycle':<40}: {weight_ratio_cycle_ILP[Iteration-Num_paths-1]:.4f}\n")
                #f.write(f"{'weight ratio of ILP cycle':<40}: {weight_ratio_cycle_ILP[Iteration-1]:.4f}\n")

                f.write(f"{'length of ILP cycle excluding some segments':<40}: {length_cycle_ILP_excluding_some_segments:.2f}\n")
                f.write(f"{'length weighted copy number of ILP cycle excluding some segments':<40}: {weight_cycle_ILP_excluding_some_segments:.2f}\n")
                f.write(f"{'length weighted copy number of (original) graph excluding some segments':<40}: {weight_graph_excluding_some_segments:.2f}\n")
                f.write(f"{'weight ratio of ILP cycle excluding some segments':<40}: {weight_ratio_cycle_ILP_excluding_some_segments[Iteration-Num_paths-1]:.4f}\n")
                #f.write(f"{'weight ratio of ILP cycle excluding some segments':<40}: {weight_ratio_cycle_ILP_excluding_some_segments[Iteration-1]:.4f}\n")

                f.write("\n")
        else:
            with open(metrics_output_path, "a") as f:
                f.write("\n")
                f.write(f"{'ILP Cycle ':}{Iteration}: \n")
                f.write(f"{'Copy count':<40}: {F:.2f}\n")
                f.write(f"{'length of ILP cycle':<40}: {length_cycle_ILP:.2f}\n")
                f.write(f"{'length weighted copy number of ILP cycle':<40}: {weight_cycle_ILP:.2f}\n")
                f.write(f"{'length weighted copy number of (original) graph':<40}: {weight_graph:.2f}\n")
                f.write(f"{'weight ratio of ILP cycle':<40}: {weight_ratio_cycle_ILP[Iteration-Num_paths-1]:.4f}\n")
                #f.write(f"{'weight ratio of ILP cycle':<40}: {weight_ratio_cycle_ILP[Iteration-1]:.4f}\n")

                f.write(f"{'length of ILP cycle excluding some segments':<40}: {length_cycle_ILP_excluding_some_segments:.2f}\n")
                f.write(f"{'length weighted copy number of ILP cycle excluding some segments':<40}: {weight_cycle_ILP_excluding_some_segments:.2f}\n")
                f.write(f"{'length weighted copy number of (original) graph excluding some segments':<40}: {weight_graph_excluding_some_segments:.2f}\n")
                f.write(f"{'weight ratio of ILP cycle excluding some segments':<40}: {weight_ratio_cycle_ILP_excluding_some_segments[Iteration-Num_paths-1]:.4f}\n")
                #f.write(f"{'weight ratio of ILP cycle excluding some segments':<40}: {weight_ratio_cycle_ILP_excluding_some_segments[Iteration-1]:.4f}\n")

                f.write("\n")
        #f.write("\n")
######## Update the graph
def update_the_graph( concordant_edges, f_concordant_edge, X_concordant_edge, capacity_concordant_edges,
                     discordant_edges, f_discordant_edge, X_discordant_edge, capacity_discordant_edges,
                     sequence_edges, f_sequence_edge, X_sequence_edge, capacity_sequence_edges):
    # # collect all edges from the Walk into a set for fast lookup
    # walk_edges = set()
    # for edge_type, edge in Walk:
    #     if edge_type in ('concordant edge', 'discordant edge', 'sequence edge'):
    #         walk_edges.add(tuple(sorted(edge)))

    for (i,j) in concordant_edges:
        edge = tuple(sorted((i, j)))
        if f_concordant_edge[edge] > 0 and X_concordant_edge[edge] > 0 : #and edge in walk_edges:
            capacity_concordant_edges[edge] -= f_concordant_edge[edge]
            
            
    for (i,j) in discordant_edges:
        edge = tuple(sorted((i, j)))
        if f_discordant_edge[edge] > 0 and X_discordant_edge[edge] > 0 : #and edge in walk_edges:
            capacity_discordant_edges[edge] -= f_discordant_edge[edge]
        
    for (i,j) in sequence_edges:
        edge = tuple(sorted((i, j)))
        if f_sequence_edge[edge] > 0 and X_sequence_edge[edge] > 0: # and edge in walk_edges:
            capacity_sequence_edges[edge] -= f_sequence_edge[edge]
    return capacity_concordant_edges, capacity_discordant_edges, capacity_sequence_edges
############ Remove s and t nodes and the their edges
def remove_s_t_from_graph(
    G,
    capacity_discordant_edges,
    capacity_sequence_edges,
    K_discordant_edges,
    length_sequence_edges
):
    """
    Remove nodes 's' and 't' and all their incident edges from the graph G.
    Also remove the associated entries from the capacity and length dictionaries.
    """

    # Make fresh dictionaries so original ones stay untouched
    new_capacity_discordant_edges = dict(capacity_discordant_edges)
    new_capacity_sequence_edges   = dict(capacity_sequence_edges)
    new_K_discordant_edges        = dict(K_discordant_edges)
    new_length_sequence_edges     = dict(length_sequence_edges)

    # 1. Collect all edges to 's' and 't'
    s_edges = [(u, v) for u, v in G.edges() if 's' in (u, v)]
    t_edges = [(u, v) for u, v in G.edges() if 't' in (u, v)]

    # 2. Remove these edges from G
    for u, v in s_edges + t_edges:
        if G.has_edge(u, v):
            G.remove_edge(u, v)

    # 3. Remove the nodes 's' and 't' from G
    if 's' in G:
        G.remove_node('s')
    if 't' in G:
        G.remove_node('t')

    # 4. Remove dictionary entries for edges connected to 's' or 't'
    for d in [new_capacity_discordant_edges, new_capacity_sequence_edges,
              new_K_discordant_edges, new_length_sequence_edges]:
        for edge in list(d.keys()):
            if 's' in edge or 't' in edge:
                d.pop(edge, None)

    # 5. Recompute sets of edges after removal
    discordant_edges = {tuple(sorted((u, v))) for u, v, attr in G.edges(data=True)
                        if attr.get('type') == 'discordant'}
    sequence_edges   = {tuple(sorted((u, v))) for u, v, attr in G.edges(data=True)
                        if attr.get('type') == 'sequence'}

    return G, G.nodes, discordant_edges, new_capacity_discordant_edges, sequence_edges, \
           new_capacity_sequence_edges, new_K_discordant_edges, new_length_sequence_edges
########## creat graph with nodes "s" and "t"
def creat_s_t_graph(G,capacity_discordant_edges,K_discordant_edges):
    G.add_node('s')
    G.add_node('t')
    # Make a copy so we don’t mutate caller dict
    new_capacity_discordant_edges = dict(capacity_discordant_edges)               
    for node in G.nodes():
        if node in {'s', 't'} or (isinstance(node, str) and node.isdigit()):
            continue
    
        
        dc_cns = [
            data.get('capacity', 0)
            for _, _, data in get_incident_edges(G, node, edge_type_filter='discordant')
        ]
        
        
        if dc_cns:  # prefer discordant/concordant if present
            capacity = max(dc_cns)
        else:
            seq_cns = [
                data.get('capacity', 0)
                for _, _, data in get_incident_edges(G, node, edge_type_filter='sequence')
            ]
            capacity = max(seq_cns) 
        
        G.add_edge('s', node, type='discordant', capacity=capacity, Read_Count=1000, K=1)
        G.add_edge('t', node, type='discordant', capacity=capacity, Read_Count=1000, K=1)
        
        
    
        new_capacity_discordant_edges[tuple(sorted(('s', node)))] = capacity
        new_capacity_discordant_edges[tuple(sorted(('t', node)))] = capacity
    
    discordant_edges = {
        tuple(sorted((u, v))) 
        for u, v, attr in G.edges(data=True) 
        if attr.get('type') == 'discordant'}
    
    # capacity_discordant_edges = {
    #     tuple(sorted((u, v))): attr['capacity']
    #     for u, v, attr in G.edges(data=True)
    #     if attr.get('type') == 'discordant'}
    # KEEP updated capacities, don’t rebuild from G
    capacity_discordant_edges = new_capacity_discordant_edges
    
    K_discordant_edges = {
        tuple(sorted((u, v))): attr["K"]
        for u, v, key, attr in G.edges(keys=True, data=True)
        if attr.get("type") == "discordant" and "K" in attr}
    return G,G.nodes,discordant_edges,capacity_discordant_edges,K_discordant_edges
####### create s t graph like Coral. Connect s and t to the beginning and end of intervals only
def creat_s_t_graph_Coral(G, capacity_discordant_edges, K_discordant_edges, sequence_edges):
    """
    Add source 's' and sink 't' nodes connected to the start/end of intervals
    derived from connected sequence edges.
    """
    G.add_node('s')
    G.add_node('t')
    new_capacity_discordant_edges = dict(capacity_discordant_edges)

    # --- Detect intervals ---
    # Convert edges to list of (chrom, start, end) tuples
    seq_list = []
    for u, v in sequence_edges:
        chrom_u, start_u = u.split(":")
        chrom_v, end_v = v.split(":")
        start_u = int(start_u.rstrip("-+"))
        end_v = int(end_v.rstrip("-+"))
        seq_list.append((chrom_u, start_u, end_v, u, v))

    # Sort by chromosome, then start position
    seq_list.sort(key=lambda x: (x[0], x[1]))

    intervals = []
    if seq_list:
        chrom, start, end, u_node, v_node = seq_list[0]
        interval_start_node = u_node
        interval_end_node = v_node
        prev_chrom = chrom
        prev_end = end

        for chrom, start, end, u_node, v_node in seq_list[1:]:
            if chrom == prev_chrom and start == prev_end + 1:
                # extend current interval
                interval_end_node = v_node
                prev_end = end
            else:
                # save previous interval
                intervals.append((interval_start_node, interval_end_node))
                # start new interval
                interval_start_node = u_node
                interval_end_node = v_node
                prev_end = end
                prev_chrom = chrom
        intervals.append((interval_start_node, interval_end_node))  # last interval

    # --- Add s/t edges for each interval ---
    for start_node, end_node in intervals:
        # s -> start
        dc_cns = [
            data.get('capacity', 0)
            for _, _, data in get_incident_edges(G, start_node, edge_type_filter='discordant')
        ]
        if dc_cns:
            capacity = max(dc_cns)
        else:
            seq_cns = [
                data.get('capacity', 0)
                for _, _, data in get_incident_edges(G, start_node, edge_type_filter='sequence')
            ]
            capacity = max(seq_cns) if seq_cns else 0
        G.add_edge('s', start_node, type='discordant', capacity=capacity, Read_Count=1000, K=1)
        new_capacity_discordant_edges[tuple(sorted(('s', start_node)))] = capacity
        G.add_edge('t', start_node, type='discordant', capacity=capacity, Read_Count=1000, K=1)
        new_capacity_discordant_edges[tuple(sorted(('t', start_node)))] = capacity

        # end -> t
        dc_cns = [
            data.get('capacity', 0)
            for _, _, data in get_incident_edges(G, end_node, edge_type_filter='discordant')
        ]
        if dc_cns:
            capacity = max(dc_cns)
        else:
            seq_cns = [
                data.get('capacity', 0)
                for _, _, data in get_incident_edges(G, end_node, edge_type_filter='sequence')
            ]
            capacity = max(seq_cns) if seq_cns else 0
        G.add_edge(end_node, 's', type='discordant', capacity=capacity, Read_Count=1000, K=1)
        new_capacity_discordant_edges[tuple(sorted((end_node, 's')))] = capacity
        G.add_edge(end_node, 't', type='discordant', capacity=capacity, Read_Count=1000, K=1)
        new_capacity_discordant_edges[tuple(sorted((end_node, 't')))] = capacity

    # --- Update sets ---
    discordant_edges = {
        tuple(sorted((u, v))) 
        for u, v, attr in G.edges(data=True) 
        if attr.get('type') == 'discordant'
    }
    capacity_discordant_edges = new_capacity_discordant_edges
    K_discordant_edges = {
        tuple(sorted((u, v))): attr["K"]
        for u, v, key, attr in G.edges(keys=True, data=True)
        if attr.get("type") == "discordant" and "K" in attr
    }

    return G, G.nodes, discordant_edges, capacity_discordant_edges, K_discordant_edges

###########
def build_path_model(
                G,
                nodes,
                concordant_edges,
                capacity_concordant_edges,
                discordant_edges,
                capacity_discordant_edges,
                read_count_discordant_edges,
                K_discordant_edges,
                sequence_edges,
                capacity_sequence_edges,
                length_sequence_edges,
                p_sequence_edges,
                p_concordant_edges,
                p_discordant_edges,
                path_constraints,
                enforce_connectivity,
                gamma):
    model = Model("ILP_CoRAL_cycle")
    # Decision variables   
    X_concordant_edge = model.addVars(concordant_edges, vtype=GRB.BINARY, name="X_concordant_edge")
    X_discordant_edge = model.addVars(discordant_edges, vtype=GRB.BINARY, name="X_discordant_edge")
    X_sequence_edge = model.addVars(sequence_edges, vtype=GRB.BINARY, name="X_sequence_edge")   
    f_concordant_edge = model.addVars(concordant_edges, vtype=GRB.CONTINUOUS, lb=0, name="f_concordant_edge")
    f_discordant_edge = model.addVars(discordant_edges, vtype=GRB.CONTINUOUS, lb=0, name="f_discordant_edge")
    f_sequence_edge = model.addVars(sequence_edges, vtype=GRB.CONTINUOUS, lb=0, name="f_sequence_edge")    
    F = model.addVar(name="f_min", lb=0, vtype=GRB.CONTINUOUS)  
    P = model.addVars(len(path_constraints), vtype=GRB.BINARY, name="P_path_constraint")
    y_discordant = {}
    for i, j in discordant_edges:
        edge_key = tuple(sorted((i, j)))
        for m in range(1, K_discordant_edges[edge_key] + 1):
            y_discordant[edge_key + (m,)] = model.addVar(name=f"y_discordant_{edge_key[0]}_{edge_key[1]}_{m}", vtype=GRB.BINARY)
    z_discordant = {}
    for i, j in discordant_edges:
        edge_key = tuple(sorted((i, j)))
        for m in range(1, K_discordant_edges[edge_key] + 1):
            z_discordant[edge_key + (m,)] = model.addVar(lb=0, name=f"z_discordant_{edge_key[0]}_{edge_key[1]}_{m}", vtype=GRB.CONTINUOUS)

    
    ### Objective function
    if len(path_constraints)>0:
        multiplier= gamma* sum(capacity_sequence_edges[tuple(sorted((i, j)))] * length_sequence_edges[tuple(sorted((i, j)))] for (i, j) in sequence_edges)/ len(path_constraints)
    else:
        multiplier=0
    
   ### Objective: Maximize total flow on sequence edges weighted by length and number of path constraints satisfied
    model.setObjective(
        sum(f_sequence_edge[tuple(sorted((i, j)))] * length_sequence_edges[tuple(sorted((i, j)))] for (i, j) in sequence_edges)+ multiplier*sum(P[k] for k in range(len(path_constraints))),
        GRB.MAXIMIZE
    )
    # Constraint 1: Copy number of an edge ≤ capacity × selection variable (for undirected edges)
    
    for (i, j) in concordant_edges:
        edge = tuple(sorted((i, j)))
        model.addConstr(
            f_concordant_edge[edge] <= capacity_concordant_edges[edge] * X_concordant_edge[edge],
            name=f"ConcordantEdgeCapacity_{edge[0]}_{edge[1]}")    
    for (i, j) in discordant_edges: 
        edge = tuple(sorted((i, j)))
        model.addConstr(
            f_discordant_edge[edge] <= capacity_discordant_edges[edge] * X_discordant_edge[edge],
            name=f"DiscordantEdgeCapacity_{edge[0]}_{edge[1]}")    
    for (i, j) in sequence_edges: 
        edge = tuple(sorted((i, j)))
        model.addConstr(
            f_sequence_edge[edge] <= capacity_sequence_edges[edge] * X_sequence_edge[edge],
            name=f"SequenceEdgeCapacity_{edge[0]}_{edge[1]}")
    ### Constraint 3: copy number sequence edge incident = copy number discordant/concordant edge incident 
    ##########################################################################    
    for k in nodes:
        if k in {"s", "t"}:
            continue
    
        # Define undirected incident edges
        incident_edges = {
          tuple(sorted((i, k)))
          for i in nodes
          if tuple(sorted((i, k))) in concordant_edges | discordant_edges | sequence_edges and i != k}
        # Flow into k (sequence edges)
        flow_in_sequence = sum(
            f_sequence_edge[edge]
            for edge in incident_edges
            if edge in sequence_edges and k in edge)
    
        # Flow out of k (concordant + discordant)
        flow_out_concordant = sum(
            f_concordant_edge[edge]
            for edge in incident_edges
            if edge in concordant_edges and k in edge)
        flow_out_discordant = sum(
            f_discordant_edge[edge]
            for edge in incident_edges
            if edge in discordant_edges and k in edge)
    
        model.addConstr(
            flow_in_sequence == flow_out_concordant + flow_out_discordant,
            name=f"FlowConservation2_{k}")
    ### Constraint 4 "s" and "t" edges and the copy number of their incident edges
    flow_out_s_discordant = sum(
            f_discordant_edge[tuple(sorted(("s", k)))]
            for k in nodes
            if tuple(sorted(("s", k))) in discordant_edges
        )
        
        ### Total flow into sink 't' on discordant edges
    flow_in_t_discordant = sum(
            f_discordant_edge[tuple(sorted((k, "t")))]
            for k in nodes
            if tuple(sorted((k, "t"))) in discordant_edges
        )
        
        ### Enforce the flow match
    model.addConstr(
            flow_out_s_discordant == flow_in_t_discordant,
            name="TotalFlowFromS_Equals_FlowIntoT"
        )
        
    model.addConstr(
            sum(
                X_discordant_edge[tuple(sorted(("s", k)))]
                for k in nodes
                if tuple(sorted(("s", k))) in discordant_edges
            ) == 1,
            name="OneEdgeOutOfS"
        )

    model.addConstr(
            sum(
                X_discordant_edge[tuple(sorted((k, "t")))]
                for k in nodes
                if tuple(sorted((k, "t"))) in discordant_edges
            ) == 1,
            name="OneEdgeIntoT"
        )
    ### Constraint 7 #Big-M constant (ensure it's large enough)    
    M = 1000000    
    for edge in concordant_edges:
        i, j = edge
        edge_key = tuple(sorted((i, j)))
        model.addConstr(
            F <= f_concordant_edge[edge_key] + M * (1 - X_concordant_edge[edge_key]),
            name=f"ConcordantEdgeCapacity_{i}_{j}")    
    for edge in discordant_edges:
        i, j = edge
        edge_key = tuple(sorted((i, j)))
        model.addConstr(
            F <= f_discordant_edge[edge_key] + M * (1 - X_discordant_edge[edge_key]),
            name=f"DiscordantEdgeCapacity_{i}_{j}")
    for i, j in discordant_edges:
        edge_key = tuple(sorted((i, j)))
        model.addConstr(
            sum(y_discordant[edge_key + (m,)] for m in range(1, K_discordant_edges[edge_key] + 1)) == X_discordant_edge[edge_key],
            name=f"one_multiple_discordant_{edge_key[0]}_{edge_key[1]}")
    for i, j in discordant_edges:
        edge_key = tuple(sorted((i, j)))
        for m in range(1, K_discordant_edges[edge_key] + 1):
            z = z_discordant[edge_key + (m,)]
            y = y_discordant[edge_key + (m,)]
            model.addConstr(z <= F, name=f"z_leq_fmin_discordant_{edge_key[0]}_{edge_key[1]}_{m}")
            model.addConstr(z <= M * y, name=f"z_leq_My_discordant_{edge_key[0]}_{edge_key[1]}_{m}")
            model.addConstr(z >= F - M * (1 - y), name=f"z_geq_fmin_minus_M_discordant_{edge_key[0]}_{edge_key[1]}_{m}")
            model.addConstr(z >= 0, name=f"z_geq_0_discordant_{edge_key[0]}_{edge_key[1]}_{m}")
    for i, j in discordant_edges:
        edge_key = tuple(sorted((i, j)))
        model.addConstr(
            f_discordant_edge[edge_key] == sum(m * z_discordant[edge_key + (m,)] for m in range(1, K_discordant_edges[edge_key] + 1)),
            name=f"flow_def_discordant_{edge_key[0]}_{edge_key[1]}")
    #### Constarints 9 Path constraints
    for (i, j, k) in p_sequence_edges:
        edge_key = tuple(sorted((i, j)))
        model.addConstr(
            X_sequence_edge[edge_key] >= P[k] * p_sequence_edges[(i, j, k)])
    for (i, j, k) in p_concordant_edges:
        edge_key = tuple(sorted((i, j)))
        model.addConstr(
            X_concordant_edge[edge_key] >= P[k] * p_concordant_edges[(i, j, k)])
    for (i, j, k) in p_discordant_edges:
        edge_key = tuple(sorted((i, j)))
        model.addConstr(
            X_discordant_edge[edge_key] >= P[k] * p_discordant_edges[(i, j, k)])
        
    if enforce_connectivity:
        #####################################################    
        ##### Ensure Connectivity Variables
        #####################################################
        all_edges = concordant_edges | discordant_edges | sequence_edges
        # Generate all directed versions of directed variables on undirected edges
        directed_sequence_edges = {(u, v) for u, v in sequence_edges} | {(v, u) for u, v in sequence_edges}
        directed_concordant_edges = {(u, v) for u, v in concordant_edges} | {(v, u) for u, v in concordant_edges}
        directed_discordant_edges = {(u, v) for u, v in discordant_edges} | {(v, u) for u, v in discordant_edges}
        c_seq  = model.addVars(directed_sequence_edges,   vtype=GRB.BINARY, name="c_seq")
        c_conc = model.addVars(directed_concordant_edges, vtype=GRB.BINARY, name="c_conc")
        c_disc = model.addVars(directed_discordant_edges, vtype=GRB.BINARY, name="c_disc")
        d = model.addVars(G.nodes, vtype=GRB.INTEGER, lb=0, ub=len(G.nodes), name="d")
        
    
        
        ######################################################
        ######### Connectivity constraint 0 #################
        ###### Just about 's' which is the start node #######
        ##### d(s)=1, y_plus ('s', only one node)=1 that is where X_...=1
        ##### y_minu('s', every other node)=0, y_plus(every other node,'s')=0
        ###################################################### 
       
        for u, v in directed_discordant_edges:
            if u != 's':
                continue
            edge = tuple(sorted((u, v)))  # Get the undirected version
        
            # Only add constraint if this edge exists in X variables
            
            if edge in X_discordant_edge:
                x_var = X_discordant_edge[edge]
            else:
                continue  # Skip if edge not in any X set
        
            model.addConstr(
                c_disc[u, v] ==x_var,
                name=f"y_dir_limit_{u}_{v}"
            )  
            model.addConstr(
                c_disc[v,u] ==0,
                name=f"y_dir_limit_{u}_{v}"
            ) 
            
            
        model.addConstr(d['s'] == 1, name='d_start_node_one')
        ######################################################
        ######### Connectivity constraint 1 #################
        ######################################################
        num_nodes = len(G.nodes)
    
        for u in G.nodes:
            incident_x_vars = []
            for v in G.neighbors(u):
                edge = tuple(sorted((u, v)))
                if edge in sequence_edges:
                    incident_x_vars.append(X_sequence_edge[edge])
                if edge in concordant_edges:
                    incident_x_vars.append(X_concordant_edge[edge])
                if edge in discordant_edges:
                    incident_x_vars.append(X_discordant_edge[edge])
            model.addConstr(
                d[u] <= num_nodes * sum(incident_x_vars),
                name=f"only_assign_d_if_used_{u}"
            )
            
        ###################################################
        ###### Connectivity constraint 2 ##################
        ###################################################
                            
        for u, v in sequence_edges:
            edge = tuple(sorted((u, v)))
            model.addConstr(c_seq[u, v] + c_seq[v, u] <= X_sequence_edge[(u, v)], name=f"dir_seq_use_{u}_{v}")
        
        for u, v in concordant_edges:
            edge = tuple(sorted((u, v)))
            model.addConstr(c_conc[u, v] + c_conc[v, u] <= X_concordant_edge[(u, v)], name=f"dir_conc_use_{u}_{v}")
        
        for u, v in discordant_edges:
            edge = tuple(sorted((u, v)))
            model.addConstr(c_disc[u, v] + c_disc[v, u] <= X_discordant_edge[(u, v)], name=f"dir_disc_use_{u}_{v}")
        
        ###################################################
        ######### Connectivity constraint 3 ###############
        ###################################################
        all_nodes = list(G.nodes)
        for u in all_nodes:
            if u == 's':
                continue  # skip source
            incoming_c = []
    
            for v in all_nodes:
                if (v, u) in c_seq:
                    incoming_c.append(c_seq[v, u])
                if (v, u) in c_conc:
                    incoming_c.append(c_conc[v, u])
                if (v, u) in c_disc:
                    incoming_c.append(c_disc[v, u])
        
            incoming_sum = sum(incoming_c)
        
            # Compute total x variables that involve u (undirected edge participation)
            x_incident = []
            for v in G.neighbors(u):
                edge = tuple(sorted((u, v)))  # undirected
                if edge in X_sequence_edge:
                    x_incident.append(X_sequence_edge[edge])
                if edge in X_concordant_edge:
                    x_incident.append(X_concordant_edge[edge])
                if edge in X_discordant_edge:
                    x_incident.append(X_discordant_edge[edge])
        
            x_total = sum(x_incident)
            print('x_total for', u, 'is', x_total)
        
            # Constraint: if node u is in the solution, it must have at least one incoming edge
            model.addConstr(
                incoming_sum >= (1/num_nodes)*x_total,
                name=f"incoming_constraint_{u}"
            )
        ###################################################
        ######### Connectivity constraint 4 ###############
        ###################################################
                        
        for u, v in directed_sequence_edges:
            if (u, v) in c_seq:
                model.addConstr(
                    d[v] - d[u] + num_nodes * (1 - c_seq[u, v]) >= 1,
                    name=f"d_order_sequence_{u}_{v}"
                )
        for u, v in directed_concordant_edges:
            if (u, v) in c_conc:
                model.addConstr(
                    d[v] - d[u] + num_nodes * (1 - c_conc[u, v]) >= 1,
                    name=f"d_order_concordant_{u}_{v}"
                )
        for u, v in directed_discordant_edges:
            if (u, v) in c_disc:
                model.addConstr(
                    d[v] - d[u] + num_nodes * (1 - c_disc[u, v]) >= 1,
                    name=f"d_order_discordant_{u}_{v}"
                )
        ###################################################
        ######### End of constraints ###############
        ###################################################
    return model, {
        "X_concordant_edge": X_concordant_edge,
        "X_discordant_edge": X_discordant_edge,
        "X_sequence_edge": X_sequence_edge,
        "f_concordant_edge": f_concordant_edge,
        "f_discordant_edge": f_discordant_edge,
        "f_sequence_edge": f_sequence_edge,
        "F": F,
        "P": P,
        "y_discordant": y_discordant,
        "z_discordant": z_discordant,
    }
########## Heuristic to create paths
def create_s_t_walks(sequence_edges, f_sequence_edge, discordant_edges, f_discordant_edge, 
                     concordant_edges, f_concordant_edge,max_attempts,segments, F):
    error_in_closed_walks=0 # by default we assume we have a simple path
    node_to_segment_id = {}
    for seg_id, info in segments.items():
        node_start = f"{info['chrom']}:{info['start']}-"
        node_end = f"{info['chrom']}:{info['end']}+"
        node_to_segment_id[node_start] = seg_id
        node_to_segment_id[node_end] = seg_id
    def edge_segment_rank(edge):
        return min(
            node_to_segment_id.get(edge[0], float('inf')),
            node_to_segment_id.get(edge[1], float('inf'))
        )
    # Filter and normalize undirected sequence edges
    solution_sequence_edges = []
    for (i, j) in sequence_edges:
        # if i == "t" or j == "t":
        #     continue
        if f_sequence_edge[i, j] > 0:
            multiplicity = round(f_sequence_edge[i, j] / F)
            solution_sequence_edges.extend([tuple(sorted((i, j)))] * multiplicity)
    
    solution_discordant_edges = []
    for (i, j) in discordant_edges:
        # if i == "s" or j == "s":
        #     continue
        if f_discordant_edge[i, j] > 0:
            multiplicity = round(f_discordant_edge[i, j] / F)
            solution_discordant_edges.extend([tuple(sorted((i, j)))] * multiplicity)
    
    solution_concordant_edges = []
    for (i, j) in concordant_edges:
        if f_concordant_edge[i, j] > 0:
            multiplicity = round(f_concordant_edge[i, j] / F)
            solution_concordant_edges.extend([tuple(sorted((i, j)))] * multiplicity)
    
    remaining_sequence_edges   = Counter(solution_sequence_edges)
    remaining_discordant_edges = Counter(solution_discordant_edges)
    remaining_concordant_edges = Counter(solution_concordant_edges)
    # Step 1: build s→t walk
    first_walk_edges = []
    current_node = 's'
    
    # Start with a discordant edge from s
    start_edge = next(
        e for e in sorted(remaining_discordant_edges, key=edge_segment_rank)
        if 's' in e
    )
    use_edge(remaining_discordant_edges, start_edge)
    first_walk_edges.append(("node", current_node))
    first_walk_edges.append(("discordant edge", start_edge))
    current_node = start_edge[1] if start_edge[0] == current_node else start_edge[0]
    first_walk_edges.append(("node", current_node))
    
    edge_type = "sequence"
    
    while True:
        if edge_type == "sequence":
            match = next((e for e in sorted(remaining_sequence_edges, key=edge_segment_rank) if current_node in e), None)
            if not match:
                break
            use_edge(remaining_sequence_edges, match)
            next_node = match[1] if match[0] == current_node else match[0]
            first_walk_edges.append(("sequence edge", match))
            first_walk_edges.append(("node", next_node))
            current_node = next_node
            edge_type = "non-sequence"
    
        elif edge_type == "non-sequence":
            match = next((e for e in sorted(remaining_discordant_edges, key=edge_segment_rank) if current_node in e), None)
            if match is None:
                match = next((e for e in sorted(remaining_concordant_edges, key=edge_segment_rank) if current_node in e), None)
            if not match:
                break
            edge_label = "discordant edge" if match in remaining_discordant_edges else "concordant edge"
            use_edge(remaining_discordant_edges if edge_label == "discordant edge" else remaining_concordant_edges, match)
            next_node = match[1] if match[0] == current_node else match[0]
            first_walk_edges.append((edge_label, match))
            first_walk_edges.append(("node", next_node))
            current_node = next_node
            edge_type = "sequence"
    
        # Stop condition for s→t walk
        if current_node == 't': #and edge_label == "discordant edge":
            break
    
    s_t_Walk = [first_walk_edges]
    ### if there are still edges remaining we should see if it is possible to merge them in the s-t  path or not
    ### if it is not possible to merge those remaining edges in the s-t path, then we have a disconnected soution which should be a cycle
    Number_Of_Remaining_edges = (
        sum(remaining_sequence_edges.values()) +
        sum(remaining_discordant_edges.values()) +
        sum(remaining_concordant_edges.values())
    )
    if Number_Of_Remaining_edges > 0:
        Closed_Walks_Of_Path, error_in_closed_walks_of_path=create_closed_walks(remaining_sequence_edges, f_sequence_edge, remaining_discordant_edges, f_discordant_edge, 
                             remaining_concordant_edges, f_concordant_edge, max_attempts,segments, F)
        s_t_Walk.extend(Closed_Walks_Of_Path)
        
        for walk in s_t_Walk[1:]:   # skip the first one which is s-t path and is not close
            if any(not test_closed_walk(walk) for walk in Closed_Walks_Of_Path[1:]):
                error_in_closed_walks = 1
            else:
                error_in_closed_walks = 0
        
    return s_t_Walk, error_in_closed_walks

############# convert s-t walk to segment string
def convert_s_t_walk_to_segment_string(walk, segments, copy_count):
    segment_path = []

    for i in range(len(walk)):
        entry = walk[i]
        if entry[0] != "sequence edge":
            continue

        node1, node2 = entry[1]

        if i == 0 or i + 1 >= len(walk):
            continue
        if walk[i - 1][0] != "node" or walk[i + 1][0] != "node":
            continue

        prev_node = walk[i - 1][1]
        next_node = walk[i + 1][1]

        chrom1, pos1 = parse_node_colon_format(node1)
        chrom2, pos2 = parse_node_colon_format(node2)

        if chrom1 != chrom2:
            continue

        #for seg_id, (seg_chrom, seg_start, seg_end) in segments.items():
        for seg_id, seg_info in segments.items():
            seg_chrom = seg_info["chrom"]
            seg_start = seg_info["start"]
            seg_end   = seg_info["end"]
            if seg_chrom != chrom1:
                continue

            # Match either forward or reverse segment definition
            if {seg_start, seg_end} == {pos1, pos2}:
                # Determine direction based on strand
                from_strand = prev_node[-1]
                to_strand = next_node[-1]

                if from_strand == '+' and to_strand == '-':
                    segment_path.append(f"{seg_id}-")
                elif from_strand == '-' and to_strand == '+':
                    segment_path.append(f"{seg_id}+")
                break
    # add 0+ at the beginning and 0- at the end
    segment_path = ["0+"] + segment_path + ["0-"]
    
    return f"Path={1};Copy_count={copy_count};Segments={','.join(segment_path)}"    
####### Write paths in the cycles file
def write_paths_in_ILP_cycle_file(
        segments, new_output_file, path_string, Iteration,
        Path_Constraints, matched_index_set):
    """
    Append paths to ILP cycle/paths output file.
    If file does not yet exist, create and initialize it with
    segments + constraints + header.
    """

    # --------------------------------------------------------
    # Normalize path_string so it always begins with Path=Iteration
    # --------------------------------------------------------
    parts = path_string.split(";", 1)
    if len(parts) > 1:
        path_string = f"Path={Iteration};{parts[1]}"
    else:
        path_string = f"Path={Iteration};{path_string}"

    # ========================================================
    # CASE 1: FILE ALREADY EXISTS → append only the new path
    # ========================================================
    if os.path.exists(new_output_file) and Iteration > 1:

        with open(new_output_file, "a") as out_f:
            out_f.write(path_string)
            path_str = ",".join(str(i) for i in sorted(matched_index_set))
            out_f.write(f";Path_constraints_satisfied={path_str}\n")

        return   # done for this case

    # ========================================================
    # CASE 2: FILE DOES NOT EXIST → create and print headers
    # ========================================================
    with open(new_output_file, "w") as out_f:

        # ------- Print segments -------
        out_f.write("List of cycle segments\n")
        for seg_id in sorted(segments.keys()):
            seg = segments[seg_id]
            out_f.write(
                f"Segment\t{seg['id']}\t{seg['chrom']}\t{seg['start']}\t{seg['end']}\n"
            )

        # ------- Print path constraints -------
        out_f.write("List of longest subpath constraints\n")

        for pc in Path_Constraints:
            status = "Satisfied" if pc["index"] in matched_index_set else "Unsatisfied"
            out_f.write(
                f"Path constraint\t{pc['index']}\t{pc['constraint_string']}"
                f"\tSupport={pc['support']}\t{status}\n"
            )

        # ------- Section header for paths -------
        out_f.write("List of extracted paths\n")

        # ------- Write first path -------
        out_f.write(path_string)
        path_str = ",".join(str(i) for i in sorted(matched_index_set))
        out_f.write(f";Path_constraints_satisfied={path_str}\n")

############# write the LWCN of paths in LWCN output file
def write_LWCN_path_in_file(metrics_output_path
                             ,weight_path_ILP,weight_graph,weight_ratio_path_ILP, length_path_ILP
                             ,weight_path_ILP_excluding_some_segments,weight_graph_excluding_some_segments,weight_ratio_path_ILP_excluding_some_segments, length_path_ILP_excluding_some_segments
                             ,F,Iteration,Num_cycles):
    with open(metrics_output_path, "a") as f:
        # Only write Coral cycle info for the first iteration
        if Iteration == 1:
            with open(metrics_output_path, "w") as f:
                i=0
                # for i, cycle in enumerate(cycle_weights, start=1):
                #     f.write(f"Reconstructed Cycle {i}:\n")
                #     f.write(f"{'Copy count':<35}: {cycle['copy_count']:.4f}\n")
                #     f.write(f"{'Total segment length':<35}: {cycle['total_length']}\n")
                #     f.write(f"{'Length weighted copy number':<35}: {cycle['weight']:.2f}\n")
                #     f.write(f"{'Weight ratio of reconstructed cycle':<35}: {cycle['weight_ratio']:.4f}\n")
                #     f.write("\n")
                # for j, path in enumerate(path_weights, start=1):
                #     f.write(f"Reconstructed path {j+i}:\n")
                #     f.write(f"{'Copy count':<35}: {path['copy_count']:.4f}\n")
                #     f.write(f"{'Total segment length':<35}: {path['total_length']}\n")
                #     f.write(f"{'Length weighted copy number':<35}: {path['weight']:.2f}\n")
                #     f.write(f"{'Weight ratio of reconstructed path':<35}: {path['weight_ratio']:.4f}\n")
                #     f.write("\n")
                f.write("\n")
                f.write(f"{'ILP Path ':}{Iteration}: \n")
                f.write(f"{'Copy count':<40}: {F:.2f}\n")
                f.write(f"{'length of ILP path':<40}: {length_path_ILP:.2f}\n")
                f.write(f"{'length weighted copy number of ILP path':<40}: {weight_path_ILP:.2f}\n")
                f.write(f"{'length weighted copy number of (original) graph':<40}: {weight_graph:.2f}\n")
                f.write(f"{'weight ratio of ILP path':<40}: {weight_ratio_path_ILP[Iteration-Num_cycles-1]:.4f}\n")
                f.write(f"{'length of ILP path excluding some segments':<40}: {length_path_ILP_excluding_some_segments:.2f}\n")
                f.write(f"{'length weighted copy number of ILP path excluding some segments':<40}: {weight_path_ILP_excluding_some_segments:.2f}\n")
                f.write(f"{'length weighted copy number of (original) graph excluding some segments':<40}: {weight_graph_excluding_some_segments:.2f}\n")
                f.write(f"{'weight ratio of ILP path excluding some segments':<40}: {weight_ratio_path_ILP_excluding_some_segments[Iteration-Num_cycles-1]:.4f}\n")
                f.write("\n")
        else:
            with open(metrics_output_path, "a") as f:
                f.write("\n")
                f.write(f"{'ILP Path ':}{Iteration}: \n")
                f.write(f"{'Copy count':<40}: {F:.2f}\n")
                f.write(f"{'length of ILP path':<40}: {length_path_ILP:.2f}\n")
                f.write(f"{'length weighted copy number of ILP path':<40}: {weight_path_ILP:.2f}\n")
                f.write(f"{'length weighted copy number of (original) graph':<40}: {weight_graph:.2f}\n")
                f.write(f"{'weight ratio of ILP path':<40}: {weight_ratio_path_ILP[Iteration-Num_cycles-1]:.4f}\n")
                f.write(f"{'length of ILP path excluding some segments':<40}: {length_path_ILP_excluding_some_segments:.2f}\n")
                f.write(f"{'length weighted copy number of ILP path excluding some segments':<40}: {weight_path_ILP_excluding_some_segments:.2f}\n")
                f.write(f"{'length weighted copy number of (original) graph excluding some segments':<40}: {weight_graph_excluding_some_segments:.2f}\n")
                f.write(f"{'weight ratio of ILP path excluding some segments':<40}: {weight_ratio_path_ILP_excluding_some_segments[Iteration-Num_cycles-1]:.4f}\n")
                f.write("\n")

######## Print the cycles and paths in order of LWCN excluding some segments
def write_all_cycles_and_paths(segments, new_output_file, all_cycles_ordered, all_paths_ordered, path_constraints):
    """
    Writes all cycles and all paths at once, in the order provided.
    If the file doesn't exist, writes the list of segments and path constraints first.
    """

    # import os
    # first_write = not os.path.exists(new_output_file)

    # Open file in write mode if new, else append
   # with open(new_output_file, "w" if first_write else "a") as out_f:
    with open(new_output_file, "w") as out_f:
        # --- Write segments and path constraints if first write ---
       # if first_write:
            # List of cycle segments
        out_f.write("List of cycle segments\n")
        for seg_id in sorted(segments.keys()):
            seg = segments[seg_id]
            out_f.write(
                f"Segment\t{seg['id']}\t{seg['chrom']}\t{seg['start']}\t{seg['end']}\n"
            )

        # List of path constraints
        out_f.write("List of longest subpath constraints\n")
        if path_constraints:
            for pc in path_constraints:
                # For now, no cycle-specific match info; just list constraints
                out_f.write(
                    f"Path constraint\t{pc['index']}\t{pc['constraint_string']}"
                    f"\tSupport={pc['support']}\tN/A\n"
                )

        # Section header for cycles/paths
        out_f.write("List of extracted cycles/paths\n")

        # --- Write all cycles ---
        for cycle in all_cycles_ordered:
            seg_str = ",".join(cycle["Segments"])
            matched_index_set = cycle.get("Matched_Index_Set", set())
            path_str = ",".join(str(i) for i in sorted(matched_index_set))

            out_f.write(
                f"Cycle={cycle['New_Order_Of_Cycle_Number']};"
                f"Copy_count={cycle['CopyNumber']};"
                f"Segments={seg_str};"
                f"Path_constraints_satisfied={path_str}\n"
            )

        # --- Write all paths ---
        for path in all_paths_ordered:
            seg_str = ",".join(path["Segments"])
            matched_index_set = path.get("Matched_Index_Set", set())
            path_str = ",".join(str(i) for i in sorted(matched_index_set))

            out_f.write(
                f"Path={path['New_Order_Of_Path_Number']};"
                f"Copy_count={path['CopyNumber']};"
                f"Segments={seg_str};"
                f"Path_constraints_satisfied={path_str}\n"
            )
#####################################################################
################## The code starts here #############################
#####################################################################

__version__ = "1.0.0"

#ENFORCE_CONNECTIVITY = False #True False
TIME_LIMIT_SECONDS = 2*60*60 #2 * 60 * 60  # 2 hours
delta_F=0
# test_dir = Path("/Users/Mahsa/Desktop/test/GBM39EC")
# graph_file_path = Path("/Users/Mahsa/Desktop/test/GBM39EC/GBM39EC_EGFR_graph.txt")


# -----------------------------------------
# Parse command-line arguments
# -----------------------------------------
parser = argparse.ArgumentParser(description="Run Circular Element ILP Solver.")
parser.add_argument("--graph", required=True, help="Path to input graph file")
parser.add_argument("--enforce-connectivity",action="store_true",help="Enable connectivity-enforced ILP (default = off)")
parser.add_argument("--gamma",type=float,default=0.01,help="Gamma value (default = 0.01)")
parser.add_argument("--version",action="version",version=f"CE version {__version__}",help="Print version and exit")
parser.add_argument("--sort-by",choices=["CopyNumber", "LWCN"],default="CopyNumber",help="How to sort cycles/paths: CopyNumber (default) or LWCN")
parser.add_argument("--s-t-strategy",choices=["all_nodes", "intervals"],default="all_nodes",help="s/t connection strategy: 'all_nodes' (default) or 'intervals' (only interval start/end)")
args = parser.parse_args()

# -----------------------------------------
# Resolve paths
# -----------------------------------------
graph_file_path = Path(args.graph).resolve()
test_dir = graph_file_path.parent   # directory containing the graph
ilp_test_dir = test_dir
ENFORCE_CONNECTIVITY = args.enforce_connectivity  # <--- HERE
print(f"Using graph file: {graph_file_path}")
print(f"Test/output directory: {test_dir}")

#gamma=0.01
gamma = args.gamma
print(f"Gamma value set to: {gamma}")
print(f"CE version: {__version__}")
# if args.version:
#     print(f"CE version {__version__}")
#     sys.exit(0)

ilp_test_dir=test_dir


#ILP_cycle_file_path = ilp_test_dir/"ILP_cycles_Not_Ordered.txt"
amplicon_file_name = graph_file_path.name  # e.g., "amplicon1_graph.txt"
amplicon_name = amplicon_file_name.replace("_graph.txt", "")  # e.g., "amplicon1"           
#Messages_file_path = test_dir / "Messages.txt"
#LWCN_output_path = ilp_test_dir / "LWCN.txt"
# Extract amplicon number using regex
# match = re.search(r"amplicon(\d+)", amplicon_name.lower())
# amplicon_num = match.group(1) # e.g., "1"
if ENFORCE_CONNECTIVITY:
    Cycles_file_path = ilp_test_dir / "ILP_Connectivity_cycles.txt"
else:
    Cycles_file_path = ilp_test_dir / "ILP_cycles.txt"
if ENFORCE_CONNECTIVITY:
    log_file = test_dir / "ILP_Connectivity_log.txt"
else:
    log_file = test_dir / "ILP_log.txt"

# Redirect stdout to both console and file
class Logger:
    def __init__(self, log_file):
        self.terminal = sys.stdout
        self.log = open(log_file, "w")
        self.buffer = ""

    def write(self, message):
        self.buffer += message
        while "\n" in self.buffer:
            line, self.buffer = self.buffer.split("\n", 1)
            if line.strip():
                timestamp = datetime.now().strftime("[%H:%M:%S.%f ]  ")[:-3]
                self.log.write(f"{timestamp}{line}\n")
                self.terminal.write(f"{line}\n")
            else:
                self.log.write("\n")
                self.terminal.write("\n")

    def flush(self):
        if self.buffer:
            self.log.write(self.buffer)
            self.terminal.write(self.buffer)
            self.buffer = ""
        self.log.flush()
        self.terminal.flush()
# Open file in write mode (overwrite) or append mode ('a')
sys.stdout = Logger(log_file)
####################### 
test_name = test_dir.name

ilp_test_dir.mkdir(exist_ok=True)  # Create test* subdir in ILP if needed 
print(f"CE version: {__version__}")
print(f"Gamma value set to: {gamma}")
print(f"Using graph file: {graph_file_path}")
print(f"Test/output directory: {test_dir}")

#gamma=0.01
gamma = args.gamma


if ENFORCE_CONNECTIVITY:
    print("Connectivity is enforced.")
else:
    print("Connectivity is NOT enforced.")
#print("start time is: " , time.time())
print("Creating the graph:")
# Create the graph with integer nodes for foldbacks
(
    G,
    nodes,
    concordant_edges,
    capacity_concordant_edges,
    discordant_edges,
    capacity_discordant_edges,
    read_count_discordant_edges,
    K_discordant_edges,
    sequence_edges,
    capacity_sequence_edges,
    length_sequence_edges,
    p_sequence_edges,
    p_concordant_edges,
    p_discordant_edges,
    path_constraints
) = Create_graph(graph_file_path) 
if ENFORCE_CONNECTIVITY:
    G, nodes, discordant_edges, capacity_discordant_edges, sequence_edges, capacity_sequence_edges, K_discordant_edges, length_sequence_edges= creat_s_t_graph_for_connected_cycles(G,capacity_discordant_edges,capacity_sequence_edges,K_discordant_edges,length_sequence_edges)
capacity_sequence_edges_original=capacity_sequence_edges.copy()
total_solution_weight=0
weight_graph=compute_Length_weighted_copy_number_graph(capacity_sequence_edges,length_sequence_edges)
weight_graph_excluding_some_segments=compute_LWCN_graph_excluding_some_edges(capacity_sequence_edges_original,length_sequence_edges)
p_ij = {}        
weight_ratio_ILP_cycle=[]
weight_ratio_ILP_cycle_excluding_some_segments=[]
weight_ratio_ILP_path=[]
weight_ratio_ILP_path_excluding_some_segments=[]
Iteration=1
Cycle_Or_Path_Number=1
Num_paths=0
Num_cycles=0
All_Cycles = []
All_Paths = []
Disconnected_Solutions = []
print("Start cycle extraction")
while True: #while Iteration==1:
    ###################################################
    # ########### Solve ###############################
    ###################################################
    
    ####### Cycle solution
    ILP_Cycle_model, variables= build_cycle_model(
                    G,
                    nodes,
                    concordant_edges,
                    capacity_concordant_edges,
                    discordant_edges,
                    capacity_discordant_edges,
                    read_count_discordant_edges,
                    K_discordant_edges,
                    sequence_edges,
                    capacity_sequence_edges,
                    length_sequence_edges,
                    p_sequence_edges,
                    p_concordant_edges,
                    p_discordant_edges,
                    path_constraints,
                    ENFORCE_CONNECTIVITY,
                    gamma)
    
    ILP_Cycle_model.setParam('TimeLimit', TIME_LIMIT_SECONDS)
    print(f"Start extracting cycle {Cycle_Or_Path_Number}")
    start_time = time.time()
    ILP_Cycle_model.optimize()
    end_time = time.time()
    
    X_concordant_edge = variables["X_concordant_edge"]
    X_discordant_edge= variables["X_discordant_edge"]
    X_sequence_edge = variables["X_sequence_edge"]
    f_concordant_edge = variables["f_concordant_edge"]
    f_discordant_edge = variables["f_discordant_edge"]
    f_sequence_edge = variables["f_sequence_edge"]
    F = variables["F"]
    P = variables["P"]
    y_discordant = variables["y_discordant"]
    z_discordant = variables["z_discordant"]
    
    X_concordant_edge_val  = {k: float(v.X) for k, v in X_concordant_edge.items()}
    X_discordant_edge_val  = {k: float(v.X) for k, v in X_discordant_edge.items()}
    X_sequence_edge_val    = {k: float(v.X) for k, v in X_sequence_edge.items()}
    f_concordant_edge_val  = {k: float(v.X) for k, v in f_concordant_edge.items()}
    f_discordant_edge_val  = {k: float(v.X) for k, v in f_discordant_edge.items()}
    f_sequence_edge_val    = {k: float(v.X) for k, v in f_sequence_edge.items()}
    F_val = float(F.X)
    P_val = {k: float(v.X) for k, v in P.items()}
    
    if ILP_Cycle_model.status != GRB.OPTIMAL  or   total_solution_weight > 0.9*weight_graph or F_val<1 :
        
        Num_cycles=Cycle_Or_Path_Number-1
        break
    max_random_search_heuristic_attempts = 100
    # if ENFORCE_CONNECTIVITY:
    #     ILP_output_path = test_dir / f"ILP_Connectivityoutput_iter{Iteration}.txt"
    # else: 
    #     ILP_output_path = test_dir / f"ILP_output_iter{Iteration}.txt"

    # write_ILP_output(ILP_output_path,ILP_Cycle_model,sequence_edges,f_sequence_edge_val,X_sequence_edge_val,concordant_edges,f_concordant_edge_val,
    #                  X_concordant_edge_val,discordant_edges,f_discordant_edge_val, X_discordant_edge_val,P_val,path_constraints)                                                       

    segments, Path_Constraints = parse_graph_file(graph_file_path)
    all_closed_walks, error_in_closed_walks = create_closed_walks(sequence_edges, f_sequence_edge_val, discordant_edges, f_discordant_edge_val, 
                             concordant_edges, f_concordant_edge_val,max_random_search_heuristic_attempts,segments, F_val)
    Final_closed_walk, not_merged_walks = merge_all_closed_walks(all_closed_walks)
    Final_closed_walk_test, test_message = test_closed_walk(Final_closed_walk)
    Final_closed_walk_completion_test = test_merged_closed_walk(Final_closed_walk, sequence_edges, f_sequence_edge_val,
               discordant_edges, f_discordant_edge_val,
               concordant_edges, f_concordant_edge_val, F_val,Path=False)
    delta_F=0
    if not Final_closed_walk_completion_test:
        delta_F=Increase_CN(Final_closed_walk, F_val, capacity_sequence_edges, capacity_concordant_edges, capacity_discordant_edges)
        #print("F =", F)
        #print("Cycle 1 = ", Final_closed_walk)
        #print("delta_F1 = " ,delta_F)
        #print("f_seq_edge before =" , f_sequence_edge_val)
        if delta_F>0:
            start_time_to_update = time.time()
            #print("f_sequence_before =" , f_sequence_edge_val)
            walk_edges = {
                edge for edge_type, edge in Final_closed_walk
                if edge_type in ("sequence edge", "concordant edge", "discordant edge")
            }
            f_sequence_edge_val = {
                edge: val + (val / F_val) * delta_F if edge in walk_edges else val
                for edge, val in f_sequence_edge_val.items()
            }
            f_concordant_edge_val = {
                edge: val + (val / F_val) * delta_F if edge in walk_edges else val
                for edge, val in f_concordant_edge_val.items()
            }
            f_discordant_edge_val = {
                edge: val + (val / F_val) * delta_F if edge in walk_edges else val
                for edge, val in f_discordant_edge_val.items()
            }
            # f_sequence_edge_val = {edge: val + val / F_val * delta_F
            #    for edge, val in f_sequence_edge_val.items()}
            # #print("f_sequence_after =" , f_sequence_edge_val)
            # f_concordant_edge_val = {edge: val + val / F_val * delta_F
            #      for edge, val in f_concordant_edge_val.items()}
            # f_discordant_edge_val = {edge: val + val / F_val * delta_F
            #      for edge, val in f_discordant_edge_val.items()}
            F_val += delta_F
            end_time_to_update = time.time()
            
    final_segment_str = convert_walk_to_segment_string(Final_closed_walk, segments, copy_count=F_val)
    final_segment_str = final_segment_str.replace("Cycle=1", f"Cycle={len(all_closed_walks)+1}")
    Final_Closed_Walk_segments = final_segment_str.split("Segments=")[1].split(',')
    print(f"Finish extracting cycle {Cycle_Or_Path_Number}")
    print("Metric Calculation and Print files")
    matched_pcs = match_e_constraints(Final_Closed_Walk_segments, Path_Constraints)
    matched_index_set = set(matched_pcs)
    # write_ILP_cycle_file(segments, ILP_cycle_file_path,final_segment_str, Cycle_Or_Path_Number,
    #                          Final_Closed_Walk_segments, Final_closed_walk_test, Final_closed_walk_completion_test, error_in_closed_walks,
    #                          Path_Constraints, matched_index_set,max_random_search_heuristic_attempts)
    Solution_Type="Cycle"
    # write_Messages_file(Solution_Type,Messages_file_path, Final_closed_walk_test, Final_closed_walk_completion_test, error_in_closed_walks,max_random_search_heuristic_attempts,Cycle_Or_Path_Number)

    weight_ILP_cycle=compute_Length_weighted_copy_number_solution(sequence_edges,f_sequence_edge_val,length_sequence_edges, Final_closed_walk)                
    #weight_ILP_cycle = as_float(weight_ILP_cycle_raw)   # get numeric value
    length_ILP_cycle = weight_ILP_cycle/F_val
    weight_ratio_ILP_cycle.append(weight_ILP_cycle/weight_graph)    
    weight_ILP_cycle_excluding_some_segments=compute_LWCN_solution_excludeing_some_segments(capacity_sequence_edges_original,sequence_edges,f_sequence_edge_val,length_sequence_edges, Final_closed_walk)                

    length_ILP_cycle_excluding_some_segments = weight_ILP_cycle_excluding_some_segments / F_val
    weight_ratio_ILP_cycle_excluding_some_segments.append(weight_ILP_cycle_excluding_some_segments/weight_graph_excluding_some_segments)
    total_solution_weight += weight_ILP_cycle_excluding_some_segments
    # write_LWCN_cycle_in_file(LWCN_output_path
    #                          ,weight_ILP_cycle,weight_graph,weight_ratio_ILP_cycle,length_ILP_cycle
    #                          ,weight_ILP_cycle_excluding_some_segments,weight_graph_excluding_some_segments,weight_ratio_ILP_cycle_excluding_some_segments, length_ILP_cycle_excluding_some_segments
    #                          ,F_val, Cycle_Or_Path_Number,Num_paths)
    Cycle_Info={"Cycle": Cycle_Or_Path_Number,
                "CopyNumber": F_val,
                "Segments": Final_Closed_Walk_segments,
                "Length": length_ILP_cycle,
                "Length_excluding_segments": length_ILP_cycle_excluding_some_segments,
                "LWCN": weight_ILP_cycle,  # length-weighted copy number
                "LWCN_excluding_segments": weight_ILP_cycle_excluding_some_segments,
                "LWCNR": weight_ILP_cycle/weight_graph,  # length-weighted copy number
                "LWCNR_excluding_segments": weight_ILP_cycle_excluding_some_segments/weight_graph_excluding_some_segments,
                "Matched_Index_Set": matched_index_set, # This is for the path constraints satisifed
                }
    if not Final_closed_walk_completion_test :
      Disconnected_Solutions.append(Cycle_Info)
      if delta_F>0:
         Disconnected_Solutions[-1]["Delta_F"] = delta_F
         Disconnected_Solutions[-1]["CopyNumber_After_Increase"] = F_val
         Disconnected_Solutions[-1]["Increase_Time"] = end_time_to_update - start_time_to_update
      else:
         Disconnected_Solutions[-1]["Delta_F"] = 0.0
         Disconnected_Solutions[-1]["CopyNumber_After_Increase"] = F_val
         Disconnected_Solutions[-1]["Increase_Time"] = 0.0
    All_Cycles.append(Cycle_Info) 
    
    while not_merged_walks != []: # This while loop is for situation of disconnected cycle where they do not get merged
        Cycle_Or_Path_Number +=1
        all_closed_walks = not_merged_walks
        Final_closed_walk, not_merged_walks = merge_all_closed_walks(all_closed_walks)
        Final_closed_walk_test, test_message = test_closed_walk(Final_closed_walk)
        Final_closed_walk_completion_test = test_merged_closed_walk(Final_closed_walk, sequence_edges, f_sequence_edge_val,
                   discordant_edges, f_discordant_edge_val,
                   concordant_edges, f_concordant_edge_val, F_val, Path=False)
        F_val -= delta_F
        delta_F=0
        if not Final_closed_walk_completion_test:
            #Disconnected_Solutions.append(Cycle_Info) 
            delta_F=Increase_CN(Final_closed_walk, F_val, capacity_sequence_edges, capacity_concordant_edges, capacity_discordant_edges)
            # print("F=", F)
            # print("Cycle 2 = ", Final_closed_walk)
            # print("delta_F2 = " ,delta_F)
            # print("f_sequence_before =" , f_sequence_edge_val)
            if delta_F>0:
                start_time_to_update = time.time()
                walk_edges = {
                    edge for edge_type, edge in Final_closed_walk
                    if edge_type in ("sequence edge", "concordant edge", "discordant edge")
                }
                f_sequence_edge_val = {
                    edge: val + (val / F_val) * delta_F if edge in walk_edges else val
                    for edge, val in f_sequence_edge_val.items()
                }
                f_concordant_edge_val = {
                    edge: val + (val / F_val) * delta_F if edge in walk_edges else val
                    for edge, val in f_concordant_edge_val.items()
                }
                f_discordant_edge_val = {
                    edge: val + (val / F_val) * delta_F if edge in walk_edges else val
                    for edge, val in f_discordant_edge_val.items()
                }
                # f_sequence_edge_val = {edge: val + val / F_val * delta_F
                #    for edge, val in f_sequence_edge_val.items()}
                # print("f_sequence_after =" , f_sequence_edge_val)
                # f_concordant_edge_val = {edge: val + val / F_val * delta_F
                #      for edge, val in f_concordant_edge_val.items()}
                # f_discordant_edge_val = {edge: val + val / F_val * delta_F
                #      for edge, val in f_discordant_edge_val.items()}
                F_val += delta_F
                end_time_to_update = time.time()
        # Write Final_closed_walk as list of segments
        final_segment_str = convert_walk_to_segment_string(Final_closed_walk, segments, copy_count=F_val)
        final_segment_str = final_segment_str.replace("Cycle=1", f"Cycle={len(all_closed_walks)+1}")
        Final_Closed_Walk_segments = final_segment_str.split("Segments=")[1].split(',')
        #matched_index_set = set()
        matched_pcs = match_e_constraints(Final_Closed_Walk_segments, Path_Constraints)
        matched_index_set = set(matched_pcs)
        # write_ILP_cycle_file(segments, ILP_cycle_file_path,final_segment_str, Cycle_Or_Path_Number,
        #                          Final_Closed_Walk_segments, Final_closed_walk_test, Final_closed_walk_completion_test, error_in_closed_walks,
        #                          Path_Constraints, matched_index_set,max_random_search_heuristic_attempts)
        Solution_Type="Cycle"
        # write_Messages_file(Solution_Type,Messages_file_path, Final_closed_walk_test, Final_closed_walk_completion_test, error_in_closed_walks,max_random_search_heuristic_attempts,Cycle_Or_Path_Number)
        weight_ILP_cycle=compute_Length_weighted_copy_number_solution(sequence_edges,f_sequence_edge_val,length_sequence_edges, Final_closed_walk)                
        #weight_ILP_cycle = as_float(weight_ILP_cycle_raw)   # get numeric value
        length_ILP_cycle = weight_ILP_cycle/F_val
        weight_ratio_ILP_cycle.append(weight_ILP_cycle/weight_graph)    
        weight_ILP_cycle_excluding_some_segments=compute_LWCN_solution_excludeing_some_segments(capacity_sequence_edges_original,sequence_edges,f_sequence_edge_val,length_sequence_edges, Final_closed_walk)                

        length_ILP_cycle_excluding_some_segments = weight_ILP_cycle_excluding_some_segments / F_val
        weight_ratio_ILP_cycle_excluding_some_segments.append(weight_ILP_cycle_excluding_some_segments/weight_graph_excluding_some_segments)
        total_solution_weight += weight_ILP_cycle_excluding_some_segments
        # write_LWCN_cycle_in_file(LWCN_output_path
        #                          ,weight_ILP_cycle,weight_graph,weight_ratio_ILP_cycle,length_ILP_cycle
        #                          ,weight_ILP_cycle_excluding_some_segments,weight_graph_excluding_some_segments,weight_ratio_ILP_cycle_excluding_some_segments, length_ILP_cycle_excluding_some_segments
        #                          ,F_val, Cycle_Or_Path_Number,Num_paths)
        Cycle_Info={"Cycle": Cycle_Or_Path_Number,
                    "CopyNumber": F_val,
                    "Segments": Final_Closed_Walk_segments,
                    "Length": length_ILP_cycle,
                    "Length_excluding_segments": length_ILP_cycle_excluding_some_segments,
                    "LWCN": weight_ILP_cycle,  # length-weighted copy number
                    "LWCN_excluding_segments": weight_ILP_cycle_excluding_some_segments,
                    "LWCNR": weight_ILP_cycle/weight_graph,  # length-weighted copy number
                    "LWCNR_excluding_segments": weight_ILP_cycle_excluding_some_segments/weight_graph_excluding_some_segments,
                    "Matched_Index_Set": matched_index_set, # This is for the path constraints satisifed
                    }
        if not Final_closed_walk_completion_test :
             Disconnected_Solutions.append(Cycle_Info)
             if delta_F>0:
                Disconnected_Solutions[-1]["Delta_F"] = delta_F
                Disconnected_Solutions[-1]["CopyNumber_After_Increase"] = F_val
                Disconnected_Solutions[-1]["Increase_Time"] = end_time_to_update - start_time_to_update
             else:
                Disconnected_Solutions[-1]["Delta_F"] = 0.0
                Disconnected_Solutions[-1]["CopyNumber_After_Increase"] = F_val
                Disconnected_Solutions[-1]["Increase_Time"] = 0.0
        All_Cycles.append(Cycle_Info)
        
       
    capacity_concordant_edges, capacity_discordant_edges, capacity_sequence_edges = update_the_graph(concordant_edges, f_concordant_edge_val, X_concordant_edge_val, capacity_concordant_edges,
                             discordant_edges, f_discordant_edge_val, X_discordant_edge_val, capacity_discordant_edges,
                             sequence_edges, f_sequence_edge_val, X_sequence_edge_val, capacity_sequence_edges)

    Iteration += 1 
    Cycle_Or_Path_Number +=1

######## Path solution
if ENFORCE_CONNECTIVITY:
    (
        G,
        nodes,
        discordant_edges,
        capacity_discordant_edges,
        sequence_edges,
        capacity_sequence_edges,
        K_discordant_edges,
        length_sequence_edges
    ) = remove_s_t_from_graph(
        G,
        capacity_discordant_edges,
        capacity_sequence_edges,
        K_discordant_edges,
        length_sequence_edges
    )
        
#G,G.nodes,discordant_edges,capacity_discordant_edges,K_discordant_edges = creat_s_t_graph(G,capacity_discordant_edges,K_discordant_edges)
#G,G.nodes,discordant_edges,capacity_discordant_edges,K_discordant_edges = creat_s_t_graph_Coral(G,capacity_discordant_edges,K_discordant_edges, sequence_edges)

if args.s_t_strategy == "all_nodes":
    G, nodes, discordant_edges, capacity_discordant_edges, K_discordant_edges = \
        creat_s_t_graph(G, capacity_discordant_edges, K_discordant_edges)
else:
    sequence_edges_left = {
        edge for edge, cap in capacity_sequence_edges.items()
        if cap > 0
    }
    G, nodes, discordant_edges, capacity_discordant_edges, K_discordant_edges = \
        creat_s_t_graph_Coral(G, capacity_discordant_edges, K_discordant_edges, sequence_edges_left)

print("Start path extraction")    
start=time.time()
delta_F=0
while True: # while Iteration==1:
    ILP_Path_model, variables= build_path_model(
                    G,
                    nodes,
                    concordant_edges,
                    capacity_concordant_edges,
                    discordant_edges,
                    capacity_discordant_edges,
                    read_count_discordant_edges,
                    K_discordant_edges,
                    sequence_edges,
                    capacity_sequence_edges,
                    length_sequence_edges,
                    p_sequence_edges,
                    p_concordant_edges,
                    p_discordant_edges,
                    path_constraints,
                    ENFORCE_CONNECTIVITY,
                    gamma) 
    ILP_Path_model.setParam('TimeLimit', TIME_LIMIT_SECONDS)
    print(f"Start running Iteration {Iteration}")
    start_time = time.time()
    ILP_Path_model.optimize()
    end_time = time.time()
    print(f"Finish running Iteration {Iteration}")

    X_concordant_edge = variables["X_concordant_edge"]
    X_discordant_edge= variables["X_discordant_edge"]
    X_sequence_edge = variables["X_sequence_edge"]
    f_concordant_edge = variables["f_concordant_edge"]
    f_discordant_edge = variables["f_discordant_edge"]
    f_sequence_edge = variables["f_sequence_edge"]
    F = variables["F"]
    P = variables["P"]
    y_discordant = variables["y_discordant"]
    z_discordant = variables["z_discordant"]
    
    X_concordant_edge_val  = {k: float(v.X) for k, v in X_concordant_edge.items()}
    X_discordant_edge_val  = {k: float(v.X) for k, v in X_discordant_edge.items()}
    X_sequence_edge_val    = {k: float(v.X) for k, v in X_sequence_edge.items()}
    f_concordant_edge_val  = {k: float(v.X) for k, v in f_concordant_edge.items()}
    f_discordant_edge_val  = {k: float(v.X) for k, v in f_discordant_edge.items()}
    f_sequence_edge_val    = {k: float(v.X) for k, v in f_sequence_edge.items()}
    F_val = float(F.X)
    P_val = {k: float(v.X) for k, v in P.items()}
    
    if ILP_Path_model.status != GRB.OPTIMAL or total_solution_weight > 0.9*weight_graph or F_val<1 :
        Num_paths=Cycle_Or_Path_Number-Num_cycles-1
        break
    max_random_search_heuristic_attempts = 100       
    # if ENFORCE_CONNECTIVITY:
    #     ILP_output_path = ilp_test_dir / f"ILP_Connectivity_output_iter{Iteration}.txt"
    # else:
    #     ILP_output_path = ilp_test_dir / f"ILP_output_iter{Iteration}.txt"

    # write_ILP_output(ILP_output_path,ILP_Cycle_model,sequence_edges,f_sequence_edge_val,X_sequence_edge_val,concordant_edges,f_concordant_edge_val,
    #                  X_concordant_edge_val,discordant_edges,f_discordant_edge_val, X_discordant_edge_val,P_val,path_constraints)                                                       
    
    segments, Path_Constraints = parse_graph_file(graph_file_path)
    s_t_Walk, error_in_closed_walks = create_s_t_walks(sequence_edges, f_sequence_edge_val, discordant_edges, f_discordant_edge_val, 
                         concordant_edges, f_concordant_edge_val,max_random_search_heuristic_attempts,segments, F_val)

    merged_s_t_walk, not_merged_closed_walks_to_s_t_path = merge_all_closed_walks(s_t_Walk)
    final_segment_str = convert_s_t_walk_to_segment_string(s_t_Walk[0][2:-2], segments, copy_count=F_val)
    final_segment_str = re.sub(r"Path=\d+", f"Path={Cycle_Or_Path_Number}", final_segment_str)
    parts = final_segment_str.split(';', 1)  # split at first semicolon
    Final_s_t_Walk = final_segment_str.split("Segments=")[1].split(',')
    Final_s_t_Walk_segments = final_segment_str.split("Segments=")[1].split(',')
    #matched_index_set = set()
    matched_pcs = match_e_constraints(Final_s_t_Walk_segments, Path_Constraints)
    matched_index_set = set(matched_pcs)
    Final_s_t_walk_completion_test = test_merged_closed_walk(s_t_Walk[0], sequence_edges, f_sequence_edge_val,
               discordant_edges, f_discordant_edge_val,
               concordant_edges, f_concordant_edge_val, F_val,Path=True)
    F_val -= delta_F
    delta_F=0
    if not Final_s_t_walk_completion_test:
        delta_F=Increase_CN(s_t_Walk[0], F_val, capacity_sequence_edges, capacity_concordant_edges, capacity_discordant_edges)
        # print("F=", F)
        # print("Cpath 3 = ", Final_closed_walk)
        # print("delta_F3 = " ,delta_F)
        # print("f_sequence_before =" , f_sequence_edge_val)
        if delta_F>0:
            start_time_to_update = time.time()
            walk_edges = {
                edge for edge_type, edge in Final_closed_walk
                if edge_type in ("sequence edge", "concordant edge", "discordant edge")
            }
            f_sequence_edge_val = {
                edge: val + (val / F_val) * delta_F if edge in walk_edges else val
                for edge, val in f_sequence_edge_val.items()
            }
            f_concordant_edge_val = {
                edge: val + (val / F_val) * delta_F if edge in walk_edges else val
                for edge, val in f_concordant_edge_val.items()
            }
            f_discordant_edge_val = {
                edge: val + (val / F_val) * delta_F if edge in walk_edges else val
                for edge, val in f_discordant_edge_val.items()
            }
            # f_sequence_edge_val = {edge: val + val / F_val * delta_F
            #    for edge, val in f_sequence_edge_val.items()}
            # print("f_sequence_after =" , f_sequence_edge_val)
            # f_concordant_edge_val = {edge: val + val / F_val * delta_F
            #      for edge, val in f_concordant_edge_val.items()}
            # f_discordant_edge_val = {edge: val + val / F_val * delta_F
            #      for edge, val in f_discordant_edge_val.items()}
            F_val += delta_F
            end_time_to_update = time.time() 
    # write_paths_in_ILP_cycle_file(original_cycle_file, new_output_file,final_segment_str, 
    #                          Final_s_t_Walk_segments, Path_Constraints, matched_index_set,Cycle_Or_Path_Number,Final_s_t_walk_completion_test)
    # write_paths_in_ILP_cycle_file(segments, ILP_cycle_file_path, final_segment_str, Iteration,Path_Constraints, matched_index_set)
    Final_closed_walk_test=True
    Solution_Type="Path"
    # write_Messages_file(Solution_Type,Messages_file_path, Final_closed_walk_test, Final_s_t_walk_completion_test, error_in_closed_walks,max_random_search_heuristic_attempts,Cycle_Or_Path_Number)
    weight_ILP_path =compute_Length_weighted_copy_number_solution(sequence_edges,f_sequence_edge_val,length_sequence_edges, s_t_Walk[0])                 
    length_ILP_path = weight_ILP_path/F_val    
    weight_ratio_ILP_path.append(weight_ILP_path/weight_graph)
    weight_ILP_path_excluding_some_segments=compute_LWCN_solution_excludeing_some_segments(capacity_sequence_edges_original,sequence_edges,f_sequence_edge_val,length_sequence_edges, s_t_Walk[0])
    length_ILP_path_excluding_some_segments = weight_ILP_path_excluding_some_segments /F_val
    weight_ratio_ILP_path_excluding_some_segments.append(weight_ILP_path_excluding_some_segments/weight_graph_excluding_some_segments)
    total_solution_weight += weight_ILP_path_excluding_some_segments
    # write_LWCN_path_in_file(LWCN_output_path
    #                          ,weight_ILP_path,weight_graph,weight_ratio_ILP_path,length_ILP_path
    #                          ,weight_ILP_path_excluding_some_segments,weight_graph_excluding_some_segments,weight_ratio_ILP_path_excluding_some_segments, length_ILP_path_excluding_some_segments
    #                          ,F_val,Cycle_Or_Path_Number,Num_cycles)
    Path_Info={"Path": Cycle_Or_Path_Number,
               "CopyNumber": F_val,
                "Segments": Final_s_t_Walk_segments,
                "Length": length_ILP_path,
                "Length_excluding_segments": length_ILP_path_excluding_some_segments,
                "LWCN": weight_ILP_path,  # length-weighted copy number
                "LWCN_excluding_segments": weight_ILP_path_excluding_some_segments,
                "LWCNR": weight_ILP_path/weight_graph,  # length-weighted copy number
                "LWCNR_excluding_segments": weight_ILP_path_excluding_some_segments/weight_graph_excluding_some_segments,
                "Matched_Index_Set": matched_index_set, # This is for the path constraints satisifed
                }
    All_Paths.append(Path_Info) 
    
    # if not Final_closed_walk_completion_test :
    #     Disconnected_Solutions.append(Path_Info) 
    
    Num_paths+=1
    Cycle_Or_Path_Number +=1
    while not_merged_closed_walks_to_s_t_path != []: # This while loop is for situation of disconnected cycle where they do not get merge to the s-t path because they are disconnected and there is not shared node or segment between the s-t path and the cycle
        #### Whatever we get in this loop is a cycle because we already got the s-t path from the solution and saved and printed it 
        Disconnected_Solutions.append(Path_Info) 
        
        all_closed_walks = not_merged_closed_walks_to_s_t_path
        Final_closed_walk, not_merged_closed_walks_to_s_t_path = merge_all_closed_walks(all_closed_walks)
        Final_closed_walk_test, test_message = test_closed_walk(Final_closed_walk)
        Final_closed_walk_completion_test = test_merged_closed_walk(Final_closed_walk, sequence_edges, f_sequence_edge_val,
                   discordant_edges, f_discordant_edge_val,
                   concordant_edges, f_concordant_edge_val, F_val, Path=False)
        F_val -= delta_F
        delta_F=0
        if not Final_closed_walk_completion_test:
            #Disconnected_Solutions.append(Cycle_Info)
            delta_F=Increase_CN(Final_closed_walk, F_val, capacity_sequence_edges, capacity_concordant_edges, capacity_discordant_edges)
            # print("F=", F)
            # print("Cycle 4 = ", Final_closed_walk)
            # print("delta_F4 = " ,delta_F)
            # print("f_sequence_before =" , f_sequence_edge_val)
            if delta_F>0:
                start_time_to_update = time.time()
                walk_edges = {
                    edge for edge_type, edge in Final_closed_walk
                    if edge_type in ("sequence edge", "concordant edge", "discordant edge")
                }
                f_sequence_edge_val = {
                    edge: val + (val / F_val) * delta_F if edge in walk_edges else val
                    for edge, val in f_sequence_edge_val.items()
                }
                f_concordant_edge_val = {
                    edge: val + (val / F_val) * delta_F if edge in walk_edges else val
                    for edge, val in f_concordant_edge_val.items()
                }
                f_discordant_edge_val = {
                    edge: val + (val / F_val) * delta_F if edge in walk_edges else val
                    for edge, val in f_discordant_edge_val.items()
                }
                # f_sequence_edge_val = {edge: val + val / F_val * delta_F
                #    for edge, val in f_sequence_edge_val.items()}
                # print("f_sequence_after =" , f_sequence_edge_val)
                # f_concordant_edge_val = {edge: val + val / F_val * delta_F
                #      for edge, val in f_concordant_edge_val.items()}
                # f_discordant_edge_val = {edge: val + val / F_val * delta_F
                #      for edge, val in f_discordant_edge_val.items()}
                F_val += delta_F
                end_time_to_update = time.time()
                
        # Write Final_closed_walk as list of segments
        final_segment_str = convert_walk_to_segment_string(Final_closed_walk, segments, copy_count=F_val)
        final_segment_str = final_segment_str.replace("Cycle=1", f"Cycle={len(all_closed_walks)+1}")
        Final_Closed_Walk_segments = final_segment_str.split("Segments=")[1].split(',')
        #matched_index_set = set()
        matched_pcs = match_e_constraints(Final_Closed_Walk_segments, Path_Constraints)
        matched_index_set = set(matched_pcs)
        # write_ILP_cycle_file(segments, ILP_cycle_file_path,final_segment_str, Cycle_Or_Path_Number,
        #                          Final_Closed_Walk_segments, Final_closed_walk_test, Final_closed_walk_completion_test, error_in_closed_walks,
        #                          Path_Constraints, matched_index_set,max_random_search_heuristic_attempts)
        Solution_Type="Cycle"
        # write_Messages_file(Solution_Type,Messages_file_path, Final_closed_walk_test, Final_closed_walk_completion_test, error_in_closed_walks,max_random_search_heuristic_attempts, Cycle_Or_Path_Number)
        weight_ILP_cycle=compute_Length_weighted_copy_number_solution(sequence_edges,f_sequence_edge_val,length_sequence_edges, Final_closed_walk)                

        length_ILP_cycle = weight_ILP_cycle/F_val
        weight_ratio_ILP_cycle.append(weight_ILP_cycle/weight_graph)
        weight_ILP_cycle_excluding_some_segments=compute_LWCN_solution_excludeing_some_segments(capacity_sequence_edges_original,sequence_edges,f_sequence_edge_val,length_sequence_edges, Final_closed_walk)                

        length_ILP_cycle_excluding_some_segments = weight_ILP_cycle_excluding_some_segments/F_val
        weight_ratio_ILP_cycle_excluding_some_segments.append(weight_ILP_cycle_excluding_some_segments/weight_graph_excluding_some_segments)

        total_solution_weight += weight_ILP_cycle_excluding_some_segments
        # write_LWCN_cycle_in_file(LWCN_output_path
        #                          ,weight_ILP_cycle,weight_graph,weight_ratio_ILP_cycle,length_ILP_cycle
        #                          ,weight_ILP_cycle_excluding_some_segments,weight_graph_excluding_some_segments,weight_ratio_ILP_cycle_excluding_some_segments, length_ILP_cycle_excluding_some_segments
        #                          ,F_val, Cycle_Or_Path_Number,Num_paths)
        Cycle_Info={"Cycle": Cycle_Or_Path_Number,
                    "CopyNumber": F_val,
                    "Segments": Final_Closed_Walk_segments,
                    "Length": length_ILP_cycle,
                    "Length_excluding_segments": length_ILP_cycle_excluding_some_segments,
                    "LWCN": weight_ILP_cycle,  # length-weighted copy number
                    "LWCN_excluding_segments": weight_ILP_cycle_excluding_some_segments,
                    "LWCNR": weight_ILP_cycle/weight_graph,  # length-weighted copy number
                    "LWCNR_excluding_segments": weight_ILP_cycle_excluding_some_segments/weight_graph_excluding_some_segments,
                    "Matched_Index_Set": matched_index_set, # This is for the path constraints satisifed
                    }
        if not Final_closed_walk_completion_test :
            Disconnected_Solutions.append(Cycle_Info)
            if delta_F>0:
               Disconnected_Solutions[-1]["Delta_F"] = delta_F
               Disconnected_Solutions[-1]["CopyNumber_After_Increase"] = F_val
               Disconnected_Solutions[-1]["Increase_Time"] = end_time_to_update - start_time_to_update
            else:
               Disconnected_Solutions[-1]["Delta_F"] = 0.0
               Disconnected_Solutions[-1]["CopyNumber_After_Increase"] = F_val
               Disconnected_Solutions[-1]["Increase_Time"] = 0.0
        All_Cycles.append(Cycle_Info) 
        Cycle_Or_Path_Number +=1
        Num_cycles += 1
        
    capacity_concordant_edges, capacity_discordant_edges, capacity_sequence_edges = update_the_graph(concordant_edges, f_concordant_edge_val, X_concordant_edge_val, capacity_concordant_edges,
                         discordant_edges, f_discordant_edge_val, X_discordant_edge_val, capacity_discordant_edges,
                         sequence_edges, f_sequence_edge_val, X_sequence_edge_val, capacity_sequence_edges)
    Iteration += 1 
# All_Cycles_Ordered = sorted(All_Cycles, key=lambda x: x["LWCN_excluding_segments"], reverse=True)
# All_Paths_Ordered = sorted(All_Paths, key=lambda x: x["LWCN_excluding_segments"], reverse=True)
# All_Cycles_Ordered = sorted(All_Cycles, key=lambda x: x["CopyNumber"], reverse=True)
# All_Paths_Ordered = sorted(All_Paths, key=lambda x: x["CopyNumber"], reverse=True)

#Determine sorting key
if args.sort_by == "CopyNumber":
    cycle_sort_key = lambda x: x["CopyNumber"]
else:  # LWCN
    cycle_sort_key = lambda x: x["LWCN_excluding_segments"]

# Sort cycles and paths
All_Cycles_Ordered = sorted(All_Cycles, key=cycle_sort_key, reverse=True)
All_Paths_Ordered = sorted(All_Paths, key=cycle_sort_key, reverse=True)

for i, cycle in enumerate(All_Cycles_Ordered, start=1):
    cycle["New_Order_Of_Cycle_Number"] = i
for i, path in enumerate(All_Paths_Ordered, start=1):
    path["New_Order_Of_Path_Number"] = i+Num_cycles

# Add new order numbers to disconnected solutions (keep the same dicts)
for item in Disconnected_Solutions:
    if "Cycle" in item:
        # Find the cycle in All_Cycles_Ordered that matches the original number
        for cycle in All_Cycles_Ordered:
            if cycle["Cycle"] == item["Cycle"]:
                item["New_Order_Of_Cycle_Number"] = cycle["New_Order_Of_Cycle_Number"]
                break
    elif "Path" in item:
        for path in All_Paths_Ordered:
            if path["Path"] == item["Path"]:
                item["New_Order_Of_Path_Number"] = path["New_Order_Of_Path_Number"]
                break
            
write_all_cycles_and_paths(segments, Cycles_file_path, All_Cycles_Ordered, All_Paths_Ordered, Path_Constraints)
    
    
