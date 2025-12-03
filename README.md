# CE — Cycle Extractor
**CE (Cycle Extractor)** is a tool that takes a breakpoint graph (generated from short-read or long-read sequencing) and identifies candidate ecDNA structures by extracting cycles and paths. CE runs on most operating systems, including: Linux and macOS. 

# Requirements: 
**python>=3.10** 

**Gurobi optimizer** (free academic license): https://www.gurobi.com/ Python packages: pip install numpy networkx (gurobipy is installed automatically when you install Gurobi.) 

# Usage
**Basic Run:** python3 CE.py --graph <path_to_graph_file> --output <path_to_output(cycles)_file>

# Optional Parameters: 

**1. Set Gamma Value** (default = 0.01) 
python3 CE.py --graph sample_graph.txt --gamma <The value for gamma>

**2. Check Version**
python3 CE.py --version

**3. Sort Cycles/Paths** (default = CopyNumber)
Sort extracted cycles/paths by either CopyNumber or LWCN_excluding_segments:
python3 CE.py --graph file.txt --sort-by LWCN

**4. s/t Node Connection Strategy** (default = all_nodes)
Choose whether to add source/sink nodes to all nodes or only to interval start/end nodes:
python3 CE.py --graph mygraph.txt --s-t-strategy intervals

**5. Enforce Connectivity** 
By default, CE does not enforce connectivity constraints (faster). To enforce connectivity:
python3 CE.py --graph <path_to_graph_file> --output <path_to_output_file> --enforce-connectivity
