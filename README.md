# CE — Cycle Extractor
**CE (Cycle Extractor)** is a tool that takes a breakpoint graph (generated from short-read or long-read sequencing) and identifies candidate ecDNA structures by extracting cycles and paths. CE runs on most operating systems, including: Linux and macOS. 

# Requirements: 
**python>=3.10** 

**Gurobi optimizer** (free academic license): https://www.gurobi.com/ Python packages: pip install numpy networkx (gurobipy is installed automatically when you install Gurobi.) 

# Usage
**Basic Run:** python3 CE.py --graph <path_to_graph_file> --output <path_to_output(cycles)_file>

**Optional: Enforce connectivity** By default, CE does not enforce connectivity constraints (this keeps the model faster). To enforce connectivity: python3 CE.py --graph <path_to_graph_file> --output <path_to_output(cycles)_file> --enforce-connectivity
