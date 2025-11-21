# CE — Cycle Extractor
**CE (Cycle Extractor)** is a computational tool for identifying candidate ecDNA structures from breakpoint graphs.
It accepts graphs generated from short-read or long-read sequencing and extracts cycles and paths that represent possible circular amplicons.


# Requirements: 

**python>=3.10**
**Gurobi optimizer** (free academic license): https://www.gurobi.com/
Python packages:
pip install numpy networkx
(Note:gurobipy is installed automatically when you install Gurobi.)

# Usage
**Basic Run:**

python3 CE.py --graph <path_to_graph_file>

**Optional: Enforce connectivity**

By default, CE does not enforce connectivity constraints (this keeps the model faster).
To enforce connectivity:

python3 CE.py --graph <path_to_graph_file> --enforce-connectivity
