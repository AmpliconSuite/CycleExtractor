# CE — Cycle Extractor
**CE (Cycle Extractor)**

CE (Cycle Extractor) is a tool that takes a breakpoint graph (generated from short-read or long-read sequencing) and identifies candidate ecDNA structures by extracting cycles and paths. CE runs on most modern Unix-like operating systems, including: .........

Requirements: 

python>=..... 
Gurobi optimizer licence (free for academic use at: https://www.gurobi.com/).

Usage
Basic Run:

python3 CE.py --graph <path_to_graph_file>

Optional: Enforce connectivity

By default, CE does not enforce connectivity constraints (this keeps the model faster).
To enforce connectivity:

python3 CE.py --graph <path_to_graph_file> --enforce-connectivity
