# CE — Cycle Extractor
**CE (Cycle Extractor)** is a tool that takes a breakpoint graph (either generated from short-read or long-read sequencing) and identifies candidate ecDNA structures by extracting cycles and paths. CE runs on most operating systems, including: Linux and macOS. 

# Requirements: 
**python>=3.10** 

**Gurobi optimizer** (free academic license): https://www.gurobi.com/ Python packages: pip install numpy networkx (gurobipy is installed automatically when you install Gurobi.) 

# Usage
**Basic Run:** python3 CE.py --graph <path_to_graph_file> --output <path_to_output(cycles)_file> 

# Optional Parameters: 

**1. Set Gamma Value** (default = 0.01)

python3 CE.py --graph sample_graph.txt --gamma <value>

**2. Check Version**
python3 CE.py --version

**3. Sort Cycles/Paths** (default = CopyNumber)
By default, CE sorts and prints cycles by their copy number. Here we give the option to the user to sort by the Length Weighted Copy Number of the cycle:

python3 CE.py --graph file.txt --sort-by LWCN

**4. s/t Node Connection Strategy** (default = all_nodes)
Choose whether to add source/sink nodes to all nodes or only to interval start/end nodes:

python3 CE.py --graph mygraph.txt --s-t-strategy intervals

**5. Enforce Connectivity** 
By default, CE does not enforce connectivity constraints (faster). To enforce connectivity:

python3 CE.py --graph <path_to_graph_file> --output <path_to_output_file> --enforce-connectivity


# Breakpoint Graph file and cycles file

The input of CE is the breakpoint graph which can be generate either from Amplicon Architech (if data is short read sequencing) or from CoRAL (if data is long read sequencing). The straucture of the graph that is accepted to the current code of CE the structure of the breakpoint graph of CoRAL. 
The cycles files of AA have different segments of AA graph file, but the cycles files of CoRAL have the same segments of CoRAL graph file. IT should be noted that CE cycles files have the same segments as the input graph file (euther from AA or CoRAL). The user should pay attention to this fact when comapring AA cycles files and CE cycles file.

# Starting from bam files

We understand that the user may start from the bam files which is not th input of CE. In that case if the user has long read data, they should follow the steps of CoRAL to get the breakpoint graph https://github.com/AmpliconSuite/CoRAL and if they have short read data, they should follow the Amplicon Architech https://github.com/virajbdeshpande/AmpliconArchitect .

