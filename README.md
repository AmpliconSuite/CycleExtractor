# Cycle-Extractor (CE)
**Cycle-Extractor (CE)** is a fast and accurate tool for reconstructing extrachromosomal DNA (ecDNA) sequences from *amplicon graphs* (edge-weighted breakpoint graph). Given an amplicon graph derived from short-read (e.g., Illumina) or long-read (e.g., ONT, PacBio) sequencing data, CE identifies candidate ecDNA structures by extracting the cycles (and optionally paths) that best explain the amplified region. CE runs on most operating systems, including: Linux and macOS. 

If you use CE in your research, please cite:

- Faizrahnemoon, M., Luebeck, J., Hung, K. L., Rao, S., Prasad, G., Wong, I. T. L., Jones, M. G., Mischel, P. S., Chang, H. Y., Zhu, K., & Bafna, V. Fast and accurate resolution of ecDNA sequence using Cycle-Extractor. bioRxiv (2026): 2026-03. [https://doi.org/10.64898/2026.03.10.710955](https://doi.org/10.64898/2026.03.10.710955)
  
# Requirements
- **Python >= 3.10**
- **Gurobi Optimizer** — a free academic license is available [here](https://www.gurobi.com/). The ```gurobipy``` package is installed automatically when you install Gurobi.
- Python packages:
  - numpy
  - networkx

These dependencies can be installed with:
```bash
pip install numpy networkx
```

# Installation
Clone the repository:
```bash
git clone https://github.com/AmpliconSuite/CycleExtractor.git
cd CycleExtractor
```
Then Python script `CE.py` can be run directly (For now, build step is not required).

# Quick Start
Run CE on an amplicon graph (```*.graph```) file:
```bash
python3 CE.py --graph <path_to_graph_file> --output <path_to_output_file>
```
A sample graph is provided in the `sample_data/` directory, so you can test your installation with something like:
```bash
python3 CE.py --graph sample_data/<sample_graph.txt> --output sample_cycles.txt
```

# Input Format
CE takes the standard [AmpliconArchitect](https://github.com/AmpliconSuite/AmpliconArchitect) (paired-end short read) or [CoRAL](https://github.com/AmpliconSuite/CoRAL) (long read) amplicon graph (```*.graph```) as input. The ```*.graph``` file is plain text and contains three sections:
- `SequenceEdge:` - unbroken genomic segments, with predicted copy numbers, coverage, length, and read counts.
- `BreakpointEdge:` - `concordant` edges (adjacent segments on the same chromosome) and `discordant` edges (novel breakpoint junctions), each with a predicted copy number.
- `PathConstraint:` (optional, only from CoRAL output) - long-read derived subwalk constraints used to disambiguate cycle traversals.

An example short read graph:
```
SequenceEdge: StartPosition, EndPosition, PredictedCopyCount, AverageCoverage, Size, NumberReadsMapped
sequence	chr7:54659673-	chr7:54763279+	3.9605242507835285	9.77902004356503	103606	13333
sequence	chr7:54763280-	chr7:55127267+	118.29946643491265	299.6266250528023	363987	1435195
sequence	chr7:55127268-	chr7:55155019+	3.2711053904267375	8.351454287880353	27751	3050
sequence	chr7:55155020-	chr7:55669499+	118.29946643491255	300.26766128934145	514479	2032921
sequence	chr7:55669500-	chr7:55760699+	102.34823701028765	258.56853962048655	91199	310323
sequence	chr7:55760700-	chr7:56049369+	118.6521371601259	297.4308907534588	288669	1129877
sequence	chr7:56049370-	chr7:56114999+	4.313194975996781	11.322663422843723	65629	9779
sequence	chr7:56115000-	chr7:56149666+	3.1992170423544346	8.294523139918248	34666	3784
BreakpointEdge: StartPosition->EndPosition, PredictedCopyCount, NumberOfReadPairs, HomologySizeIfAvailable(<0ForInsertions), Homology/InsertionSequence
source	chr7:-1-->chr7:54659673-	3.960524250783531	10	None	None
concordant	chr7:54763279+->chr7:54763280-	3.9605242507835308	2	None	None
concordant	chr7:55127267+->chr7:55127268-	3.2711053904267375	0	None	None
concordant	chr7:55155019+->chr7:55155020-	3.271105390426737	1	None	None
source	chr7:-1-->chr7:55669499+	15.951229424626579	5	None	None
concordant	chr7:55669499+->chr7:55669500-	102.34823701028759	86	None	None
source	chr7:-1-->chr7:55760700-	16.303900149839915	9	None	None
concordant	chr7:55760699+->chr7:55760700-	102.34823701028766	81	None	None
concordant	chr7:56049369+->chr7:56049370-	4.313194975996781	3	None	None
source	chr7:-1-->chr7:56114999+	1.1139779336439413	0.0001	None	None
concordant	chr7:56114999+->chr7:56115000-	3.1992170423544333	3	None	None
source	chr7:-1-->chr7:56149666+	3.1992170423544333	3	None	None
discordant	chr7:55155020-->chr7:55127267+	115.02836104448744	68	None	None
discordant	chr7:56049369+->chr7:54763280-	114.33894218413072	48	None	None
```
and [long read graph](https://github.com/AmpliconSuite/CycleExtractor/blob/main/sample_data/GBM39EC_EGFR_graph.txt), with subwalk constraints:
```
SequenceEdge: StartPosition, EndPosition, PredictedCN, AverageCoverage, Size, NumberOfLongReads
sequence	chr7:54659673-	chr7:54763281+	4.227734401796591	45.90736325994846	103609	576
sequence	chr7:54763282-	chr7:55127266+	91.00208132312208	1052.714361855571	363985	40637
sequence	chr7:55127267-	chr7:55155020+	2.8965470408747023	32.72955249693738	27754	172
sequence	chr7:55155021-	chr7:55609190+	91.00208132312208	1013.1828566395843	454170	49675
sequence	chr7:55609191-	chr7:55610094+	2.747919212161151	31.027654867256636	904	915
sequence	chr7:55610095-	chr7:56049369+	91.00208132312208	1023.280632860964	439275	49106
sequence	chr7:56049370-	chr7:56149664+	4.22773440179659	49.62389949648537	100295	562
BreakpointEdge: StartPosition->EndPosition, PredictedCN, NumberOfLongReads
concordant	chr7:54763281+->chr7:54763282-	4.227734401796591	26
concordant	chr7:55127266+->chr7:55127267-	2.8965470408747023	36
concordant	chr7:55155020+->chr7:55155021-	2.8965470408747023	32
concordant	chr7:55609190+->chr7:55609191-	2.747919212161151	38
concordant	chr7:55610094+->chr7:55610095-	2.747919212161151	41
concordant	chr7:56049369+->chr7:56049370-	4.22773440179659	45
discordant	chr7:55155021-->chr7:55127266+	88.10553428224738	978
discordant	chr7:55610095-->chr7:55609190+	88.25416211096093	869
discordant	chr7:56049369+->chr7:54763282-	86.77434692132549	981
PathConstraint: Path, Support
path_constraint	e2+:1,c2-:1,e3+:1,c3-:1,e4+:1	6
path_constraint	e4+:1,c4-:1,e5+:1,c5-:1,e6+:1	34
```
# Output
CE writes the extracted cycles representing potential ecDNA (and s,t-walks, representing linear genomes, if applicable) to the file specified by `--output`. Each cycle includes:
- cycle ID
- its assigned copy number,
- the segments (and orientations) in order composing the cycle, and
- the list of subwalk constraints it can satisfy.

By default, cycles (and walks) are sorted by their copy numbers (see the `--sort-by` option to set up alternative orders by length weighted copy numbers). An example [cycles file](https://github.com/AmpliconSuite/CycleExtractor/blob/main/sample_data/GBM39EC_EGFR_cycles.txt):
```
List of cycle segments
Segment	1	chr7	54659673	54763281
Segment	2	chr7	54763282	55127266
Segment	3	chr7	55127267	55155020
Segment	4	chr7	55155021	55609190
Segment	5	chr7	55609191	55610094
Segment	6	chr7	55610095	56049369
Segment	7	chr7	56049370	56149664
List of longest subpath constraints
Path constraint	1	e2+:1,c2-:1,e3+:1,c3-:1,e4+:1	Support=6	N/A
Path constraint	2	e4+:1,c4-:1,e5+:1,c5-:1,e6+:1	Support=34	N/A
List of extracted cycles/paths
Cycle=1;Copy_count=85.18981770204846;Segments=2+,4+,6+;Path_constraints_satisfied=
```

# Command-Line Options
|Option|Description|Default|
| ------------- | ------------- |------------- |
|`--graph`|	Path to the input breakpoint graph file. **(required)** | -	
|`--output`|	Path to the output file where cycles/paths will be written.	**(required)** | -
|`--gamma`|	Weight balancing number of subwalk constraint satisfied vs. length weighted copy number included ($\sum_{(u,v)\in E}CN(u,v)\cdot l_{uv}$) in MILP objective, larger value favors number of subwalks satisfied.	|`0.01`
|`--sort-by <CopyNumber\|LWCN>`|	Sort output cycles by their copy numbers or length weighted copy number.	|`CopyNumber`
|`--s-t-strategy <all_nodes\|intervals>`| Choose whether to connect source/sink nodes to all nodes or only to interval start/end nodes.	| `all_nodes`
|`--enforce-connectivity`|	CEc mode, which outputs a single connected cycle/s,t-walk at each iteration. CE simply returns the cycle (or walk) with maximum CN when multiple disjoint cycles are produced due to the lack of connectivity constraints, making it potentially faster than CEc. |	off
|`--version`|	Print CE version and exit.|	-

# Starting from bam files

We understand that the user may start from the bam files which is not th input of CE. In that case if the user has long read data, they should follow the steps of CoRAL to get the breakpoint graph https://github.com/AmpliconSuite/CoRAL and if they have short read data, they should follow the Amplicon Architech https://github.com/virajbdeshpande/AmpliconArchitect .

