# GSGA
Genetic Sensor generation by Genetic Algorithm
GSGA is is abbreviated as Genetic Sensor generation by Genetic Algorithm. GSGA is the software package for the design of a synthetic oligonucleotide as a tandem array of separate transcription factor binding sites (TFBSs or cis-regulatory elements) from DNA sequences 
# Common description
The GSGA program complex (PC) is designed to select synthetic polymers (tandem arrays) of from DNA sequence monomer units. These units are separate TFBSs of length from 7-8 to 20-30 base pairs (bp). They are also called cis-regulatory elements or individual binding sites (BS) of the target transcription factor (TF). Target TF implies the TF specific DNA binding to which is expected in the experiment. PC requires a library of DNA motifs [D'haeseleer, 2006](https://doi.org/10.1038/nbt0406-423) representing BS models of different TFs from public databases of TFBS motifs. Currently PC is adopted only for motifs of *A.thaliana* TFs from [Plant Cistrome](http://neomorph.salk.edu/dap_web/pages/index.php) database, it represents motifs derived by the *in vitro* DAP-seq technology ([O’Malley et al., 2016](https://doi.org/10.1016/j.cell.2016.08.063). PC can be easily adopted for TFs from other taxa, e.g. mammals (human/mouse) and insects (drosophila), e.g. from [Hocomoco](https://hocomoco13.autosome.org/) [Vorontsov et al., 2024](https://doi.org/10.1093/nar/gkad1077) or [JASPAR](https://jaspar.elixir.no/) [Rauluseviciute et al., 2024](https://doi.org/10.1093/nar/gkad1059). PC applies the stannard motif model of Position Weight Matrix [(Wasserman & Sandelin, 2004)](https://doi.org/10.1038/nrg1315). PC searches for polymers in which the BS quality scores of the target TF are maximized and the BS scores of all non-target TFs are minimized. PC applies two consecutive blocks. The first block creates a polymer, i.e. select participating monomers and defines their order. The second block improves a polymer by introducing certain single-nucleotide substitutions (SNS). PC generates a population of solutions defined for each block. For the first/second blocks these are a set of polymers and a set of SNS patterns, correpondingly. PC evaluates the solutions using a rank-based weighting scheme for motifs of potential TFBSs, and applies genetic algorithms to optimize the population of individuals. PC can be used to design the structure of a polymer and its subsequent incorporation into a transgene promoter for specific expression under the action of a target TF. Default polymer contains TOut =10 monomer units, each of them has length about 20-25 bp, so that central 5-10 bp in each unit are essential for TF binding (core) and its 5'/3' flanking 7-10 bp allow more variability (non-cores).

# Requirements
The source code is written in C++ language. To compile exetubables from the source code you need:

* In Linux system, C++ compiler, e.g. [GCC](https://gcc.gnu.org/) compiler 
* In Windows system any VC++ package, e.g. [Microsoft Visual Studio Community](https://visualstudio.microsoft.com/vs/community/)

# Repository structure
Folder [**src**](https://github.com/parthian-sterlet/gsga/tree/main/src) contains two major C++ source code files.  

Folder [**run**](https://github.com/parthian-sterlet/gsga/tree/main/run) contains three main scripts and the corresponding three main command line examples, implementing different steps of the polymer search.

Folder [**examples**](https://github.com/parthian-sterlet/gsga/tree/main/examples) contains the functional examples of the three main steps of the polymer search.

Folder [**library**](https://github.com/parthian-sterlet/gsga/tree/main/library) represents the library of DNA motifs, it contains two types of files: (1) files of PWM motifs as matrices of log-odds weights computed as desribed previously ([Levitsky et al., 2007](https://doi.org/10.1186/1471-2105-8-481); [Levitsky et al., 2019](https://doi.org/10.1093/nar/gkz800)) and (2) files of lists of recognition thresholds computed with the whole-genome set of promoters of protein coding genes, see Levitsky et al. ([2019](https://doi.org/10.1093/nar/gkz800)) and Tsukanov et al. ([2022](https://doi.org/10.3389/fpls.2022.938545)); to make the recognition for all DNA motifs uniform, each threshold respects the expected recognition rate (ERR) value, it means the probability recognition per one bp position in the promoters ([Levitsky et al., 2024](https://doi.org/10.18699/vjgb-24-90)); hence, each file of the list contains two columns: (a) recognition threshold and (b) ERR value as LOG10(probability), where the probability is ratio of the positions with recognized BSs (in either DNA strand) to the total number of tested positions in promoters.

# Input data description
). As input data, GA1 uses (1) a set of N sequences (monomers) containing the target TF binding site, and (2) a set of matrices for off-target TFs with (3) a preliminarily computed list of thresholds and respective FPRs for each matrix as described in (Levitsky et al., 2019). 

# Pipeline description

# Common input data for all blocks
(1) a set of weight matrices for off-target TFs; (2) a set of preliminarily computed lists of thresholds and respective ERRs for each matrix ([Tsukanov et al. 2022](https://doi.org/10.3389/fpls.2022.938545); [Levitsky et al., 2019](https://doi.org/10.1093/nar/gkz800); [Levitsky et al., 2024](https://doi.org/10.18699/vjgb-24-90)).

## First block
To generate a new polymer, the default number of TOut = 10 monomers are stacked in a polymer. Note the defthat number of distinct input monomers TIn should be higher, TOut >= TIn = 10, to support the polymer specificity. Note that these TIn monomers are presumed to be the native DNA sequences supporte by ChIP-seq/RNA-seq etc. experimental edidence of specific binding of the target TF. The task of the first block is dual: (1) to select exact output Tout monomers among the total TIn provided in input data; (2) to denote the exact order of TOut selected monomers. For example, let we have 20 input monomers {T1, T2, ... T20}, then the version of the ouput order is {T17, T2, T5, T13, T4, T1, T18, T9, T15, T11}. To find an optimal combination, a genetic algorithm of the first block (GA1) selects the multiple versions of polymers with the least susceptibility to off-target TFs binding. 

## Second block
The second block the second genetic algorithm (GA2) selects appropriate SNS outside the essential positions the target TF binding in each monomer. Hence, GA2 requires (1) a polymer assembled from the units comprising the target TF binding site (the essential core) flanked by several nucleotides (less essential flanks) on 5' and 3' sides; (2) a matrix for the target TF; (3) a list of thresholds and respective ERRs for this matrix; (4) a list of positions in the assembled sequence designating spacers between the essential cores  and non-core elements flanking the assembled sequence; (5) the probability p of nucleotide substitutions  (SNS) within designated elements. 

# Command line arguments

## 1. Genomic background sequence generation approach
The [major](https://github.com/parthian-sterlet/antinoise/blob/main/src/background_genome_mono.cpp) propgram of this tool finds the genomic background sequences for a particular genome (hg38, mm10, tair10, etc.). The background sequences match almost perfectly A/T content. The background sequences either match exactly the length of DNA sequences from the foreground set, or user defines the same length for all background sequences using a special parameters.
## 2. Synthetic background sequence generation approach
The alternative program [mix0.cpp](https://github.com/parthian-sterlet/antinoise/blob/master/src/mix0.cpp) generates synthetic background sequences that exactly match the nucleotide content of the foreground sequences.
