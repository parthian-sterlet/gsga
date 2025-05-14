# GSGA
Genetic Sensor generation by Genetic Algorithm
GSGA is is abbreviated as Genetic Sensor generation by Genetic Algorithm. GSGA is the software package for the design of a synthetic oligonucleotide as a tandem array of separate transcription factor binding sites (TFBSs or cis-regulatory elements) from DNA sequences 
# Common description
The GSGA program complex (PC) is designed to select synthetic polymers (tandem arrays) of from DNA sequence monomer units. These units are separate TFBSs of length from 7-8 to 20-30 base pairs (bp). They are also called cis-regulatory elements or individual binding sites (BS) of the target transcription factor (TF). Target TF implies the TF specific DNA binding to which is expected in the experiment. PC requires a library of DNA motifs ([D'haeseleer, 2006](https://doi.org/10.1038/nbt0406-423)) representing BS models of different TFs from public databases of TFBS motifs. Currently PC is adopted only for motifs of *A.thaliana* TFs from [Plant Cistrome](http://neomorph.salk.edu/dap_web/pages/index.php) database, it represents motifs derived by the *in vitro* DAP-seq technology ([O’Malley et al., 2016](https://doi.org/10.1016/j.cell.2016.08.063)). PC can be easily adopted for TFs from other taxa, e.g. mammals (human/mouse) and insects (drosophila), e.g. from [Hocomoco](https://hocomoco13.autosome.org/) ([Vorontsov et al., 2024](https://doi.org/10.1093/nar/gkad1077)) or [JASPAR](https://jaspar.elixir.no/) ([Rauluseviciute et al., 2024](https://doi.org/10.1093/nar/gkad1059)). PC applies the stannard motif model of Position Weight Matrix [(Wasserman & Sandelin, 2004)](https://doi.org/10.1038/nrg1315). PC searches for polymers in which the BS quality scores of the target TF are maximized and the BS scores of all non-target TFs are minimized. PC applies two consecutive blocks. The first block creates a polymer, i.e. select participating monomers and defines their order. The second block improves a polymer by introducing certain single-nucleotide substitutions (SNS). PC generates a population of solutions defined for each block. For the first/second blocks these are a set of polymers and a set of SNS patterns, correpondingly. PC evaluates the solutions using a rank-based weighting scheme for motifs of potential TFBSs, and applies genetic algorithms to optimize the population of individuals. PC can be used to design the structure of a polymer and its subsequent incorporation into a transgene promoter for specific expression under the action of a target TF. Default polymer contains TOut =10 monomer units, each of them has length about 20-25 bp, so that central 5-12 bp in each unit are essential for TF binding (core) and its 5'/3' flanking 7-10 bp allow more variability (non-cores).

# Requirements
The source code is written in C++ language. To compile exetubables from the source code you need:

* In Linux system, C++ compiler, e.g. [GCC](https://gcc.gnu.org/) compiler 
* In Windows system any VC++ package, e.g. [Microsoft Visual Studio Community](https://visualstudio.microsoft.com/vs/community/)

# Repository structure
Folder [**src**](https://github.com/parthian-sterlet/gsga/tree/main/src) contains two major C++ source code files.  

Folder [**run**](https://github.com/parthian-sterlet/gsga/tree/main/run) contains the main [script](https://github.com/parthian-sterlet/gsga/blob/main/run/build.sh) for source code compliation, the three main scripts and the corresponding three command line examples, implementing different steps of the polymer search.

Folder [**examples**](https://github.com/parthian-sterlet/gsga/tree/main/examples) contains the functional examples of the three main steps of the polymer search.

Folder [**library**](https://github.com/parthian-sterlet/gsga/tree/main/library) represents the library of DNA motifs, it contains two types of files: (1) files of PWM motifs as matrices of log-odds weights computed as desribed previously ([Levitsky et al., 2007](https://doi.org/10.1186/1471-2105-8-481); [Levitsky et al., 2019](https://doi.org/10.1093/nar/gkz800)) and (2) files of lists of recognition thresholds computed with the whole-genome set of promoters of protein coding genes, see Levitsky et al. ([2019](https://doi.org/10.1093/nar/gkz800)) and Tsukanov et al. ([2022](https://doi.org/10.3389/fpls.2022.938545)); to make the recognition for all DNA motifs uniform, each threshold respects the expected recognition rate (ERR) value, it means the probability recognition per one bp position in the promoters ([Levitsky et al., 2024](https://doi.org/10.18699/vjgb-24-90)); hence, each file of the list contains two columns: (a) recognition threshold and (b) ERR value as Log10(probability), where the probability is ratio of the positions with recognized BSs (in either DNA strand) to the total number of tested positions in promoters.

# Pipeline description

# Common input data for all blocks
* a set of weight matrices for off-target TFs;
* a set of preliminarily computed lists of thresholds and respective ERRs for each matrix ([Tsukanov et al. 2022](https://doi.org/10.3389/fpls.2022.938545); [Levitsky et al., 2019](https://doi.org/10.1093/nar/gkz800); [Levitsky et al., 2024](https://doi.org/10.18699/vjgb-24-90)).

## First block, GA1, Ante mare undae
To generate a new polymer, the default number of TOut = 10 monomers are stacked in a polymer. Note the defthat number of distinct input monomers TIn should be higher, TOut >= TIn = 10, to support the polymer specificity. Note that these TIn monomers are presumed to be the native DNA sequences supporte by ChIP-seq/RNA-seq etc. experimental edidence of specific binding of the target TF. The task of the first block is dual: 
1. to select exact output Tout monomers among the total TIn provided in input data; 
2. to denote the exact order of TOut selected monomers.

For example, let we have 20 input monomers {T1, T2, ... T20}, then the version of the ouput order is {T17, T2, T5, T13, T4, T1, T18, T9, T15, T11}. To find an optimal combination, a genetic algorithm of the first block (GA1) selects the multiple versions of polymers with the least susceptibility to the non-target TFs binding. As input data, GA1 requires:
* a set of TIn sequences (monomers) containing individual binding sites of the target TF; 
* the number of TOut monomer units of a polymer; 
* a matrix for the target TF, and a list of its recognition thresholds and respective ERRs for this matrix.

## Second block, GA2, Non est terminus ad perfectionem
The second block the second genetic algorithm (GA2) selects appropriate SNS outside the essential positions the target TF binding in each monomer. Hence, GA2 requires: 
* a polymer assembled from the units comprising the target TF binding site (the essential core) flanked by several nucleotides on 5' and 3' sides (less essential flanks, non-cores); 
* a matrix for the target TF, and a list of its recognition thresholds and respective ERRs for this matrix; 
* a list of positions in the polymer (the assembled sequence from the first block) designating spacers between the essential cores and non-core elements, and flanking regions before/after the first/last monomer units of the polymer; 
* the probability P of nucleotide substitutions (SNS) within non-core elements. 

## Third block, GA3, Caedite eos. Novit enim Dominus qui sunt eius
The third block is another application of approach developped for the preceedin second block. Here the same source code is applied to destroy any DNA binding motif, hence it is not important here BSs of which TF to exclude. Hence, BSs of neither target nor non-target TFs are now undesirable. Hence, GA3 requires:
* a polymer assembled from the units comprising the target TF binding site (the essential core) flanked by several nucleotides on 5' and 3' sides (less essential flanks, non-cores), this polymer may be the result of either the first or second block;
* a list of positions in the polymer designating the essential cores between the non-core regions and flanking regions before/after the first/last essential cores of the polymer;
* the probability P of nucleotide substitutions (SNS) within core elements.

# How to compile
* In Linux system: 

git clone https://github.com/parthian-sterlet/gsga

cd gsga/run

chmod a+x build.sh

./build.sh

* In Windiws system:

separate compilation of all source files in VC++

# Command line arguments

## 1. Selection of a subset of total native units and definition of their order in a polymer
[Order](https://github.com/parthian-sterlet/gsga/blob/main/src/genosensor_seq_order_ga.cpp) propgram defines the composition of units and their order. 
1. path to files of (a) the target TF DNA motif and (b) its threshold list, this file contains the list of pairs {Threshold, -10(ERR)} values. The last symbol of path must be '/' and '\' for Linux and Windows OS, respectively.
2. path to files of (a) all non-target TFs DNA motifs and (b) their threshold lists, these files contain the lists of pairs {Threshold, -10(ERR)} values. The last symbol of path must be '/' and '\' for Linux and Windows OS, respectively.
3. input file in FASTA format of total amount of monomer units that can be used to generate momomer.
4. integer value, count of selected monomer units in a polymer, default number is 10.
5. integer value, count of motifs in library, default number is 528, it implies DNA motifs of *A.thaliana* TFs from [Plant Cistrome](http://neomorph.salk.edu/dap_web/pages/index.php) database, from DAP-seq experoment ([O’Malley et al., 2016](https://doi.org/10.1016/j.cell.2016.08.063)
6. char name of motif file, the default value "dapseq" means (a) for non-target TFs: the motif files dapseq1.pwm, dapseq2.pwm, etc. up to dapseq528.pwm, and threshold list files dapseq1.dist, dapseq2.dist, etc. up to dapseq528.dist, (b) for the target TF the motif file. .dapseq0.pwm and the threshold list file dapseq0.dist.
7. output file listing results, i.e. the multiple solutions in the FASTA format in the descending order of the qulity.
8. output log file showing the progress in calculation.

## 2. Improvement of a polymer of native units by single nucleotide mutations within non-core regions of polymer
[Improvement](https://github.com/parthian-sterlet/antinoise/blob/master/src/mix0.cpp) program introduce mutation of nucleotides in non-core regions of a polymer.
1. path to files of (a) the target TF DNA motif and (b) its threshold list, this file contains the list of pairs {Threshold, -Log10(ERR)} values. The last symbol of path must be '/' and '\' for Linux and Windows OS, respectively.
2. path to files of (a) all non-target TFs DNA motifs and (b) their threshold lists, these files contain the lists of pairs {Threshold, -Log10(ERR)} values. The last symbol of path must be '/' and '\' for Linux and Windows OS, respectively.
3. input file in FASTA format with a DNA sequence of polymer selected by the previous analysis step, the first block, [Order](https://github.com/parthian-sterlet/gsga/blob/main/src/genosensor_seq_order_ga.cpp)
4. file of tab-delimited table, this table marks positions in the polymer non-cores regions (spacers between the essential cores) and flanking sequences before/after the first/last monomers of the polymer
5. integer value, count of motifs in library, default number is 529, it implies 528 DNA motifs of *A.thaliana* TFs from [Plant Cistrome](http://neomorph.salk.edu/dap_web/pages/index.php) database, from DAP-seq experoment ([O’Malley et al., 2016](https://doi.org/10.1016/j.cell.2016.08.063))
6. char name of motif file, the default value "dapseq" means (a) for non-target TFs: the motif files dapseq1.pwm, dapseq2.pwm, etc. up to dapseq528.pwm, and threshold list files dapseq1.dist, dapseq2.dist, etc. up to dapseq528.dist, (b) for the target TF the motif file dapseq0.pwm and the threshold list file dapseq0.dist.
7. output file, log file listing results, i.e. the solution in the descending order of the qulity.
8. integer value, the anchor mode. The values 1 or 0 mean that a specific matrix of a target TF is used / not used for optimization.
9. double value, the probability P of nucleotide substitutions (SNS) within designated elements, P value is equal to the ratio between the number of mutation and sequence length, the number of substitutions is the same for each non-core spacer between two neighbor core regions.
10. output log file showing the progress in calculation.

# Examples scripts:

These scripts implement various blocks for Linux OS:
1. [Order](https://github.com/parthian-sterlet/gsga/blob/master/src/order) - First block defines the composition and monomer order of native units
2. [Improvement](https://github.com/parthian-sterlet/gsga/blob/master/src/imrovemnet) - Second block improves a polymer of native units to a polymer of synthetic units by single nucleotide mutations within the non-core regions
3. [Degradation](https://github.com/parthian-sterlet/gsga/blob/master/src/degradation) - Third block destroys a polymer of synthetic units by single nucleotide mutations withiin the core regions
