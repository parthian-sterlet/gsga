# GSGA
Genetic Sensor generation by Genetic Algorithm
GSGA is is abbreviated as Genetic Sensor generation by Genetic Algorithm. GSGA is the software package for the design of a synthetic oligonucleotide as a tandem array of separate transcription factor binding sites (TFBSs or cis-regulatory elements) from DNA sequences 
# Description
The GSGA program complex (PC) is designed to select synthetic polymers (tandem arrays) of from DNA sequence monomer units. These units are separate TFBSs of length from 7-8 to 20-30 base pairs (bp). They are also called cis-regulatory elements or individual binding sites (BS) of the target transcription factor (TF). Target TF implies the TF specific DNA binding to which is expected in the experiment. PC requires a library of DNA motifs [D'haeseleer, 2006](https://doi.org/10.1038/nbt0406-423) representing BS models of different TFs from public databases of TFBS motifs such as [JASPAR](https://jaspar.elixir.no/), [Hocomoco](http://hocomoco13.autosome.ru/) or [Plant Cistrome](http://neomorph.salk.edu/dap_web/pages/index.php). PC applies the stannard motif model of Position Weight Matrix [(Wasserman & Sandelin, 2004)](https://doi.org/10.1038/nrg1315). PC searches for polymers in which the BS quality scores of the target TF are maximized and the BS scores of all non-target TFs are minimized. PC applies two consecutive blocks. The first block creates a polymer, i.e. select participating monomers and defines their order. The second block improves a polymer by introducing certain single-nucleotide substitutions (SNS). PC generates a population of solutions defined for each block. For the first/second blocks these are a set of polymers and a set of SNS patterns, correpondingly. PC evaluates the solutions using a rank-based weighting scheme for motifs of potential TFBSs, and applies genetic algorithms to optimize the population of individuals. PC can be used to design the structure of a polymer and its subsequent incorporation into a transgene promoter for specific expression under the action of a target TF. Default polymer contains ten monomer units, each of them has length about 20-25 bp, so that central 5-10 bp in each unit are essential for TF binding and its 5'/3' flanking 7-10 bp allow more variability.

# Requirements
The source code is written in C++ language. To compile exetubables from the source code you need:

* In Linux system, C++ compiler, e.g. [GCC](https://gcc.gnu.org/) compiler 
* In Windows system any VC++ package, e.g. [Microsoft Visual Studio Community](https://visualstudio.microsoft.com/vs/community/)

# Repository structure
Folder [**src**](https://github.com/parthian-sterlet/gsga/tree/main/src) contains two major C++ source code files.  

Folder [**run**](https://github.com/parthian-sterlet/gsga/tree/main/run) contains three main scripts and the corresponding three main command line examples, implementing different steps of the polymer search.

Folder [**examples**](https://github.com/parthian-sterlet/gsga/tree/main/examples) contains the functional examples of the three main steps of the polymer search.

Folder [**library**](https://github.com/parthian-sterlet/gsga/tree/main/library) contains PWM motifs files and their lists of recognition thresholds 

[Tsukanov et al., 2022](https://doi.org/10.3389/fpls.2022.938545)
[Levitsky et al., 2019](https://doi.org/10.1093/nar/gkz800)
