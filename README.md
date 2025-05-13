# gsga
Genetic Sensor generation by Genetic Algorithm
Genetic Sensor generation by Genetic Algorithm or GSGA is the software package for the design of a synthetic oligonucleotide as a tandem array of separate transcription factor binding sites (TFBSs or cis-regulatory elements) from DNA sequences 
Description
The GSGA program complex (PC) is designed for the design of synthetic polymers (tandem arrays) of from DNA sequence monomer units. These units are separate TFBSs, or cis-regulatory elements as individual binding sites (BS) of the target transcription factor (TF). 
The PC requires a library of BS models of different TFs from public databases such as JASPAR, Hocomoco or Plant Cistrome. The PC searches for polymers in which the BS quality scores of the target TF are maximized and the BS scores of all non-target TFs are minimized. 
The PC applies two consecutive blocks. The first block creates a polymer, i.e. select participating monomers and defines their order. The second block improves a polymer by introducing certain single-nucleotide substitutions (SNS). 
The PC generates a population of solutions defined for each block. For the first/second blocks these are a set of polymers and a set of SNS patterns, correpondingly. 
The PC evaluates the solutions using a rank-based weighting scheme for motifs of potential TFBSs, and applies genetic algorithms to optimize the population of individuals. 
PC can be used to design the structure of a polymer and its subsequent incorporation into a transgene promoter for specific expression under the action of a target TF.
