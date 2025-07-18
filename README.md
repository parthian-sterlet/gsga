# GSGA
GSGA is is abbreviated as Genetic Sensor generation by Genetic Algorithm. GSGA is the software package for the design of a synthetic oligonucleotide as a tandem array of separate transcription factor binding sites (TFBSs or cis-regulatory elements) from DNA sequences 

# Description
The GSGA program complex (PC) is designed to select synthetic polymers (tandem arrays) of DNA sequence monomer units. These units are separate TFBSs of lengths from 7-8 to 20-30 base pairs (bp). They are also called cis-regulatory elements or individual binding sites (BS) of the target transcription factor (TF). Target TF implies the TF for which specific DNA binding is expected in the experiment. PC requires the collection of DNA motifs ([D'haeseleer, 2006](https://doi.org/10.1038/nbt0406-423)) representing BS models of different TFs from public databases of TFBS motifs. Currently PC uses BS motifs of *A. thaliana* TFs from [Plant Cistrome](http://neomorph.salk.edu/dap_web/pages/index.php) database, it represents motifs derived by the *in vitro* DAP-seq technology ([O’Malley et al., 2016](https://doi.org/10.1016/j.cell.2016.08.063)). PC can be easily adopted for TFs from other taxa, e.g. mammals (human/mouse) and insects (drosophila), e.g. from [Hocomoco](https://hocomoco13.autosome.org/) ([Vorontsov et al., 2024](https://doi.org/10.1093/nar/gkad1077)) or [JASPAR](https://jaspar.elixir.no/) ([Rauluseviciute et al., 2024](https://doi.org/10.1093/nar/gkad1059)). PC applies the standard motif model of position weight matrix (PWM) [(Wasserman & Sandelin, 2004)](https://doi.org/10.1038/nrg1315). PC searches for polymers in which the recognition scores of the given motif of the target TF BSs is maximized and the recognition scores of all other motifs defining BSs of non-target TFs are minimized. The default polymer contains T<sub>OUT</sub> = 10 monomer units, each unit has length of about 20-25 bp, it is presumed that the central 5-12 bp of each unit are essential for TF binding (core region) and its 5'/3' flanking 7-10 bp allow more variability (non-cores regions). Hence, PC can be used to design the structure of a polymer as an oligonucleotide of about 200-250 bp in length.

For polymer definition PC implements the approach of genetic algorithm (GA). PC applies three consecutive steps:
* the first step (GA<sub>1</sub>) creates the polymer, i.e. it selects participating units and defines their exact order;
* the second step (GA<sub>2</sub>) improves the polymer by introducing single-nucleotide substitutions (SNSs) within non-core regions;
* the third step (GA<sub>3</sub>) destroys  the polymer given as a result of the first or second step by introducing SNSs within the core regions.

Generally GAs were developed according to the principles published earlier ([Levitsky et al., 2007](https://doi.org/10.1186/1471-2105-8-481); [Tsukanov et al. 2022](https://doi.org/10.3389/fpls.2022.938545)). To start each step, each GA generates a population of randomly chosen solutions and gradually improves them to obtain the final set of solutions. For the first and second/third steps the solutions are a set of polymers and a set of SNS patterns for the certain input polymer, correspondingly. PC evaluates the solutions using a rank-based weighting scheme for motifs of potential TFBSs, and applies GAs to optimize the population of individuals, i.e. distinct polymers and SNS patterns for the first and second/third step, respectively. 

The main PC output is an oligonucleotide, it is called genetic sensor or biosensor. The term biosensor here refers to a stably or transiently expressed DNA construct that allow researchers to visualize and quantify certain biologically important ligand (e.g. a phytohormone) biosynthesis, binding, signaling or response; these biosensors consist of a sensory module that responds to the ligand of interest fused to a reporter gene that produces an easy-to-read signal ([Fernandez‐Moreno & Stepanova, 2020](https://doi.org/10.1002/smtd.201900260)). The biosensor incorporation into a transgene promoter for specific expression under the action of a target TF can be used as a marker of presence of a target TF mediating the ligand action. 

# Requirements
The source code is written in C++ language. To compile exetubables from the source code it is required:

* In Linux system, C++ compiler, e.g. [GCC](https://gcc.gnu.org/) compiler 
* In Windows system any VC++ package, e.g. [Microsoft Visual Studio Community](https://visualstudio.microsoft.com/vs/community/)

# Repository structure
Folder [**src**](https://github.com/parthian-sterlet/gsga/tree/main/src) contains two major C++ source code files, and two supplementary C++ source code files are provided in [MCOT repository](https://github.com/parthian-sterlet/mcot-kernel).  

Folder [**run**](https://github.com/parthian-sterlet/gsga/tree/main/run) contains the script [script](https://github.com/parthian-sterlet/gsga/blob/main/run/build.sh) for source code compliation, and three scripts implementing three command line examples, corresponding three steps.

Folder [**examples**](https://github.com/parthian-sterlet/gsga/tree/main/examples) contains input files as examples of the three main steps.

Folder [**matrices**](https://github.com/parthian-sterlet/gsga/tree/main/matrices) represents the collection of PWMs, it contains two types of files: 
* files of PWMs as matrices of position-specific log-odds weights computed as desribed previously (Levitsky et al., [2007](https://doi.org/10.1186/1471-2105-8-481), [2019](https://doi.org/10.1093/nar/gkz800)), and
* files of lists of recognition thresholds computed with the whole-genome set of promoters of protein coding genes, see Levitsky et al. ([2019](https://doi.org/10.1093/nar/gkz800)) and Tsukanov et al. ([2022](https://doi.org/10.3389/fpls.2022.938545)).

To make the recognition for all DNA motifs uniform, each threshold respects the expected recognition rate (ERR) value, this value means the probability of site recognition per one tested bp position in all promoters ([Levitsky et al., 2024](https://doi.org/10.18699/vjgb-24-90)); hence, each file of the list contains two columns: (a) recognition threshold and (b) ERR value as Log<sub>10</sub>(probability), where the probability is ratio of the number of positions with recognized BSs to the total number of tested positions in promoters.

# Pipeline description

# Common input data for all three steps
* a set of weight matrices for non-target TFs;
* a set of lists of recognition thresholds and respective ERR values for weight matrices ([Tsukanov et al., 2022](https://doi.org/10.3389/fpls.2022.938545); Levitsky et al., [2019](https://doi.org/10.1093/nar/gkz800), [2024](https://doi.org/10.18699/vjgb-24-90)).

## First step, GA<sub>1</sub>
To generate a new polymer, the default number of T<sub>OUT</sub> = 10 monomers are required in a polymer. The total number of distinct input monomer units T<sub>IN</sub> should be higher, T<sub>IN</sub> >= T<sub>OUT</sub> = 10, to support the polymer specificity. These T<sub>IN</sub> monomers are presumed to be the native DNA sequences supported by ChIP-seq/RNA-seq etc. experimental edidence of specific binding of the target TF. The task of the first step is dual: 
1. to select exact output T<sub>OUT</sub> monomers among the total T<sub>IN</sub> provided in input data; 
2. to denote the exact order of these T<sub>OUT</sub> selected monomers.

For example, if we have 20 input monomers {M<sub>1</sub>, M<sub>2</sub>, ... M<sub>20</sub>}, then the example version of the ouput order of the selected top-scored ten monomers is {M<sub>17</sub>, M<sub>2</sub>, M<sub>5</sub>, M<sub>13</sub>, M<sub>4</sub>, M<sub>1</sub>, M<sub>18</sub>, M<sub>9</sub>, M<sub>15</sub>, M<sub>11</sub>}. To find an optimal combination, GA of the first step (GA<sub>1</sub>) selects the multiple versions of polymers with the least susceptibility to the non-target TFs binding. As input data, GA<sub>1</sub> requires:
* the set of T<sub>IN</sub> monomer units containing individual binding sites of the target TF; 
* the number of T<sub>OUT</sub> monomer units of a polymer; 
* the weight matrix for the target TF, and a list of its recognition thresholds and respective ERRs for this matrix.

## Second step, GA<sub>2</sub>
The second step (GA<sub>2</sub>) selects appropriate SNS outside the essential positions the target TF binding in each monomer. Hence, GA<sub>2</sub> requires: 
* the polymer assembled from the units comprising the target TF binding site, each site contains the essential core flanked by several nucleotides on 5' and 3' sides (non-cores); 
* the weight matrix for the target TF, and the list of its recognition thresholds and respective ERR values; 
* the list of positions in the polymer designating non-core elements, and flanking regions before/after the first/last monomer units of the polymer; 
* the probability of SNSs within non-core elements. 

## Third step, GA<sub>3</sub>
The third step (GA<sub>3</sub>) is another application of approach developped for the preceeding second step (GA<sub>2</sub>). Here the same source code is applied to destroy any DNA binding motif, hence it is not important here BSs of which TF to exclude. Hence, BSs of neither target nor non-target TFs are now undesirable. Hence, GA<sub>3</sub> requires:
* the polymer assembled from the units comprising the target TF binding site, each site contains the essential core flanked by several nucleotides on 5' and 3' sides (non-cores), this polymer may be the result of either the first or second step;
* the list of positions in the polymer designating the cores regions;
* the probability of SNSs within core elements.

# How to compile

## In Linux system: 

```
git clone https://github.com/parthian-sterlet/gsga

cd gsga/run

chmod a+x build.sh

./build.sh
```

## In Windiws system:

separate compilation of all source files in VC++

# Command line arguments

## 1. First step: Selection of a subset of total native units and definition of their order in a polymer
[Order](https://github.com/parthian-sterlet/gsga/blob/main/src/genosensor_seq_order_ga.cpp) propgram defines the composition of units and their order. 

1. path to files of (a) weight matrix for the target TF DNA motif, [see example matrix](https://github.com/parthian-sterlet/gsga/blob/main/examples/order/dapseq0.pwm) and (b) its threshold list, this file contains the list of pairs {Threshold, -Log<sub>10</sub>(ERR)} values, see [example distribution](https://github.com/parthian-sterlet/gsga/blob/main/examples/order/dapseq0.dist). The last symbol of path must be '/' and '\' for Linux and Windows OS, respectively.
2. path to files of (a) all non-target TFs DNA motifs and (b) their threshold lists, these files contain the lists of pairs {Threshold, -Log<sub>10</sub>(ERR)} values. The last symbol of path must be '/' and '\' for Linux and Windows OS, respectively.
3. input file in FASTA format representing the alignment of proven target TF BSs, these are total amount of monomer units that are used to generate momomer.
4. integer value, count of selected monomer units in a polymer, T<sub>OUT</sub>, default number is 10.
5. integer value, count of motifs in the collection, default number is 528, it implies DNA motifs of *A.thaliana* TFs from [Plant Cistrome](http://neomorph.salk.edu/dap_web/pages/index.php) database, from DAP-seq experoment ([O’Malley et al., 2016](https://doi.org/10.1016/j.cell.2016.08.063))
6. char* value, name of motif file, the default value "dapseq" means (a) for the non-target TFs: the weight matrix files are dapseq1.pwm, dapseq2.pwm, etc. up to dapseq528.pwm, and threshold list files are dapseq1.dist, dapseq2.dist, etc. up to dapseq528.dist, (b) for the target TF the weight matrix file is dapseq0.pwm and the threshold list file is dapseq0.dist.
7. output file listing results, i.e. the multiple solutions in the FASTA format in the descending order of the qulity.
8. output log file showing the progress in calculation.

## 2. Second & Third steps: Improve a polymer of native units by SNSs within the non-core regions of polymer & Destroy a polymer of native or synthetic units by SNSs withiin the core regions
[Improve and Destroy](https://github.com/parthian-sterlet/gsga/blob/master/src/genosensor_seq_ga.cpp) program implements two distict tasks by introducing mutations of nucleotides: 
* Second step: SNSs are allowed only in the non-core regions of a polymer thereby improving the binding of the target TF and restricting the binding of non-target TFs, or
* Third step: SNSs are allowed only in the non-core regions of a polymer thereby restricting the binding of any TFs either the target or non-target. 

1. path to files of (a) weight matrix for the target TF DNA motif, see [example matrix](https://github.com/parthian-sterlet/gsga/blob/main/examples/improve/dapseq0.pwm) and (b) its threshold list, this file contains the list of pairs {Threshold, -Log<sub>10</sub>(ERR)} values, [example distribution](https://github.com/parthian-sterlet/gsga/blob/main/examples/improve/dapseq0.dist). The last symbol of path must be '/' and '\' for Linux and Windows OS, respectively.
2. path to files of (a) weight matrices for all non-target TFs and (b) their threshold lists, these files contain the lists of pairs {Threshold, -Log<sub>10</sub>(ERR)} values. The last symbol of path must be '/' and '\' for Linux and Windows OS, respectively.
3. input file in FASTA format with a DNA sequence of polymer selected by the previous analysis step, the first step, [Order](https://github.com/parthian-sterlet/gsga/blob/main/src/genosensor_seq_order_ga.cpp)
4. file of tab-delimited table, this table marks positions in the polymer non-cores regions (spacers between the essential cores) and flanking sequences before/after the first/last monomers of the polymer
5. integer value, count of motifs in the collection, default number is 529, it implies 528 DNA motifs of *A.thaliana* TFs from [Plant Cistrome](http://neomorph.salk.edu/dap_web/pages/index.php) database, from DAP-seq experoment ([O’Malley et al., 2016](https://doi.org/10.1016/j.cell.2016.08.063))
6. char* value, name of motif file, the default value "dapseq" means (a) for the non-target TFs: the weight matrix files are dapseq1.pwm, dapseq2.pwm, etc. up to dapseq528.pwm, and threshold list files are dapseq1.dist, dapseq2.dist, etc. up to dapseq528.dist, (b) for the target TF the weight matrix file is dapseq0.pwm and the threshold list file is dapseq0.dist.
7. output file, log file listing results, i.e. the multiple solutions in the descending order of the qulity.
8. integer value, the anchor mode. The values 1 or 0 mean the Improve/Destroy option respecting the Second/Third steps of analysis. In these cases a specific weight matrix of a target TF is opposed / is not opposed to matrices of all non-target TFs.
9. double value, the probability P of SNSs within designated elements, P value is equal to the ratio between the number of SNSs and the total sequence length allowed for these SNSs, the number of SNSs is the same for each non-core/core region between two neighbor core/non-core regions for Improve/Desrtoy steps.
10. output log file showing the progress in calculation.

## Supplementary programs
These programs are described in [MCOT repository](https://github.com/parthian-sterlet/mcot-kernel/tree/master#generation-of-partner-library)
* The first program [pfm_to_pwm_mat.cpp](https://github.com/parthian-sterlet/mcot-kernel/blob/master/src/pfm_to_pwm/pfm_to_pwm_mat.cpp) converts a DNA motif (position frequency matrix, PFM) to PWM (weights used to recognize BSs). This program computes the weight matrix for given DNA motif, see  [example DNA motif](https://github.com/parthian-sterlet/gsga/blob/main/examples/order/dapseq0.motif) and [example weight matrix](https://github.com/parthian-sterlet/gsga/blob/main/examples/order/dapseq0.pwm)
* The second program [pwm_iz_pwm_thr_dist0.cpp](https://github.com/parthian-sterlet/mcot-kernel/blob/master/src/pwm_thr_err/pwm_iz_pwm_thr_dist0.cpp) computes the list consisting of pairs of values {Threshold, -Log<sub>10</sub>(ERR)} for the DNA motif defined by PFM and PWM. This program is used for the uniform recognition of DNA motifs, here the term uniform implies that recognition threshods are normilized according the ERR values, see [example distribution](https://github.com/parthian-sterlet/gsga/blob/main/examples/improve/dapseq0.dist)) for a given [example weight matrix](https://github.com/parthian-sterlet/gsga/blob/main/examples/order/dapseq0.pwm).

# Correction of the list of DNA motifs of non-target TFs according to given input DNA motif of the target TF
The presence among motifs of non-target factors of motifs very similar to the selected motive of the target factor very negatively affects the quality of results of programs of all three steps. The examples provided in this repository show the G-rich non-canonical motif of the EIN3 TF as the motif of the target TFs (see file [dapseq0.motif](https://github.com/parthian-sterlet/gsga/blob/main/examples/order/dapseq0.motif) respecting to the PWM of the example target TF [dapseq0.pwm](https://github.com/parthian-sterlet/gsga/blob/main/examples/order/dapseq0.pwm). This motif does not show the high similrity to any motif from the DAP-seq collection of motifs accepted here in analysis. To perform the correct selection of any arbitrary target TF and the provided here list of non-target TFs from DAP-seq user should test the similarity of tested motif of the target TF. The significane of similarity of the motif of the target TF should be estimated by the  stadard motif comparison tool [TomTom](https://meme-suite.org/meme/tools/tomtom) (the option Select a motif database or provide motifs to compare with = 'ARABIDOPSIS (Arabidopsis thaliana) DNA. DAP motifs (O'Malley2016)'). If user found highly simialar DAP motifs, the list of motifs of presumed non-target TFs should be corrected. [XLSX file](https://github.com/parthian-sterlet/gsga/blob/main/matrices/DAP-seq_Plant_Cistrome_528_motifs_for_GSGA.xlsx) shows the list of all motifs used in both c++ files, [Order](https://github.com/parthian-sterlet/gsga/blob/main/src/genosensor_seq_order_ga.cpp) and [Improve and Destroy](https://github.com/parthian-sterlet/gsga/blob/master/src/genosensor_seq_ga.cpp). In any of these c++ file sees the piece just after the declaration of the main function

``` int main(int argc, char *argv[]) ``` 

Next, after parsing the command line arguments there is a line to select ignored motifs of non-target TFs, this line starts with 

```int m_ignore[] = ```

but not 

``` // int m_ignore[] = ```

since the first symbols in a line '//' mark a [comment in c++ language](https://learn.microsoft.com/en-us/cpp/cpp/comments-cpp?view=msvc-170).

The default content of this line in both files [Order](https://github.com/parthian-sterlet/gsga/blob/main/src/genosensor_seq_order_ga.cpp) and [Improve and Destroy](https://github.com/parthian-sterlet/gsga/blob/master/src/genosensor_seq_ga.cpp) is 

```int m_ignore[] = { 2, 102, 183, 184, 186, 207, 212, 213, 217, 227, 286, 299, 300, 303, 439, -1 };// no anchor 2024 ```

The content of this line respects the default option in ignoring motifs, the program ignores only too short and degenerate DAP motifs, this is explained in the MCOT paper ([Levitsky et al., 2019](https://doi.org/10.1093/nar/gkz800)). The numbers notation of all motifs showing the high sinilarity to the tested motif of the target TF should be added to the list of this line. Note that the standard statistical criterion p-value < 0.05 is too mild, the Bonferroni correction is p-adjusted < 0.05 / 528 < 0.0001, at least this threshold may be applied to find the significantly similar motifs. 

# Examples command lines:

These example command lines implement various steps for Linux OS:
1. [Order](https://github.com/parthian-sterlet/gsga/blob/master/src/order) - First step defines the composition and monomer order of native units
2. [Improve](https://github.com/parthian-sterlet/gsga/blob/master/src/imrove) - Second step improves a polymer of native units to a polymer of synthetic units by SNSs within the non-core regions
3. [Destroy](https://github.com/parthian-sterlet/gsga/blob/master/src/destroy) - Third step destroys a polymer of native or synthetic units by SNSs withiin the core regions
