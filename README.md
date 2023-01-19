<img src="https://github.com/WHops/NAHRwhals/blob/main/NAHRwhals.png?raw=true">

======================================================================================

# NAHRwhals (NAHR-directed Workflow for catcHing seriAL Structural Variations)
An R package providing tools for simulation, detection and visualization of complex, NAHR-driven rearrangements. This includes:
		- Liftover of coordinates between arbitrary human DNA assemblies
		- Accurate sequence alignments of multi-MB DNA sequences
		- Dotplot visualizations 
		- BitLocus - condensed dotplots
		- A tree-based caller for complex, nested NAHR-mediated rearrangements. 

Authors: Wolfram Höps

## Installation

To install NTK from Github, follow the steps below: 

1. Clone the github repository

2. Install required side-packages:
	- minimap2

3. Change the paths in conf/config.yml to fit your system (... step will be abolished in future versions). You can use the 'which' function in the terminal to find the location of your binaries, e.g. 'which minimap2'

4. Edit the path in seqbuilder_functions.R / query_config function to point to your config.yml (... step will be abolished in future versions)

5. run the installation script
    `cd nahrchainer`
    `Rscript install_package.R`
    
6. Verify successful installation
`R` 

` library(nahrtoolkit)`

`confirm_loaded_nahr()`

`>The NAHRtoolkit is loaded and ready to go.`

7. Run an example run

    `cd nahrchainer`
    `Rscript run_a_list.R --example T`

8. Prepare running your own stuff. Whenever we run an assembly, we need to pre-compute an alignment with minimap2. Runtime for two whole human assemblies can be 30m - 2h.

    [optional: open run_a_locus.R and change parameters]

    `minimap2 -t 8 -cx asm20 /your/assembly.fa /your/reference.fa > alignment_sample_hg38.paf`

    `Rscript run_a_locus.R -f /your/reference.fa -g /your/assembly.fa -p alignment_sample_hg38.paf -o testrun_name -i chr1-10000-20000`


## Usage

Please refer to the vignette (/vignettes/) for example use cases of NTK. 
If you are looking to run NTK with many samples, check https://github.com/WHops/ntk-scan-snakemake for a snakemake integration.

## Report Errors

Nahrtoolkit is still in early stages of development, and errors in usage are to be expected. 
Please help improve the code by [reporting](https://github.com/WHops/nahrchainer/issues/new) issues you encounter.

## Correspondence

Please direct any correspondence to: 
Wolfram Höps
wolfram.hoeps@embl.de

