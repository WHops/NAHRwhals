<img src="https://github.com/WHops/NAHRwhals/blob/main/NAHRwhals.png?raw=true">


# NAHRwhals (NAHR-directed Workflow for catcHing seriAL Structural Variations)
An R package providing tools for visualization and automatic detection of complex, NAHR-driven rearrangements (few kbp to multiple Mbp) using genome assemblies. Modules include:

- Liftover of coordinates between arbitrary human DNA assemblies
- Accurate sequence alignments of multi-MB DNA sequences
- Dotplot visualizations 
- Segmented Dotplots
- A tree-based caller for complex, nested NAHR-mediated rearrangements. 


# Dependencies

NAHRwhals requires installed versions of: 
- minimap2 (link)
- bedtools (link)
- gawk (link)
- R >= 4.1.0
#  Installation

### (1) Clone the repository

`git clone --recursive https://github.com/WHops/NAHRwhals.git`

`cd NAHRwhals`

### (2) Install with install_package.R, which calls devtools_install() 

`Rscript install_package.R`

### (3) [OPTIONAL] specify minimap2 / bedtools paths
In case minimap2 or bedtools are not part of your $PATH (i.e. can not be called from the commandline via `minimap2` and `bedtools`), specify your paths in conf/config.txt: (otherwise, leave them as 'default')

```
minimap2_bin = '/your/path/to/minimap2'
bedtools_bin = '/your/path/to/bedtools'
```


#  Test & Example runs

To confirm that NAHRwhals has been correctly installed, run a testrun which should produce output files and plots in the `./res` folder: 

```Rscript nahrwhals.R --config conf/conf_default.txt```



NAHRwhals can also skip the initial search for sub-sequences, and call SVs directly on two regional fasta files by setting `compare_full_fastas = TRUE'. Run an example via: 

```Rscript nahrwhals.R --config conf/conf_fa2fa.txt```


# Basic Usage

To run your own data, exchange genome_x_fa (typically a reference, e.g. hg38), genome_y_fa (typically a genome assembly) and your coordinates of interest (seqname_x, start_x, end_x) in the config file.  

Parameters in the config file can be overwritten from the commandline, e.g. to add a gene annotation track to the plots:

```Rscript nahrwhals.R --config conf/config_examplerun.txt --params anntrack='testdata/assemblies/hg38_partial_genes.bed'```


# Report Errors

NAHRwhals is still in development, and errors in usage are to be expected. 
Please help improve the code by [reporting](https://github.com/WHops/nahrchainer/issues/new) issues you encounter.

# Citation

If you find NAHRwhals useful, please cite:

`https://www.biorxiv.org/content/10.1101/2023.03.09.531868v1`



# Correspondence


Please direct any correspondence to: Wolfram Höps (wolfram.hoeps@embl.de)