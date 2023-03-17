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
- [minimap2](https://github.com/lh3/minimap2) >= 2.24
- [bedtools](https://bedtools.readthedocs.io/en/latest/content/quick-start.html) <= 2.30.0
- [gawk](ttps://formulae.brew.sh/formula/gawk) >= 5.2.0
- [R](https://www.r-project.org/) >= 4.1.0
#  Installation


### (1) Install dependencies


 
NAHRwhals has a number of dependencies, including [minimap2](https://github.com/lh3/minimap2), [bedtools](https://bedtools.readthedocs.io/en/latest/content/quick-start.html), [gawk](https://formulae.brew.sh/formula/gawk) and a number of R packages. The easiest way to install all is by using the provided conda/mamba environment. Conda takes a while to execute this, so we recommend using mamba, which is an accelerated version. 

```
conda config --set channel_priority flexible
mamba env create --file env_nahrwhals.yml
conda activate nahrwhals
```
### (2) Clone the NAHRwhals repository


```
git clone https://github.com/WHops/NAHRwhals.git
cd NAHRwhals
```

Depending on git version, this folder may appear empty. In this case, switch to the main branch with `git checkout main`.


### (3) Install 

You can install with the following command, which simply calls devtools::install. 

`Rscript install_package.R`

Confirm successful installation using

```
R
> library(nahrwhals)
> quit()
```

### (4) [OPTIONAL] NOT NEEDED WITH CONDA/MAMBA: specify minimap2 / bedtools paths; 
In case minimap2 or bedtools are not part of your $PATH (i.e. can not be called from the commandline via `minimap2` and `bedtools`), specify your paths in conf/conf_default.txt: (otherwise, leave them as 'default')

```
minimap2_bin = '/your/path/to/minimap2'
bedtools_bin = '/your/path/to/bedtools'
```



#  Test & Example runs

To confirm that NAHRwhals has been correctly installed, run a testrun which should produce output files and plots in the `./res` folder: 

```Rscript nahrwhals.R --config conf/conf_default.txt```



NAHRwhals can also skip the initial search for sub-sequences, and call SVs directly on two regional fasta files by setting `compare_full_fastas = TRUE'. Run an example via: 

```Rscript nahrwhals.R --config conf/conf_fa2fa.txt```


# Use cases

To run your own data, exchange genome_x_fa (typically a reference, e.g. hg38), genome_y_fa (typically a genome assembly) and your coordinates of interest (seqname_x, start_x, end_x) in the config file.  

For convenience, parameters in the config file can also be overwritten from the commandline. We can use this e.g. to add a gene annotation track to our plots:

```Rscript nahrwhals.R --config conf/conf_default.txt --params anntrack='testdata/assemblies/hg38_partial_genes.bed'```

# Parameters

NAHRwhals comes with a variety of parameters specified in a config file and/or via commandline. Default values are denoted in [brackets].

## Paths and input/output

- **minimap2_bin** link to minimap2 binary. If minimap2 is in $PATH, keep value 'default'. [default]
  
- **bedtools_bin** link to bedtools binary. If bedtools is in $PATH, keep value 'default'. [default]
- **genome_x_fa** Multi- or single-contig genome sequence in fasta format, typically a genome assembly used as reference, such as hg38 or T2T. ['testdata/assemblies/hg38_partial.fa']
- **genome_y_fa** Multi- or single-contig genome sequence in fasta format, typically a genome assembly to be analysed, such as a de-novo assembled haplotype. ['testdata/assemblies/assembly_partial.fa']
- **genome_y_fa_mmi** Link to minimap2 index mmi file of genome_y_fa, which is needed for fast alignments. If the mmi does not exist, NAHRwhals will invoke minimap2 to creat it at the link location. By default, it is created in the same directly as genome_y_fa. [default]
- **anntrack** You can supply a .bed file (seqname, start, end, featurename) for annotations to appear in the dotplot. Typically used to indicate e.g. genes in a dotplot. [F]
- **logfile** Output file to which sSV calls will be written ['res/unittest.tsv']

## Basic parameters: region of interest

- **compare_full_fastas** If TRUE, seqname_x, start_x and end_x are ignored, and instead the whole genome_x_fa is aligned to the whole genome_y_fa. Use this if you have two regional fasta files that you want to compare. Do NOT use this when dealing with whole genome assemblies. [F]
- **seqname_x** Contigname of region of interest on genome_x. ['chr1_partial']
- **start_x** Start of region of interest on seqname_x. (Region of interest should typically be in a range of ~20kbp to ~5-10Mbp and not contain centromeric regions, which drive up computation time significantly.) [1700000]
- **end_x** End of region of interest on seqname_x. [3300000]
  
## Advanced parameters

- **samplename_x** Set a name for the genome_x sequence (e.g., samplename) to appear in the results. ['Fasta_x']
- **samplename_y** Set a name for the genome_y sequence (e.g., samplename) to appear in the results. ['Fasta_y']
- **plot_only** If TRUE, NAHRwhals creates only dotplots but does not attempt SV calling. [F]
- **self_plots** If TRUE, NAHRwhals automatically creates a self-dotplot (aligning sequence to itself) of both the reference and homologous region of interest. Will be saved in res/chrX-start-end/self. [T]
- **plot_xy_segmented** If TRUE, NAHRwhals outputs a plot indicating obtained segmentation. [T]
- **eval_th** Set the threshold (in percent) for when NAHRwhals considers two sequences equal. Lower thresholds lead to missed SVs, higher thresholds to more 'unexplained' loci. Consider lowering to e.g. 95% for cross-species analyses. [98]
- **depth** Maximum number of consecutive SVs to model. [3]
- **chunklen** Lenght of sequence chunks which are separately aligned to each other. By default, this is a function of the total sequence length and typically 1000 bp or 10.000 bp. It is typically not necessary to touch this parameter, as pairwise alignments are robust. [auto]
- **minlen** Minimum alignment length to consider when searching for sSVs. By default, this is a function of the total sequence length and typically 1000 bp or 10.000 bp. Lower values can lead to finer alignments but high computation times and false positive calls. [auto]
- **compression** Minimum segment length to consider when segmenting a dotplot. By default, this is a function of the total sequence length and typically 1000 bp or 10.000 bp. Lower values can lead to finer segmented dotplots but high computation times and false positive calls. [auto]
- **max_size_col_plus_rows** Maximum acceptable size for the segmented dotplot (rows + columns). Is set to prevent computation / memory overflow. If this value is overstepped, minlen and chunklen are doubled until the resulting segmented dotplot its in the maximum dimensions. [250]
- **max_n_alns** Maximum number of individual alignments that can participate in dotplot segmentation. If this value is overstepped, minlen and chunklen are doubled until the resulting segmented dotplot its in the maximum dimensions. [150]



# Report Errors

NAHRwhals is still in development, and errors in usage are to be expected. 
Please help improve the code by [reporting](https://github.com/WHops/nahrchainer/issues/new) issues you encounter.

# Citation

If you find NAHRwhals useful, please cite:

`https://www.biorxiv.org/content/10.1101/2023.03.09.531868v1`



# Correspondence


Please direct any correspondence to: Wolfram Höps (wolfram.hoeps@embl.de)