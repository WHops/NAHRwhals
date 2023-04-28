<img src="https://github.com/WHops/NAHRwhals/blob/main/NAHRwhals.png?raw=true">


# NAHRwhals (NAHR-directed Workflow for catcHing seriAL Structural Variations)
NAHRwhals is an R package providing tools for visualization and automatic detection of complex, NAHR-driven rearrangements (few kbp to multiple Mbp) using genome assemblies. Modules include:

- Liftover of coordinates between arbitrary human DNA assemblies
- Accurate sequence alignments of multi-MB DNA sequences
- Dotplot visualizations 
- Segmented Dotplots
- A tree-based caller for complex, nested NAHR-mediated rearrangements. 



#Installation (total time: ca 10 minutes)

### (1) Clone the NAHRwhals repository


```
git clone https://github.com/WHops/NAHRwhals.git
cd NAHRwhals; git checkout main
```
### (2) Install dependencies
 
NAHRwhals has a number of dependencies, including [minimap2](https://github.com/lh3/minimap2), [bedtools](https://bedtools.readthedocs.io/en/latest/content/quick-start.html), [gawk](https://formulae.brew.sh/formula/gawk) and a number of R packages. The easiest way to install all is by using the provided conda/mamba environment. Conda takes a while to execute this, so we recommend using mamba, which is an accelerated version. Installation/Download time: 5-10 minutes. 

```
conda config --set channel_priority flexible
mamba env create --file env_nahrwhals.yml
conda activate nahrwhals
```


### (3) Install NAHRwhals

You can finally install NAHRwhals with the following command. 

```
Rscript install_package.R
```



# Test & Example runs (runtime: <1 min)

To confirm that NAHRwhals has been correctly installed, run a testrun inside R which will produce output files and plots in the `./res` folder: 

```
R
> library(nahrwhals)
> nahrwhals(testrun_std=T)
```




Key results in the .res folder are:

## 1. **res/res.tsv**: The SV call output file.



| seqname| start | end | sample | res_ref | res_max | mut_max|
| ---- | ------- | ------- | ------| ------- | ------- | ----|
| chr1_partial | 1700000 | 3300000 | Fasta_x_Fasta_y | 71.378| 99.595| 4_12_inv+11_18_del|

In the example provided, the input sequence x and its counterpart on y yielded a 71.3% correspondence in reference state and a 99.6% correspondence after applying the highest-scoring mutation, 'mut_max' (Inversion between blocks 4 and 12, followed by deletion of blocks 11 to 18). The file contains additional columns used mainly for QC. Those are explained at the end of the README. 

## 2. **res/chr1_partial-1700000-3300000/Fasta_x_Fasta_y_all.pdf**: 

This PDF contains seven plots that document the NAHRwhals run:

<table>
<tr>
<td><img src="https://github.com/WHops/NAHRwhals/blob/package_cleaning/inst/extdata/output_png/Fasta_x_Fasta_y_all/Fasta_x_Fasta_y_all-1.png?raw=true" alt="Page 1"></td>
<td><img src="https://github.com/WHops/NAHRwhals/blob/package_cleaning/inst/extdata/output_png/Fasta_x_Fasta_y_all/Fasta_x_Fasta_y_all-2.png?raw=true" alt="Page 2"></td>
<td><img src="https://github.com/WHops/NAHRwhals/blob/package_cleaning/inst/extdata/output_png/Fasta_x_Fasta_y_all/Fasta_x_Fasta_y_all-3.png?raw=true" alt="Page 3"></td>
<td><img src="https://github.com/WHops/NAHRwhals/blob/package_cleaning/inst/extdata/output_png/Fasta_x_Fasta_y_all/Fasta_x_Fasta_y_all-4.png?raw=true" alt="Page 4"></td>
<td><img src="https://github.com/WHops/NAHRwhals/blob/package_cleaning/inst/extdata/output_png/Fasta_x_Fasta_y_all/Fasta_x_Fasta_y_all-5.png?raw=true" alt="Page 5"></td>
<td><img src="https://github.com/WHops/NAHRwhals/blob/package_cleaning/inst/extdata/output_png/Fasta_x_Fasta_y_all/Fasta_x_Fasta_y_all-6.png?raw=true" alt="Page 6"></td>
<td><img src="https://github.com/WHops/NAHRwhals/blob/package_cleaning/inst/extdata/output_png/Fasta_x_Fasta_y_all/Fasta_x_Fasta_y_all-7.png?raw=true" alt="Page 6"></td>

</tr>
</table>

- 2.1. Self-dotplot of the selected region on sequence x
- 2.2. Self-dotplot of the corresponding homologous region on sequence y
- 2.3. Pairwise dotplot of the two regions
- 2.4. Visualization of the obtained segments on sequence x
- 2.5. Pariwise dotplot colored on x axis by obtained segments
- 2.6. Segmented pairwise alignment
- 2.7. Segmented pairwise alignment AFTER applying the top-scoring mutation


# Quick start: your own data (runtime: <1 min)

To run your own data, exchange genome_x_fa (typically a reference, e.g. hg38), genome_y_fa (typically a genome assembly) and your coordinates of interest (seqname_x, start_x, end_x) in the config file. An example run which includes a gene annotation track (optional) could look like this: 

```
nahrwhals(genome_x_fa = system.file("extdata/assemblies", "hg38_partial.fa", package="nahrwhals"),
genome_y_fa = system.file("extdata/assemblies", "assembly_partial.fa", package="nahrwhals"),
seqname_x = 'chr1_partial',
start_x = 1700000,
end_x = 3300000,
anntrack= system.file("extdata/assemblies", "hg38_partial_genes.bed", package="nahrwhals"),
samplename_x = 'Fasta_x',
samplename_y = 'Fasta_y'
)
```

Resulting in pairwise alignments plots with gene annotations overlaid on the top: 

<img src="https://github.com/WHops/NAHRwhals/blob/main/testdata/output_png/examplerun_output_hg38_T2T.png?raw=true" width="500" height="500">


Figures and results from the NAHRwhals manuscript can be replicated by querying coordinates of interest (see supp. Tables) on a reference (parameter: genome_x_fa) against assembly fastas (parameter: genome_y_fa) (see https://doi.org/10.5281/zenodo.7635935). Note that in the majority of cases, T2T is used as reference.

# Parameters

## Required

| Variable name | Description | Default value|
|-|-|--|
| genome_x_fa | Multi- or single-contig genome sequence in fasta format, typically a genome assembly used as reference, such as hg38 or T2T.| - |
| genome_y_fa | Multi- or single-contig genome sequence in fasta format, typically a genome assembly to be analyzed, such as a de-novo assembled haplotype.| - |
| seqname_x| Contigname of the region of interest on genome_x. | - |
| start_x| Start of the region of interest on seqname_x. (Region of interest should typically be in a range of ~20kbp to ~5-10Mbp and not contain centromeric regions, which drive up computation time significantly.)| - |
| end_x| End of the region of interest on seqname_x. | - |




## Recommended parameters


| Variable name | Description | Default value |
|-------|---|-------|
| samplename_x| Set a name for the genome_x sequence (e.g., samplename) to appear in the results.| 'Fasta_x' |
| samplename_y| Set a name for the genome_y sequence (e.g., samplename) to appear in the results.| 'Fasta_y' |
| anntrack| You can supply a .bed file (seqname, start, end, featurename) for annotations to appear in the dotplot. Typically used to indicate e.g. genes in a dotplot.| F|
| logfile | Output file to which sSV calls will be written.| 'res/res.tsv' |


## Advanced parameters


| Variable name | Description | Default value |
|---|-------|-------|
| genome_y_fa_mmi | Link to minimap2 index mmi file of genome_y_fa, which is needed for fast alignments. If the mmi does not exist, NAHRwhals will invoke minimap2 to create it at the link location. | 'default' (By default, it is searched/created in the same directory as genome_y_fa). |
| plot_only | If TRUE, NAHRwhals creates only dotplots but does not attempt SV calling.| F |
| self_plots| If TRUE, NAHRwhals automatically creates a self-dotplot (aligning sequence to itself) of both the reference and homologous region of interest. Will be saved in res/chrX-start-end/self. | T |
| plot_xy_segmented | If TRUE, NAHRwhals outputs a plot indicating obtained segmentation. | T |
| eval_th | Set the threshold (in percent) for when NAHRwhals considers two sequences equal. Lower thresholds lead to missed SVs, higher thresholds to more 'unexplained' loci. Consider lowering to e.g. 95% for cross-species analyses. | 98|
| depth | Maximum number of consecutive SVs to model. | 3 |
| chunklen| Length of sequence chunks which are separately aligned to each other. By default, this is a function of the total sequence length and typically 1000 bp or 10.000 bp. It is typically not necessary to touch this parameter, as pairwise alignments are robust.| 'auto'|
| minlen| Minimum alignment length to consider when searching for sSVs. By default, this is a function of the total sequence length and typically 1000 bp or 10.000 bp. Lower values can lead to finer alignments but high computation times and false positive calls.| 'auto'|
| compression | Minimum segment length to consider when segmenting a dotplot. By default, this is a function of the total sequence length and typically 1000 bp or 10.000 bp. Lower values can lead to finer segmented dotplots but high computation times and false positive calls.| 'auto'|
| max_size_col_plus_rows| Maximum acceptable size for the segmented dotplot (rows + columns). Is set to prevent computation / memory overflow. If this value is overstepped, minlen and chunklen are doubled until the resulting segmented dotplot is in the maximum dimensions. | 250 |
| max_n_alns| Maximum number of individual alignments that can participate in dotplot segmentation. If this value is overstepped, minlen and chunklen are doubled until the resulting segmented dotplot is in the maximum dimensions. | 150 |
| testrun_std | Runs a testrun with example data. | F |
| testrun_fatofa| Runs a testrun with example data, but with the mode activated that compares full fastas to each other.| F |
| compare_full_fastas| If TRUE, seqname_x, start_x, and end_x are ignored, and instead the whole genome_x_fa is aligned to the whole genome_y_fa. Use this if you have two regional fasta files that you want to compare. Do NOT use this when dealing with whole genome assemblies. | 'F' |
 minimap2_bin| Link to the minimap2 binary. If minimap2 is in $PATH, keep the value 'default'. | default|
| bedtools_bin| Link to the bedtools binary. If bedtools is in $PATH, keep the value 'default'. | default|



# Output columns

The main output from a nahrwhals run is found in res/res.tsv. Here is a description of all columsn in this file: 


| Column Name| Description |
|----|---|
| seqname| The name of the sequence. This could be a chromosome, contig, or any other identifier for a sequence. |
| start| The start position of the sequence. |
| end| The end position of the sequence. |
| sample | The name of the run, denoted as namex_namey.|
| width_orig | The width of the input sequence. |
| xpad | The amount of padding added to the sequence (obsolete).|
| res_ref| Pairwise alignment concordance of the sequences on x and y in native state.|
| res_max| Pairwise alignment concordance after applying the highest-scoring mutation.|
| n_res_max| The number of mutation sequences obtaining res_max. If >1, there are alternative paths.|
| mut_max| The highest-scoring mutation sequence. |
| mut_simulated| The number of simulated mutations.|
| mut_tested | The number of tested mutations (=mut_simulated minus rejected). |
| search_depth | The depth of the search algorithm used.|
| grid_compression | The amount of compression used for the grid, in bp. |
| exceeds_x| QC: Whether or not the sequence exceeds the x-axis.|
| exceeds_y| QC: Whether or not the sequence exceeds the y-axis.|
| grid_inconsistency | QC: Whether or not there is inconsistency in the grid. |
| flip_unsure| QC: Whether or not the algorithm was unsure about flipping the sequence.|
| cluttered_boundaries | QC: Whether or not the sequence has cluttered boundaries.|
| mut1_start | The start position of the first mutation.|
| mut1_end | The end position of the first mutation.|
| mut1_pos_pm| The confidence interval of both breakpoints (plus/minus). |
| mut1_len | The length of the first mutation.|
| mut1_len_pm| The confidence interval of SV length (plus/minus). |
| mut2_len | The length of the second mutation. |
| mut2_len_pm| The confidence interval of SV length (plus/minus).|
| mut3_len | The length of the third mutation.|
| mut3_len_pm| The confidence interval of SV length (plus/minus). |


# NOTES


- In case minimap2 or bedtools are not part of your $PATH (i.e. can not be called from the commandline via `minimap2` and `bedtools`), specify your paths when running nahrwhals:

```
library(nahrwhals)
nahrwhals( ... , minimap2_bin = '/path/to/your/minimap2', bedtools_bin = 'path/to/your/bedtools')
```

- NAHRwhals can also skip the initial search for sub-sequences, and call SVs directly on two regional fasta files by setting `compare_full_fastas = TRUE'. Run an example via: 

```
> nahrwhals(testrun_fullfa=T)
```

# Report Errors

NAHRwhals is still in development, and errors in usage are to be expected. 
Please help improve the code by [reporting](https://github.com/WHops/nahrchainer/issues/new) issues you encounter.

# Citation

If you find NAHRwhals useful, please cite:

`https://www.biorxiv.org/content/10.1101/2023.03.09.531868v1`



# Correspondence


Please direct any correspondence to: Wolfram Höps (wolfram.hoeps@embl.de)
