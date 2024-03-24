<img src="https://github.com/WHops/NAHRwhals/blob/main/NAHRwhals.png?raw=true">


# NAHRwhals (NAHR-directed Workflow for catcHing seriAL Structural Variations)
NAHRwhals is an R package providing tools for visualization and automatic detection of complex, NAHR-driven rearrangements (few kbp to multiple Mbp) using genome assemblies. Modules include:

- Liftover of coordinates between arbitrary human DNA assemblies
- Accurate sequence alignments of multi-MB DNA sequences
- Dotplot visualizations 
- Segmented Dotplots
- A tree-based caller for complex, nested NAHR-mediated rearrangements. 



# Installation 

### (1) Clone the NAHRwhals repository


```
git clone https://github.com/WHops/NAHRwhals.git
cd NAHRwhals
```
### (2) Install dependencies
 
We recommend using mamba to install dependencies:

```
conda config --set channel_priority flexible
mamba env create --file env_nahrwhals.yml
conda activate nahrwhals

julia -e 'using Pkg; Pkg.add("DelimitedFiles"); Pkg.add("ProgressMeter"); Pkg.add("ArgParse")'
```


### (3) Install NAHRwhals

Install NAHRwhals using:

```
Rscript install_package.R
```


### (4) Test your installation

Confirm the installation with a testrun (producing output files and plots in the `./res` folder). 
```
R
> library(nahrwhals)
> nahrwhals(testrun_std=T)
```

# Usage

1) Provide region=chr:start-end coordinates to genotype single region
```
R
> library(nahrwhals)
> nahrwhals(ref_fa = 'ref.fa', 
            alt_fa = 'alt.fa,
            region = 'chr:start-end',
            outdir = 'res',
            minimap_cores = 8)
```

2) Provide regions_file to coordinate multiple regions at once
```
R
> library(nahrwhals)
> nahrwhals(ref_fa = 'ref.fa', 
            alt_fa = 'alt.fa,
            regions_file = 'regions.bed',
            outdir = 'res',
            threads = 8)
```

3) Provide no region coordinates to invoke whole genome discovery mode
```
R
> library(nahrwhals)
> nahrwhals(ref_fa = 'ref.fa', 
            alt_fa = 'alt.fa,
            outdir = 'res',
            threads = 8)
```



# Output

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
<td><img src="https://github.com/WHops/NAHRwhals/blob/main/inst/extdata/output_png/Fasta_x_Fasta_y_all/Fasta_x_Fasta_y_all-1.png?raw=true" alt="Page 1"></td>
<td><img src="https://github.com/WHops/NAHRwhals/blob/main/inst/extdata/output_png/Fasta_x_Fasta_y_all/Fasta_x_Fasta_y_all-2.png?raw=true" alt="Page 2"></td>
<td><img src="https://github.com/WHops/NAHRwhals/blob/main/inst/extdata/output_png/Fasta_x_Fasta_y_all/Fasta_x_Fasta_y_all-3.png?raw=true" alt="Page 3"></td>
<td><img src="https://github.com/WHops/NAHRwhals/blob/main/inst/extdata/output_png/Fasta_x_Fasta_y_all/Fasta_x_Fasta_y_all-4.png?raw=true" alt="Page 4"></td>
<td><img src="https://github.com/WHops/NAHRwhals/blob/main/inst/extdata/output_png/Fasta_x_Fasta_y_all/Fasta_x_Fasta_y_all-5.png?raw=true" alt="Page 5"></td>
<td><img src="https://github.com/WHops/NAHRwhals/blob/main/inst/extdata/output_png/Fasta_x_Fasta_y_all/Fasta_x_Fasta_y_all-6.png?raw=true" alt="Page 6"></td>
<td><img src="https://github.com/WHops/NAHRwhals/blob/main/inst/extdata/output_png/Fasta_x_Fasta_y_all/Fasta_x_Fasta_y_all-7.png?raw=true" alt="Page 6"></td>

</tr>
</table>

- 2.1. Self-dotplot of the selected region on sequence x
- 2.2. Self-dotplot of the corresponding homologous region on sequence y
- 2.3. Pairwise dotplot of the two regions
- 2.4. Visualization of the obtained segments on sequence x
- 2.5. Pariwise dotplot colored on x axis by obtained segments
- 2.6. Segmented pairwise alignment
- 2.7. Segmented pairwise alignment AFTER applying the top-scoring mutation




Figures and results from the NAHRwhals manuscript can be replicated by querying coordinates of interest (see supp. Tables) on a reference (parameter: genome_x_fa) against assembly fastas (parameter: genome_y_fa) (see https://doi.org/10.5281/zenodo.7635935). Note the use of T2T as reference. 

# Parameters

## Required

| Variable name | Description | Default value|
|-|-|--|
| genome_x_fa | Multi- or single-contig genome sequence in fasta format, typically a genome assembly used as reference, such as hg38 or T2T.| - |
| genome_y_fa | Multi- or single-contig genome sequence in fasta format, typically a genome assembly to be analyzed, such as a de-novo assembled haplotype.| - |
| seqname_x| Contigname of the region of interest on genome_x. | - |
| start_x| Start of the region of interest on seqname_x. (Region of interest should typically be in a range of ~20kbp to ~5-10Mbp and not contain centromeric regions, which drive up computation time significantly.)| - |
| end_x| End of the region of interest on seqname_x. | - |




Based on your updated function parameters and their descriptions, here's a revised README section:

# Parameters

## Required

| Variable name | Description | Default value|
|-|-|--|
| ref_fa | Path to the reference genome FASTA file. | - |
| asm_fa | Path to the assembly genome FASTA file for comparison. | - |
| outdir | Basedirectory to store output files and plots. | - |

## Recommended parameters

| Variable name | Description | Default value |
|-------|---|-------|
| samplename_x | Reference genome label for outputs. | 'Fasta_x' |
| samplename_y | Assembly genome label for outputs. | 'Fasta_y' |
| region | Region to genotype in the format "chr:start-end". Mutually exclusive with 'regionfile' | [optional] |
| regionfile | Path to a file listing regions to genotype in bed format. Mutually exclusive with 'region' | [optional] |
| anntrack | Includes annotation tracks (e.g., genes) in dotplot. Specify as 4-column bedfile (col 4: displayed name). | FALSE |

## Advanced parameters

| Variable name | Description | Default value |
|---|-------|-------|
| depth | Depth of the BFS mutation search. | 3 |
| eval_th | Evaluation threshold as a percentage. | 98 |
| chunklen | Sequence chunk length alignment, "default" means automatic determination based on sequence length. | auto-detection |
| minlen | Disregard alignments shorter than this from the segmentation. | 350 |
| compression | Minimum segment length. | auto-detection |
| max_size_col_plus_rows | Maximum size for the alignment matrix. | 250 |
| max_n_alns | Maximum number of alignments for a single query sequence. | 150 |
| self_plots | Generates self-comparison plots for the assembly and ref. | TRUE |
| plot_only | Skip BFS/genotyping, just create dotplots. | FALSE |
| use_paf_library | Use an external PAF library for alignments. | FALSE |
| conversionpaf_link | If use_paf_library: provide link to whole-genome paf here. | FALSE |
| maxdup | BFS: max number of duplications per chain. | 2 |
| init_width | BFS: follow up only init_width best-scoring nodes per depth. | 1000 |
| region_maxlen | Maximum length of a window that can be analysed. Larger windows will be split into overlapping fragments. | 5000000 |
| testrun_std | Test the NAHRwhals installation. | FALSE |
| threads | Number of windows to parallel genotype. | 1 |
| minimap_cores_per_thread | Number of cores to give to minimap PER THREAD. | 1 |
| genome_y_fa_mmi | Path to pre-indexed assembly genome with minimap2. | "default" |


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



# Report Errors

Please help improve the code by [reporting](https://github.com/WHops/nahrchainer/issues/new) issues you encounter.

# Citation

For more information on how NAHRwhals, check out our [preprint](https://www.biorxiv.org/content/10.1101/2023.03.09.531868v1)!


# Correspondence

Please direct any correspondence to: Wolfram Höps (wolfram.hoeps@gmail.com)
