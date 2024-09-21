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
 
We recommend using mamba to install dependencies.

#### For Mac Silicon (M1, M2, M3) Users (Untested!): 
```
conda config --set channel_priority flexible
export CONDA_SUBDIR=osx-64
mamba env create --file env_nahrwhals_apple_silicon.yml
conda activate nahrwhals
unset CONDA_SUBDIR

# Install Julia manually:
# Step 1: download a julia .dmg file:
wget https://julialang-s3.julialang.org/bin/mac/aarch64/1.10/julia-1.10.5-macaarch64.dmg
# Step 2: install by double clicking and adding to your applications
# Step 3: update your $PATH in .zshrc or .profile to point to julia: 
export PATH=/Applications/Julia-1.10.app/Contents/Resources/julia/bin:$PATH
# Step 4: Install julia packages
julia -e 'using Pkg; Pkg.add("DelimitedFiles"); Pkg.add("ProgressMeter"); Pkg.add("ArgParse")'

```


#### For all other users:

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

(The latter is equivalent to running the following command:)
> nahrwhals(ref_fa = 'inst/extdata/assemblies/hg38_partial.fa', 
            asm_fa = 'inst/extdata/assemblies/assembly_partial.fa', 
            anntrack='inst/extdata/assemblies/hg38_partial_genes.bed',
            region='chr1_partial:1700000-3300000')
```

# Usage

1) Provide region=chr:start-end coordinates to genotype a single region
```
R
> library(nahrwhals)
> nahrwhals(ref_fa = 'ref.fa', 
            asm_fa = 'asm.fa,
            region = 'chr:start-end',
            outdir = 'res',
            threads = 8)
```

2) Provide a regions.bedfile to genotype multiple regions at once
```
R
> library(nahrwhals)
> nahrwhals(ref_fa = 'ref.fa', 
            asm_fa = 'asm.fa,
            regionfile = 'regions.bed',
            outdir = 'res',
            threads = 8)
```

3) Provide no region coordinates to invoke whole genome discovery mode. If ref.fa is human (or comparably complex), you must provide a blacklist bed file to skip centromeres and acrocentric chromsome arms. Tracks for hg38 and t2t are included in the package. Resources: ~1h using 16cpus, max-memory usage ~40Gb.


```
R
> library(nahrwhals)
> nahrwhals(ref_fa = 'ref.fa', 
            asm_fa = 'asm.fa,
            outdir = 'res',
            blacklist = system.file("extdata", "blacklists", "t2t_blacklist.bed", package = "nahrwhals"),
            threads = 8)
```



## Postprocessing 
Modes 2 and 3 automatically invoke post-processing of the nahrwhals_res.tsv file. This can be also done post-hoc on any res.tsv:
```
> nahrwhals::tsv_to_bed_regional_dominance('/your/nahrwhals.tsv', '/output/nahrwhals.bed')
```

# Output

Key results in the outdir folder are:

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




Figures and results from the NAHRwhals manuscript can be replicated by querying coordinates of interest (see supp. Tables and https://doi.org/10.5281/zenodo.7635935). Note the use of T2T as reference. 


# Run Modes and Required Parameters

The script supports three distinct run modes, each with its own set of required parameters. Below is an outline of what is needed for each mode:

## Common required Parameters

These parameters are applicable across all run modes unless otherwise specified.

| Variable name | Description                                 | Default value |
|---------------|---------------------------------------------|---------------|
| ref_fa        | Path to the reference genome FASTA file.    | -     |
| asm_fa        | Path to the assembly genome FASTA file for comparison. | - |
| outdir        | Path to desired output directory                   | './res' |


## Mode-dependent required Parameters
### Single-region Mode

| Variable name | Description                                 | Default value |
|---------------|---------------------------------------------|---------------|
| region        | reference coordinates to genotype                   | - |


### Multi-region Mode

| Variable name | Description                                 | Default value |
|---------------|---------------------------------------------|---------------|
| regionfile        | path to a bedfile with reference coords to genotype                 | - |

### Whole-genome Mode

| Variable name | Description                                 | Default value |
|---------------|---------------------------------------------|---------------|
| blacklist        | path to a blacklist bedfile of centromeres and acrocentric chr arms to skip. Required to avoid explosion of runtimes in human genomes.                  | - |


# Optional Parameters

### Optional Recommended Parameters

| Variable name           | Description                                                                                  | Default value |
|-------------------------|----------------------------------------------------------------------------------------------|---------------|
| ref_name            | Label for the reference genome used in outputs.                                              | 'Fasta_ref'     |
| asm_name            | Label for the assembly genome used in outputs.                                               | 'Fasta_asm'     |
| anntrack                | Include annotation tracks (e.g., genes) in dotplot. Specify as a 4-column bedfile (col 4: displayed name). | FALSE        |

### Optional Advanced Parameters

| Variable name           | Description                                                                                  | Default value |
|-------------------------|----------------------------------------------------------------------------------------------|---------------|
| depth                   | Depth of the BFS mutation search.                                                             | 3             |
| eval_th                 | Evaluation threshold as a percentage.                                                         | 98            |
| chunklen                | Sequence chunk length for alignment. | 'default' |
| minlen                  | Minimum length for considering alignments in segmentation.                                     | 'default'    |
| compression             | Minimum segment length.                                                                       | 'default'    |
| max_size_col_plus_rows  | Maximum size for the alignment matrix.                                                        | 250           |
| max_n_alns              | Maximum number of alignments for a single query sequence.                                     | 150           |
| self_plots              | Generates self-comparison plots for the assembly and reference genomes.                      | TRUE          |
| plot_only               | Skip BFS/genotyping and only create dotplots.                                                 | FALSE         |
| use_paf_library         | Use an external PAF library for alignments.                                                   | FALSE         |
| conversionpaf_link      | If `use_paf_library` is TRUE, provide link to whole-genome PAF here.                          | FALSE         |
| maxdup                  | BFS: max number of duplications per chain.                                                    | 2             |
| init_width              | BFS: follow up only the best-scoring `init_width` nodes per depth.                            | 1000          |
| region_maxlen           | Maximum length of a window that can be analyzed. Larger windows will be split into overlapping fragments. | 5000000   |
| threads                 | In whole genome / multi-region: regions analysed simultaneously. Single-region: cores for minimap2                                                       | 1             |
| asm_fa_mmi         | Path to pre-indexed assembly genome with minimap2. Useful to avoid re-calculating.                                            | 'default'     |


# Report Errors

Please help improve the code by [reporting](https://github.com/WHops/nahrchainer/issues/new) issues you encounter.

# Citation

For more information about NAHRwhals, check out our [preprint](https://www.biorxiv.org/content/10.1101/2023.03.09.531868v1)!


# Correspondence

Please direct any correspondence to: Wolfram Höps (wolfram.hoeps@gmail.com)
