<img src="https://github.com/WHops/NAHRwhals/blob/main/NAHRwhals.png?raw=true">


# NAHRwhals (NAHR-directed Workflow for catcHing seriAL Structural Variations)
An R package providing tools for simulation, detection and visualization of complex, NAHR-driven rearrangements. This includes:
		- Liftover of coordinates between arbitrary human DNA assemblies
		- Accurate sequence alignments of multi-MB DNA sequences
		- Dotplot visualizations 
		- BitLocus - condensed dotplots
		- A tree-based caller for complex, nested NAHR-mediated rearrangements. 

## Installation & Usage

[Mon, 13th March 2023] NAHRwhals is currently undergoing code streamlining and will be fully available in few days. Thanks for your patience. 

## Report Errors

NAHRwhals is still in development, and errors in usage are to be expected. 
Please help improve the code by [reporting](https://github.com/WHops/nahrchainer/issues/new) issues you encounter.

## Correspondence

Please direct any correspondence to: 
Wolfram Höps
wolfram.hoeps@embl.de


## Example runs

# Two assemblies and coordinates on one of them
Rscript nahrwhals.R --params genome_x_fa=testdata/assemblies/hg38_partial.fa genome_y_fa=testdata/assemblies/assembly_partial.fa seqname_x=chr1_partial start_x=1700000 end_x=3300000 hltrack=testdata/assemblies/hg38_partial_genes.bed anntrack=testdata/assemblies/hg38_partial_genes.bed



# Two fastas with each other
Rscript nahrwhals.R --params genome_x_fa=testdata/extracted_fastas/sequence1.fa genome_y_fa=testdata/extracted_fastas/sequence2.fa compare_full_fastas=T


