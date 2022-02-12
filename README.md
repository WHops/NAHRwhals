<img src="https://github.com/WHops/nahrchainer/blob/main/ntk_logo-01.png?raw=true">

=========================================================================

# nahrToolkiT (NTK)
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

3. run the installation
    `cd nahrchainer`
    `Rscript install_package.R`
    
4. Verify successful installation
`R`
` library(nahrtoolkit)`
`confirm_loaded_nahr()`
`>The NAHRtoolkit is loaded and ready to go.`


## Usage

Please refer to the vignette (/vignettes/) for example use cases of NTK. 

## Report Errors

Nahrtoolkit is still in early stages of development, and errors in usage are to be expected. 
Please help improve the code by [reporting](https://github.com/WHops/nahrchainer/issues/new) issues you encounter.

## Correspondence

Please direct any correspondence to: 
Wolfram Höps
wolfram.hoeps@embl.de

