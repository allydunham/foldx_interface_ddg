# FoldX Interface &Delta;&Delta;G Calculation

Calculate protein interface &Delta;&Delta;G values for a collection of proteins complexes given a list of mutations of interest for each complex.
FoldX's AnalyseComplex command is run on each WT complex, then the required mutated models are predicted with BuildModel and analysed again with AnalyseComplex to determine the difference in binding strength.

## Requirements

* FoldX 5
* Python 3
* Snakemake
* Numpy
* Pandas

## Usage

1. Clone this repo
2. Place your complex model PDB files as `data/complex/{complex id}/model.pdb`
3. Add an `individual_list` file to each complex directory, as specified by FoldX's [BuildModel](http://foldxsuite.crg.eu/command/BuildModel)
4. Run Snakemake

## Modifications

A number of simple modifications can be easily made:

* Generate all possible interface mutants based on the output of the WT interface analysis (see my [Mutfunc: SARS-CoV-2 pipeline](https://github.com/allydunham/mutfunc_sars_cov_2))
* Increase RAM requirements for individual complexes (see `get_mut_complex_ram` in `Snakefile`)
* Analyse interfaces other than between chain A and B (see `get_complex_tsv_files` in `Snakefile`)
