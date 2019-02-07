# R Code and data sets for Birks & Simpson (2013)

This repo contains the data and R code used in the analyses presented in Birks &amp; Simpson (2013), *'Diatoms and pH reconstruction' (1990) revisited*, published in the Journal of Paleolimnology.

The `.R` scripts do the following tasks:

1. data processing
2. run the transfer function analyses
3. collect the results and prepare figures

## Data processing

* `load_data.R` loads the SWAP data set. It was intended to be called from each of the analysis scripts mentioned below.
* `load_uk_data.R` loads the UK training set. It was intended to be called from each of the analysis scripts mentioned below.

## Analyses

There are four main transfer function methods applied in the paper &mdash; MAT, WA, WA-PLS, and GLR &mdash; and the code for each method is contained in a separate script. These can be run in any order, although I used separate R sessions as there are conflicts between the *analogue* and *rioja* packages in some functions names. Each script produces a number of objects that are used later, and which are serialized to an `.rds` file, which is saved to disk at the edn of each of the analysis scripts.

Analysis scripts:

* Gaussian logit regression (GLR): `glr_main_data_analysis.R`
* Modern analogue technique (MAT): `mat_main_data_analysis.R`
* Weight averaging (WA): `wa_main_data_analysis.R`
* Weighted averaging partial least squares (WA-PLS): `wapls_main_data_analysis.R`

The wrapper/helper functions in `paper_fun.R` are loaded in each of the analysis scripts.

## Processing results and figures

Each of the `.rds` files containing the output from the analyses is processed with the `process_results.R` script. This produces the figures presented in the paper and CSV files of tables of results.

## References

Birks, H. J. B., and G. L. Simpson. 2013. “Diatoms and pH reconstruction” (1990) revisited. *J. Paleolimnol.* **49**: 363–371. doi:[10.1007/s10933-013-9697-7](https://doi.org/10.1007/s10933-013-9697-7)
