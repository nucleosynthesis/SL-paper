<sup>The SLtools repository accompanies the paper 
*The Simplified Likelihood framework*, **arXiv:1809.XXXXX**, by A. Buckley, M. Citron, S. Fichet, S. Kraml, W. Waltenberger and N. Wardle.</sup>


Reference code
-----------------------------------------------------------------------

We provide here a reference code, `simplike.py`, for the construction of the Simplified Likelihood (SL), 
eq. (2.3) in arXiv:1809.XXXXX. The expected number of backgrounds in terms of the combined nuisance parameters is parameterized as 
```math
  n_I = a_I + b_I \theta_I + c_I \theta_I^2
```
where the `$\theta_I$` are distributed as a multivariate normal. 

The reference code includes functions to calculate the coefficients $`a_I`$, $`b_I`$, $`c_I`$, and $`\rho_{IJ}`$ from provided moments $`m_{1,I}`$, $`m_{2,IJ}`$ and $`m_{3,I}`$, see eqs.(2.9)-(2.12) in the paper. It also includes an `SLParams` class which computes these and higher-level statistics such as profile likelihoods, log likelihood-ratios, and related limit-setting measures computed using observed and expected signal yields.

**Note:** the code has been written with reverse engineering and comprehensibility of the calculations explicitly in mind. While it computes likelihood statistics on a reasonable timescale, further (but less readable) optimisations can be added for production code.


Simplified likelihood demo
-----------------------------------------------------------------------

A demo of the construction of the simplified likelihood, and profiling as a function of a signal strength parameter, is given in `simplikedemo.py`. 

The first step is to read the inputs from which the simplified likelihood co-efficients are constructed, and the inputs to define the signal and the observed data. These inputs for the pseudo-search from **arXiv:1809.XXXXX** are contained in the python module `model-90_100000toys.py`. 

For searches which have provided these inputs in the form of HepData yaml tables, simple functions can be found in `convert_yaml.py` to read the simplified likelihood directly from the yaml tables. As an example, you can download the yaml files for the pseudo-search from [this HepData repository](https://www.hepdata.net/record/sandbox/1535641814 "HepData for pseudo-search") and set `USE_YAML_FILES=True` in `simplikedemo.py`. This will invoke the use of the conversion functions `yaml_table_to_pylist` and `yaml_multi_table_to_pylist` for the pseudo-search. 


Citation policy
-------------------------------------------------------------------------

When using the Simplified Likelihood framework or the reference SL code, please cite all the following papers:

   - *The Simplified Likelihood framework*, A. Buckley, M. Citron, S. Fichet, S. Kraml, W. Waltenberg, N. Wardle, arXiv:1809.XXXXX

   - *Simplified likelihood for the re-interpretation of public CMS results*, The CMS Collaboration, Tech. Rep. CMS-NOTE-2017-001

   - *Taming systematic uncertainties at the LHC with the central limit theorem*, S. Fichet, Nucl.Phys. B911 (2016) 623-637, arXiv:1603.03061
 
For convenience, the `SL_refs.bib` file in this repository provides these references in BibTeX format.

