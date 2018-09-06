The SLtools repository accompanies the paper 
*The Simplified Likelihood framework*, **arXiv:1809.XXXXX**, by A. Buckley (Un. of Glasgow), M. Citron (Un. of California), S. Fichet (Caltech and ICTP-SAIFR&IFT-UNESP), S. Kraml (LPSC/CNRS/IN2P3), W. Waltenberg (Un. of Vienna) and N. Wardle (Imperial College).



Calculation of next-to-leading order simplified likelihood coefficients
-----------------------------------------------------------------------

SL backgrounds $`b_I`$ for signal region $`I`$ are parametrised as
```math
n_{b,I} = a_I + b_I \theta_I + c_I \, \theta_I^2
```
where the $`\theta_I`$ are distributed as a multivariate normal.

This reference code has been written with reverse engineering and
comprehensibility of the calculations explicitly in mind. While it computes
likelihood statistics on a reasonable timescale, further (but less readable)
optimisations can be added for production code.

This package includes functions to calculate the SL $`a_I`$, $`b_I`$, $`c_I`$, and
$`\rho_{IJ}`$ coefficients from provided moments $`m_{1,I}`$, $`m_{2,IJ}`$ and
$`m_{3,I}`$; and an `SLParams` class which computes these and higher-level
statistics such as profile likelihoods, log likelihood-ratios, and related
limit-setting measures computed using observed and expected signal yields.

Simplified likelihood demo
-----------------------------------------------------------------------

A demo of the construction of the simplified likelihood, and profiling as a function of a signal strength parameter, is given in `simplikedemo.py`. 

The first step is to read the inputs from which the simplified likelihood co-efficients are constructed, and the inputs to define the signal and the observed data. These inputs for the pseudo-search from **arXiv:1809.XXXXX** are contained in the python module `model-90_100000toys.py`. 

For searches which have provided these inputs in the form of HepData yaml tables, simple functions can be found in `convert_yaml.py` to read the simplified likelihood directly from the yaml tables. As an example, you can download the yaml files for the pseudo-search from [HepData](https://www.hepdata.net/record/sandbox/1535641814 "HepData") and set `USE_YAML_FILES=True` in `simplikedemo.py`. This will invoke the use of the conversion functions `yaml_table_to_pylist` and `yaml_multi_table_to_pylist` for the pseudo-search. 

Citation Policy
-------------------------------------------------------------------------

When using the Simplified Likelihood framework or the reference SL code, please refer to the three following papers:

*The Simplified Likelihood framework*, A. Buckley, M. Citron, S. Fichet, S. Kraml, W. Waltenberg, N. Wardle, arXiv:1809.XXXXX

*Simplified likelihood for the re-interpretation of public CMS results*, The CMS Collaboration, Tech. Rep. CMS-NOTE-2017-001

*Taming systematic uncertainties at the LHC with the central limit theorem*, S. Fichet, Nucl.Phys. B911 (2016) 623-637, arXiv:1603.03061
 
