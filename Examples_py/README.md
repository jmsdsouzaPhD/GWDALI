# List of Examples

* **1_get_hphx.py**
    - Usage of the function **get_hphx()** which returns GW polarizations (plus/cross);
* **2_get_strains.py**
    - Usage of the function **get_strain()** which returns the detectors outputs (linear combination of the polarizations);
* **3_get_derivatives.py**
    - Usage of the function **get_derivatives()** which returns signal derivatives according to the desired order of derivation and differentiation method (numerical or automatic differentiation);
* **4_get_dali-tensors.py**
    - Usage of the function **get_tensors()** which returns the *DALI* tensors (Fisher / Doublet / Triplet);
* **5_get_SNR.py**
    - Usage of the function **get_SNR()** which returns an array of SNR's (corresponding for the output in each detector);
* **6_draw_detectors.py**
    - Usage of the function **draw_detectors()** which returns a world map figure indicating the positions, orientations and shapes of the detectors;
* **7_get_map.py**
    - Usage of the function **get_map()** which returns a world map with the projection of the detectors antenna pattern responses;
* **8_check_priors.py**
    - Usage of the function **Priors()** which returns a plot with the priors of the free parameters desired by the user;
* **9_gwdali_mcmc_Fisher-vs-Singlet.py**
    - Usage of the function **GWDALI()** for Fisher inversion (samples obtained via multivariate normal distribution) and with sampler method available in bilby package. It returns the samples obtained via MCMC/Nested_Sampling methods wich can be used to build corner plots. **Comparison between Fisher and Singlet**
* **10_gwdali_mcmc_Exact-vs-Doublet.py**
    - Usage of the function **GWDALI()** with sampler method available in bilby package. It returns the samples obtained via MCMC/Nested_Sampling methods wich can be used to build corner plots. **Comparison between Exact and Doublet**
* **11_gwdali_grid.py**
    - Usage of the function **GWDALI()** with sampler method **grid** which returns a grid the desired free parameters an its corresponding posterior values;
