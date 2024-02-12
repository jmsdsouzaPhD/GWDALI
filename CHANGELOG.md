# All notable changes to this project will be documented in this file.

## [0.1.1] - 2024-02-07
### Added
- DALI-tensors $\langle\partial_{i...}h|\partial_{j...}h\rangle$ allowed as output (with the key *Tensors*, e.g. **Fisher, Doublet, Triplet = res['Tensors']**).
- Allow users to choose between running Monte-Carlo or not (e.g. the user would be interested only in DALI-tensors components insted estimating posteriors)

### Fixed
- Parameter Displacement ($\Delta\theta^i$) written incorrectely in **GWDALI/lib/Likelihood.py**

---
## <span style="color:red">WARNING!</span>

- Error found in **GWDALI/lib/Likelihood.py** (in all versions < 0.1.1):
    - line 72: Parameter Displacement ($\Delta\theta^i$) written incorrectely: <span style="color:red">*dT = self.Theta0 - Theta*</span>
    - In the next version we are going to fix it to: <span style="color:green">*dT = Theta - self.Theta0*</span>
---

## [0.1.0] - 2024-01-19
### Added
- Priors are loaded correctlly now.

### Fixed
- Priors are not loaded correctelly!
- The **Standard Likelihood** ( $log\mathcal{L}=-\frac{1}{2}\langle h_{data}-h_{model}|h_{data}-h_{model}\rangle$ ) was not implemented correctelly to handle *Dec* in units of degree.

## [0.0.9] - 2023-09-23

### Changed
- Changing prior on Declination "Dec" from **bilby.core.prior.Cosine(...)** (with *Dec* in radian units) to **bilby.core.prior.Interped(...)** (with *Dec* in degree units)

### Bugs
- Priors are not loaded correctelly!
- The Standard Likelihood was not implemented correctelly to handle *Dec* in units of degree.

## [0.0.8] - 2023-09-22
### Added
- **Arrival Time Delay** ($\tau_a$) included in signal computations. $h_a\rightarrow h_a\cdot e^{2\pi i f \tau_a}$. 
  
### Changed
- speed of light changed from $\boxed{c=3\cdot10^8\ m/s}$ to $\boxed{c=299792458\ m/s}$
- optimizing computation of the observed angles: **sky position** $(\alpha_D,\ \beta_D)$ and **polarization angle** ($\psi_D$).

## [0.0.7] - 2023-07-17
### Added
- Allow user to choose **any** waveform from *lalsimulation*!
  
### Changed
- Excluding waveforms call from **GWDALI/lib/Load_Dictionaries.py**
- Waveforms are now directelly called from **GWDALI/lib/Waveforms.py**!

## [0.0.6] - 2023-06-21
### Changed
- Excluding block command about computing DALI tensor components (brute force method) in **GWDALI/lib/Compute_Tensors.py**

## [0.0.5] - 2023-06-12
### Changed
- setup.py to include .txt files about detector sensitivities.

### Fixed
- detectors sensitivities absent (**GWDALI/Detectors_Sensitivity/**) [.txt files]

## [0.0.4] - 2023-06-12
### Added
- Allow users to change the standard priors.
- Allow users to check/compute sources SNR before Monte-Carlo simulations.
   
### Changed
- GWFISH comment (about matrix inversion) excluded from **GWDALI.py**

### Bugs
- detectors sensitivities absent (**GWDALI/Detectors_Sensitivity/**) [.txt files]

## [0.0.3] - 2023-06-09
### Fixed
- libraries absent (**GWDALI/lib/**)

### Bugs
- detectors sensitivities absent (**GWDALI/Detectors_Sensitivity/**)

## [0.0.2] - 2023-05-30
Initial Release

### Bugs
- libraries absent (**GWDALI/lib/**)
- detectors sensitivities absent (**GWDALI/Detectors_Sensitivity/**)

