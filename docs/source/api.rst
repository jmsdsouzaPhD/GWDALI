API Reference
=============

This page documents the main public functions available in GWDALI v1.0.

Main Interface
===============

.. py:function:: GWDALI.GWDALI(GwPrms, detectors, FreeParams, approx="TaylorF2", method="Doublet", sampler="nestle", diff_method="autodiff", ...)

   Main interface of GWDALI.

   Computes:
   
   - detector strain
   - signal-to-noise ratio (SNR)
   - Fisher matrix
   - DALI tensors
   - posterior samples

   Parameters
   ----------

   GwPrms : dict
      Dictionary containing the gravitational-wave source parameters.

   detectors : list
      List containing detector configurations.

   FreeParams : list
      Parameters varied during the inference procedure.

   approx : str, optional
      Waveform approximant.

   method : str, optional
      Tensor expansion method:
      
      - ``"Fisher"``
      - ``"Doublet"``
      - ``"Triplet"``

   sampler : str, optional
      Posterior sampler backend.

   diff_method : str, optional
      Derivative computation method:
      
      - ``"autodiff"``
      - ``"numdiff"``

   Returns
   -------

   Results : object
      Posterior samples and sampler output.

   Truths : list
      Injection values for the free parameters.

   Tensors : dict
      DALI tensors.

   Fisher : ndarray
      Fisher information matrix.

   runtimes : list
      Timing diagnostics.


Waveform Functions
==================

.. py:function:: GWDALI.get_hphx(detectors, GwPrms, freq, approx, enable_jax_waveforms=True)

   Computes the plus and cross GW polarizations.

.. py:function:: GWDALI.get_strain(detectors, GwPrms, freq, approx, enable_jax_waveforms=True)

   Computes the detector strain for a detector network.

.. py:function:: GWDALI.get_SNR(detectors, GwPrms, approx, enable_jax_waveforms=True)

   Computes individual and network signal-to-noise ratios.


Derivatives
============

.. py:function:: GWDALI.get_derivatives(FreeParams, approx, GwPrms, dets, freq, diff_order="first", diff_method="numdiff", ...)

   Computes waveform derivatives.

   Parameters
   ----------

   diff_order : str
      Derivative order:
      
      - ``"first"``
      - ``"second"``
      - ``"third"``

   diff_method : str
      Derivative method:
      
      - ``"autodiff"``
      - ``"numdiff"``

   Returns
   -------

   Diff_values : ndarray
      Derivative tensors.

   time_diff : list
      Runtime information.


DALI Tensors
=============

.. py:function:: GWDALI.get_dali_tensors(GwPrms, detectors, FreeParams, method, approx, ...)

   Computes Fisher, Doublet and Triplet tensors.

   Returns
   -------

   Tensors : dict

      Dictionary containing:

      - ``Fisher``
      - ``Doublet``
      - ``Triplet``


Priors
=======

.. py:function:: GWDALI.Priors(FreeParams, name=None, new_priors=None, plot=False)

   Builds prior distributions for Bayesian inference.


Detector Utilities
===================

.. py:function:: GWDALI.get_map(detectors, plot_map=False)

   Computes detector antenna response maps.

.. py:function:: GWDALI.draw_detectors(detectors)

   Draws detector locations on the Earth map.


Waveform Backends
==================

GWDALI v1.0 supports:

- JAX waveform models
- LALSuite waveform models

Supported approximants include:

- ``TaylorF2``
- ``TaylorF2_ISCO``
- ``TaylorF2_Spinless``
- ``IMRPhenomA``
- ``IMRPhenomB``
- ``IMRPhenomC``
- ``IMRPhenomD``
- ``IMRPhenomHM``


Automatic Differentiation
==========================

GWDALI provides automatic differentiation using JAX:

- ``jax.grad()``
- ``jax.vmap()``

for efficient computation of:

- waveform derivatives
- Fisher matrices
- higher-order DALI tensors


Example
========

.. code-block:: python

   import GWDALI as gw

   Results, Truths, Tensors, Fisher, runtimes = gw.GWDALI(
       GwPrms,
       detectors,
       FreeParams,
       approx="IMRPhenomHM",
       method="Doublet",
       diff_method="autodiff"
   )