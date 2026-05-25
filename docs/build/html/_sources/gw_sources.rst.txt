GW Source Characterization
==========================

In GWDALI, a gravitational-wave source is described through a set of
free parameters defining the binary intrinsic and extrinsic properties.

The framework supports multiple parameterizations for masses,
distance, inclination, and spins. Independently of the chosen
parameterization, a valid source description must contain:

- exactly 2 mass parameters;
- exactly 1 distance parameter;
- exactly 1 inclination parameter;
- exactly 6 spin parameters.

Additionally, detector-frame waveform generation requires:

- ``RA``
- ``Dec``
- ``psi``
- ``t_coal``
- ``phi_coal``


Mass Parameterizations
----------------------

The following mass parameterizations are currently supported:

.. list-table::
   :widths: 20 60
   :header-rows: 1

   * - Parameter
     - Description

   * - ``m1``
     - Redshifted mass of the primary compact object,
       :math:`(1+z)m_1` [:math:`M_\odot`].

   * - ``m2``
     - Redshifted mass of the secondary compact object,
       :math:`(1+z)m_2` [:math:`M_\odot`].

   * - ``eta``
     - Symmetric mass ratio,

       .. math::

          \eta \equiv \frac{m_1m_2}{(m_1+m_2)^2}

   * - ``Mc``
     - Redshifted chirp mass,

       .. math::

          \mathcal{M}_c \equiv
          (1+z)\eta^{3/5}(m_1+m_2)

       [:math:`M_\odot`].

   * - ``M``
     - Redshifted total mass,

       .. math::

          M \equiv (1+z)(m_1+m_2)

       [:math:`M_\odot`].

   * - ``q``
     - Mass ratio,

       .. math::

          q \equiv \frac{m_2}{m_1},
          \qquad m_1 > m_2

   * - ``inv_eta``
     - Inverse symmetric mass ratio,
       :math:`\eta^{-1}`.

   * - ``ln_Mc``
     - Logarithm of the redshifted chirp mass,

       .. math::

          \ln\left(\frac{(1+z)M_c}{M_\odot}\right)

   * - ``ln_eta``
     - Logarithm of the symmetric mass ratio,
       :math:`\ln(\eta)`.

   * - ``deltaM``
     - Dimensionless mass difference,

       .. math::

          \delta_M \equiv \frac{m_1-m_2}{m_1+m_2}


Distance Parameterizations
--------------------------

Exactly one distance parameter must be provided.

.. list-table::
   :widths: 20 60
   :header-rows: 1

   * - Parameter
     - Description

   * - ``dL``
     - Luminosity distance :math:`d_L` [Gpc].

   * - ``inv_dL``
     - Inverse luminosity distance,
       :math:`d_L^{-1}` [:math:`\mathrm{Gpc}^{-1}`].

   * - ``inv_dL2``
     - Inverse squared luminosity distance,
       :math:`d_L^{-2}` [:math:`\mathrm{Gpc}^{-2}`].

   * - ``inv_sqrtdL``
     - Inverse square-root luminosity distance,
       :math:`d_L^{-1/2}` [:math:`\mathrm{Gpc}^{-1/2}`].

   * - ``lnDL``
     - Logarithm of the luminosity distance,

       .. math::

          \ln(d_L/\mathrm{Gpc})

   * - ``inv_lnDL``
     - Inverse logarithmic luminosity distance,

       .. math::

          \frac{1}{\ln(d_L/\mathrm{Gpc})}


Inclination Parameterizations
-----------------------------

Exactly one inclination parameter must be provided.

.. list-table::
   :widths: 20 60
   :header-rows: 1

   * - Parameter
     - Description

   * - ``iota``
     - Inclination angle :math:`\iota` [rad].

   * - ``cos_iota``
     - Cosine of the inclination angle,
       :math:`\cos(\iota)`.


Spin Parameterizations
----------------------

Exactly six spin parameters must be provided.

Cartesian Spin Coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :widths: 20 60
   :header-rows: 1

   * - Parameter
     - Description

   * - ``sx1``, ``sy1``, ``sz1``
     - Cartesian spin components of the primary compact object.

   * - ``sx2``, ``sy2``, ``sz2``
     - Cartesian spin components of the secondary compact object.


Spherical Spin Coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :widths: 20 60
   :header-rows: 1

   * - Parameter
     - Description

   * - ``chi1``
     - Dimensionless spin magnitude of the primary compact object,
       :math:`\chi_1 \equiv |\vec{S}_1|`.

   * - ``theta1``
     - Polar angle of the primary spin [rad].

   * - ``phi1``
     - Azimuthal angle of the primary spin [rad].

   * - ``chi2``
     - Dimensionless spin magnitude of the secondary compact object,
       :math:`\chi_2 \equiv |\vec{S}_2|`.

   * - ``theta2``
     - Polar angle of the secondary spin [rad].

   * - ``phi2``
     - Azimuthal angle of the secondary spin [rad].


Symmetric/Antisymmetric Spin Coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :widths: 20 60
   :header-rows: 1

   * - Parameter
     - Description

   * - ``chi_s``
     - Symmetric spin combination,

       .. math::

          \chi_s = \frac{\chi_1+\chi_2}{2}

   * - ``chi_a``
     - Antisymmetric spin combination,

       .. math::

          \chi_a = \frac{\chi_1-\chi_2}{2}


Extrinsic Parameters
--------------------

.. list-table::
   :widths: 20 60
   :header-rows: 1

   * - Parameter
     - Description

   * - ``RA``
     - Right ascension [deg].

   * - ``Dec``
     - Declination [deg].

   * - ``psi``
     - Polarization angle [rad].

   * - ``phi_coal``
     - Coalescence phase [rad].

   * - ``t_coal``
     - Coalescence time [s].

Example: Source Dictionary Construction
---------------------------------------

Below we show a complete example of a gravitational-wave source
dictionary compatible with ``get_strain()``.

This example uses:

- ``(Mc, eta)`` as mass parameterization;
- ``dL`` as distance parameterization;
- ``iota`` as inclination parameterization;
- Cartesian spin coordinates.

The resulting dictionary contains the required 15 parameters.

.. code-block:: python

   import numpy as np

   GwPrms = {}

   #========================================
   # Distance sector (1 parameter)
   #========================================

   GwPrms['dL'] = 5.0          # Luminosity distance [Gpc]

   #========================================
   # Mass sector (2 parameters)
   #========================================

   GwPrms['Mc']  = 30.0        # Chirp mass [Msun]
   GwPrms['eta'] = 0.24        # Symmetric mass ratio

   #========================================
   # Inclination sector (1 parameter)
   #========================================

   GwPrms['iota'] = 1.2        # Inclination [rad]

   #========================================
   # Extrinsic parameters
   #========================================

   GwPrms['RA']         = 120.0     # Right ascension [deg]
   GwPrms['Dec']        = -30.0     # Declination [deg]
   GwPrms['psi']        = 0.5       # Polarization angle [rad]
   GwPrms['t_coal']     = 0.0       # Coalescence time [s]
   GwPrms['phi_coal']   = 0.0       # Coalescence phase [rad]

   #========================================
   # Spin sector (6 parameters)
   #========================================

   GwPrms['sx1'] = 0.0
   GwPrms['sy1'] = 0.0
   GwPrms['sz1'] = 0.5

   GwPrms['sx2'] = 0.0
   GwPrms['sy2'] = 0.0
   GwPrms['sz2'] = -0.3

   print(GwPrms)


The source dictionary can then be passed directly to the waveform
generation routines:

.. code-block:: python

   import GWDALI as gw

   freq = np.arange(1.0, 1024.0, 0.1)

   h = gw.get_strain(
       detectors,
       GwPrms,
       freq,
       approx="IMRPhenomD"
   )
   
Sky Localization
----------------

In GWDALI, astronomical coordinates (RA, Dec) are aligned with the
geocentric coordinates (longitude ``lon`` and latitude ``lat``),
such that:

.. math::

   \mathrm{RA} \parallel \mathrm{lon},
   \qquad
   \mathrm{Dec} \parallel \mathrm{lat}

.. figure:: ./_static/geo_coords.png
   :alt: Source Coordinates
   :align: center
   :scale: 50%