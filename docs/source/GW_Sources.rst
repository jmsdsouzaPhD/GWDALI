=================================  
GW Sources Characterization
=================================

************************************
   Sky Localization
************************************

In GWDALI we deal with astronomical coordinates (RA, Dec) aligned with geocentric coordinates (longitude 'lon' and latitude 'lat' ) such that :math:`Ra//lon` and :math:`Dec//lat`:

.. figure:: ./geo_coords.png
   :alt: Source Coordinates
   :align: center
   :scale: 50%

************************************
   Free Parameters
************************************

   * m1: **Redshifted  Mass of the first object** :math:`\longrightarrow (1+z)m_1` [:math:`M_{\odot}`]
   * m2: **Redshifted  Mass of the second object** :math:`\longrightarrow (1+z)m_2` [:math:`M_{\odot}`]
   * eta: **Symetric mass ratio** :math:`\longrightarrow \eta \equiv m_1m_2/(m_1+m_2)^2`
   * Mc: **Redshifted Chirp Mass** :math:`\longrightarrow (1+z)M_c \equiv (1+z)\eta^{3/5}(m_1+m_2)` [:math:`M_{\odot}`]
   * q: **Mass Ratio** :math:`\longrightarrow q=m_1/m_2` with :math:`m_2>m_1`
   * DL: **Luminosity Distance** :math:`\longrightarrow d_L` [Gpc]
   * inv_dL: **Inverse of Luminosity Distance** :math:`\longrightarrow d_L^{-1}` [:math:`Gpc^{-1}`]
   * ln_dL: **Logarithm of Luminosity Distance** :math:`\longrightarrow ln(d_L/Gpc)`
   * RA: **Right Ascencion** :math:`\longrightarrow Ra` [deg]
   * Dec: **Declination** :math:`\longrightarrow Dec` [deg]
   * iota: **Binary Inclination** :math:`\longrightarrow \iota` [rad]
   *cos\_iota: *Cosine of Inclination* :math:`\longrightarrow cos(\iota)`
   * psi: *Polarization Angle* :math:`\longrightarrow \psi` [rad]
   * phi\_coal: *Coalescence Phase* :math:`\longrightarrow \phi_{coal}` [rad]
   * t\_coal: *Coalescence Time* :math:`\longrightarrow t_{coal}` [seconds]
   * sx1: *X-Component of First Object Spin* :math:`\longrightarrow S_{x_1}`
   * sy1: *Y-Component of First Object Spin* :math:`\longrightarrow S_{y_1}`
   * sz1: *Z-Component of First Object Spin* :math:`\longrightarrow S_{z_1}`
   * sx2: *X-Component of Second Object Spin* :math:`\longrightarrow S_{x_2}`
   * sy2: *Y-Component of Second Object Spin* :math:`\longrightarrow S_{y_2}`
   * sz2: *Z-Component of Second Object Spin* :math:`\longrightarrow S_{z_2}`
   * S1: *Magnitude of the First Object Spin* :math:`\longrightarrow |\vec{S}_1|`
   * theta_1: *Polar Angle of the First Object Spin* :math:`\longrightarrow \theta_1` [rad]
   * phi_1: *Azimuthal Angle of the First Object Spin* :math:`\longrightarrow phi_1` [rad]
   * S2 *Magnitude of the Second Object Spin* :math:`\longrightarrow |\vec{S}_2|`
   * theta_2: *Polar Angle of the Second Object Spin* :math:`\longrightarrow \theta_2` [rad]
   * phi_2: *Azimuthal Angle of the Second Object Spin* :math:`\longrightarrow phi_2` [rad]


