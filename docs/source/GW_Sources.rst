=================================  
GW Sources Characterization
=================================
In GWDALI we deal with astronomical coordinates (RA, Dec) aligned with geocentric coordinates (longitude 'lon' and latitude 'lat' ) such that :math:`Ra//lon` and :math:`Dec//lat`:

.. figure:: ./geo_coords.png
   :alt: Source Coordinates
   :align: center
   :scale: 50%

************************************
   Free Parameters
************************************

   * m1: Redshifted  Mass of the first object :math:`\longrightarrow (1+z)m_1` [:math:`M_{\odot}`]
   * m2: Redshifted  Mass of the second object :math:`\longrightarrow (1+z)m_2` [:math:`M_{\odot}`]
   * eta: Symetric mass ratio :math:`\longrightarrow \eta \equiv m_1m_2/(m_1+m_2)^2`
   * Mc: Redshifted Chirp Mass:math:`\longrightarrow (1+z)M_c \quiv (1+z)\eta^{3/5}(m_1+m_2)` [:math:`M_{\odot}`]
   * q: Mass Ratio :math:`\longrightarrow q=m_1/m_2` with :math:`m_2>m_1`
   * DL: Luminosity Distance :math:`\longrightarrow d_L` [Gpc]
   * inv_dL: Inverse of Luminosity Distance :math:`\longrightarrow d_L^{-1}` [:math:`Gpc^{-1}`]
   * ln_dL: Logarithm of Luminosity Distance :math:`\longrightarrow ln(d_L/Gpc)`


