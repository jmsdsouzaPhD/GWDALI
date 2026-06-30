=======================================
Detectors (Position/Orientation/Shape)
=======================================

To setting the properties of a given detector, create a python dictionary with the following entries:

.. py:function:: {"name","lon","lat","rot","shape"}

   :param str name: Name of detector among ("ET","CE","CE_20km","CE_40km","LIGO","Virgo","KAGRA");
   :param float lon: Longitude;
   :param float lat: Latitude;
   :param float rot: Counter-clockwise orientation angle (in degrees). For the standard L-shaped interferometer (shape=90°), rot=0° corresponds to one arm pointing South and the other East, so that the detector bisector points South-East. For arbitrary arm opening angles (shape ≠ 90°), the bisector orientation remains fixed while the two arms rotate symmetrically by ±shape/2 around it.
   :param float shape: Openning angle of the interferometer's arms.

**Example**:

.. code:: python

	detector = {"name":"ET","lon":5,"lat":10,"rot":45,"shape":60}

.. figure:: ./_static/Coords.png
   :alt: detectors
   :align: center
   :scale: 70%

************************************  
Detectors Sensitivity
************************************

.. figure:: ./_static/Sensitivity.png
   :alt: sensitivity
   :align: center
   :scale: 80%
