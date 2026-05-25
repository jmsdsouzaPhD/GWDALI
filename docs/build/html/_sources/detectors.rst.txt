=======================================
Detectors (Position/Orientation/Shape)
=======================================

To setting the properties of a given detector, create a python dictionary with the following entries:

.. py:function:: {"name","lon","lat","rot","shape"}

   :param str name: Name of detector among ("ET","CE","CE_20km","CE_40km","LIGO","Virgo","KAGRA");
   :param float lon: Longitude;
   :param float lat: Latitude;
   :param float rot: Orientation anti-clockwise starting in 0° from the detector bissetiz aligned to the Noth-South direction;
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
