=================================
Derivative Approximation for LIkelihood (DALI)
=================================

.. math::

	log\mathcal{L} \approx log\mathcal{L}_{0}&
	 -\left[\frac{1}{2}\sum_{i,j}\left\langle \partial_{i}h(\theta)|\partial_{j}h(\theta)\right\rangle _{0}\Delta\theta_{i}\Delta\theta_{j}\right] \\
	& 
	-\left[\frac{1}{2}\sum_{i,j,k}\left\langle \partial_{i}h(\theta)|\partial_{j}\partial_{k}h(\theta)\right\rangle _{0}\Delta\theta^{ijk}
	+\frac{1}{8}\sum_{i,j,k,l}\left\langle \partial_{i}\partial_{j}h(\theta)|\partial_{k}\partial_{l}h(\theta)\right\rangle _{0}\Delta\theta^{ijkl}\right] \\
	& 
		-\left[\frac{1}{6}\sum_{i,\cdots, l}\left\langle \partial_{i}h|\partial_{j}\partial_{k}\partial_{l}h\right\rangle \Delta\theta^{ijkl}
		+\frac{1}{12}\sum_{i,\cdots, m}\left\langle \partial_{i}\partial_{j}h|\partial_{k}\partial_{l}\partial_{m}h\right\rangle \Delta\theta^{ijklm} \right.\\
	& \left.
		+\frac{1}{72}\sum_{i,\cdots, n}\left\langle \partial_{i}\partial_{j}\partial_{k}h|\partial_{l}\partial_{m}\partial_{n}h\right\rangle \Delta\theta^{ijklmn}
		\right] \\
	& +\mathcal{O}(\partial^{4})

where 
	..math::	

		`\Delta\theta^i =& \theta^i-\theta_0^i`  \\ 
		`\Delta\theta^{ij} =& \Delta\theta^i\cdot\Delta\theta^j`  \\ 
		`\Delta\theta^{ij..m} =& \Delta\theta^{i}\cdot\Delta\theta^j\cdots\Delta\theta^m`.