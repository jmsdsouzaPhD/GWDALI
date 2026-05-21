from jax import jit
import jax.numpy as jnp
import numpy as np
import matplotlib.pyplot as plt

@jit
def InnerProduct(A,B):
	A = jnp.array(A)
	B = jnp.array(B)
	return jnp.sum(A*B)

@jit
def CrossProduct(A,B):
	C1 = A[1]*B[2] - A[2]*B[1]
	C2 = A[2]*B[0] - A[0]*B[2]
	C3 = A[0]*B[1] - A[1]*B[0]
	return jnp.array( [C1, C2, C3] )

rad = jnp.pi/180

@jit # Analytic Formulation
def AngTransf(psi0,RA,Dec,lon,lat,rot):
	alpha0 = RA*rad ; beta0 = (90-Dec)*rad
	phi = lon*rad ; theta = (90-lat)*rad ; ksi = rot*rad

	ca0 = jnp.cos(alpha0-phi)
	sa0 = jnp.sin(alpha0-phi)
	cb0 = jnp.cos(beta0)
	sb0 = jnp.sin(beta0)
	ct = jnp.cos(theta)
	st = jnp.sin(theta)
	cx = jnp.cos(ksi)
	sx = jnp.sin(ksi)
	
	U = sb0*sa0
	V = sb0*ct*ca0 - cb0*st
	
	SinB_SinA = -V*sx + U*cx
	SinB_CosA =  V*cx + U*sx 
	CosB = sb0*st*ca0 + cb0*ct

	alpha_det = jnp.arctan2( SinB_SinA , SinB_CosA )
	beta_det = jnp.arccos(CosB)

	SinPsi = -st*sa0
	CosPsi = ct*sb0 - st*cb0*ca0

	psi_det = psi0 + jnp.arctan2(SinPsi,CosPsi)

	return alpha_det, beta_det, psi_det
