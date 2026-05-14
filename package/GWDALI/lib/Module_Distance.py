# ===============================================================
# DISTANCE BUILDERS
# ===============================================================

import jax.numpy as jnp
from jax import jit

DISTANCE_BUILDERS = {}

def register_distance(name):
    def wrapper(fn):
        DISTANCE_BUILDERS[name] = fn
        return fn    
    return wrapper

func_dL = {}
func_dL["dL"]         = jit( lambda x: x )
func_dL["inv_dL"]     = jit( lambda y: 1./y )
func_dL["inv_dL2"]    = jit( lambda y: 1./jnp.sqrt(y) )
func_dL["inv_sqrtdL"] = jit( lambda y: 1./y**2 )
func_dL["inv_lnDL"]   = jit( lambda y: jnp.exp(1./y) )
func_dL["lnDL"]      = jit( lambda y: jnp.exp(y) )

keys_dL = func_dL.keys()

for key in keys_dL:
    @register_distance(key)
    def build_distance_dL(IDX, key=key):
        i = IDX[key]
        func = func_dL[key]
        @jit
        def transform(theta):
            return func( theta[i] )
        return transform

'''
@register_distance("ln_dL")
def build_distance_ln_dL(IDX):
    i = IDX["ln_dL"]
    @jit
    def transform(theta):
        return jnp.exp(theta[i])
    return transform

@register_distance("inv_dL")
def build_distance_inv_dL(IDX):
    i = IDX["inv_dL"]
    @jit
    def transform(theta):
        return 1.0 / theta[i]
    return transform

@register_distance("inv_dL2")
def build_distance_inv_dL2(IDX):
    i = IDX["inv_dL2"]
    @jit
    def transform(theta):
        return 1.0 / jnp.sqrt(theta[i])
    return transform

@register_distance("inv_sqrtdL")
def build_distance_inv_sqrtdL(IDX):
    i = IDX["inv_sqrtdL"]
    @jit
    def transform(theta):
        return 1.0 / theta[i]**2
    return transform

@register_distance("inv_lndL")
def build_distance_inv_lndL(IDX):
    i = IDX["inv_lndL"]
    @jit
    def transform(theta):
        return jnp.exp(1.0 / theta[i])
    return transform
'''