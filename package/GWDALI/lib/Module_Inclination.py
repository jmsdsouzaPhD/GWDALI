# ===============================================================
# IOTA BUILDERS
# ===============================================================

from jax import jit

IOTA_BUILDERS = {}

def register_iota(name):
    def wrapper(fn):
        IOTA_BUILDERS[name] = fn
        return fn
    return wrapper
    
@register_iota("iota")
def build_iota(IDX):
    i = IDX["iota"]
    @jit
    def transform(theta):
        return theta[i]
    return transform

@register_iota("cos_iota")
def build_cos_iota(IDX):
    i = IDX["cos_iota"]
    @jit
    def transform(theta):
        return jnp.arccos(theta[i])
    return transform
