# ===============================================================
# SPIN BUILDERS
# ===============================================================

from jax import jit

SPIN_BUILDERS = {}

def register_spin(name):
    def wrapper(fn):
        SPIN_BUILDERS[name] = fn
        return fn
    return wrapper

@register_spin("cartesian")
def build_cartesian_spins(IDX):
    isx1 = IDX["sx1"]
    isy1 = IDX["sy1"]
    isz1 = IDX["sz1"]

    isx2 = IDX["sx2"]
    isy2 = IDX["sy2"]
    isz2 = IDX["sz2"]

    @jit
    def transform(theta):
        return (
            theta[isx1],
            theta[isy1],
            theta[isz1],
            theta[isx2],
            theta[isy2],
            theta[isz2],
        )
    return transform

@register_spin("spherical")
def build_spherical_spins(IDX):
    iS1 = IDX["S1"]
    iT1 = IDX["theta1"]
    iP1 = IDX["phi1"]

    iS2 = IDX["S2"]
    iT2 = IDX["theta2"]
    iP2 = IDX["phi2"]

    @jit
    def transform(theta):
        S1, theta1, phi1 = theta[iS1], theta[iT1], theta[iP1]
        S2, theta2, phi2 = theta[iS2], theta[iT2], theta[iP2]
        sx1 = S1*jnp.sin(theta1)*jnp.cos(phi1)
        sy1 = S1*jnp.sin(theta1)*jnp.sin(phi1)
        sz1 = S1*jnp.cos(theta1)

        sx2 = S2*jnp.sin(theta2)*jnp.cos(phi2)
        sy2 = S2*jnp.sin(theta2)*jnp.sin(phi2)
        sz2 = S2*jnp.cos(theta2)
        return sx1,sy1,sz1,sx2,sy2,sz2
    return transform

@register_spin("chiS_chiA")
def build_chi_schi_a(IDX):
    ichi_s = IDX["chi_s"]
    ichi_a = IDX["chi_a"]
    @jit
    def transform(theta):
        chi_s = theta[ichi_s]
        chi_a = theta[ichi_a]
        sz1 = chi_s + chi_a
        sz2 = chi_s - chi_a
        return (0.0,0.0,sz1,0.0,0.0,sz2)
    return transform
