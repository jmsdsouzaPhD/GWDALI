"""
===============================================================
Differentiable GW parameter transformations in JAX
===============================================================

FINAL ARCHITECTURE

Goals:
- fully differentiable
- JAX-native
- no Python branching inside jit
- no dicts inside jit
- no dynamic tracing
- modular parameterizations
- autodiff-ready
- Fisher-ready

Pipeline:

generalized parameters
        ↓
differentiable coordinate transform
        ↓
canonical parameters
        ↓
waveform
        ↓
autodiff

"""

from typing import NamedTuple

import jax
import jax.numpy as jnp

from jax import jit
from jax import jacrev

jax.config.update("jax_enable_x64", True)
# ===============================================================
# CANONICAL PARAMETERS
# ===============================================================

class CanonicalParams(NamedTuple):

    # extrinsic
    dL: float
    iota: float
    psi: float
    phi_coal: float
    RA: float
    Dec: float
    t_coal: float

    # masses
    Mc: float
    eta: float

    # spins
    sx1: float
    sy1: float
    sz1: float

    sx2: float
    sy2: float
    sz2: float


# ===============================================================
# REGISTRIES
# ===============================================================
# @register_distance
# @register_iota
# @register_mass
# @register_spin

from .Module_Distance import *
from .Module_Inclination import *
from .Module_Mass import *
from .Module_Spins import *

# ===============================================================
# PARAMETRIZATION DETECTION
# ===============================================================

def detect_distance_kind(keys):
    for k in "dL,inv_dL,inv_dL2,inv_sqrtdL,inv_lnDL,lnDL".split(','): 
        if k in keys: return k
    raise ValueError("No distance parametrization found")

def detect_iota_kind(keys):
    if "iota" in keys:      return "iota"
    if "cos_iota" in keys:  return "cos_iota"
    raise ValueError("No iota parametrization found")

def detect_mass_kind(keys):
    Mass_keys = "Mc,eta,m1,m2,M,q,deltaM,ln_Mc,ln_eta,inv_eta".split(',')
    n = len(Mass_keys)
    for i in range(n):
        key1 = Mass_keys[i]
        for j in range(i):
            key2 = Mass_keys[j]
            if key1 in keys and key2 in keys: return f"{key2}_{key1}"
    raise ValueError("Unsupported mass parametrization")

def detect_spin_kind(keys):
    if "chi_s" in keys and "chi_a" in keys:  return "chiS_chiA"
    if all([xi in keys for xi in "sx1,sy1,sz1,sx2,sy2,sz2".split(',')]): return "cartesian"
    if all([xi in keys for xi in "S1,theta1,phi1,S1,theta2,phi2".split(',')]): return "spherical"
    raise ValueError("Unsupported spins parametrization")

# ===============================================================
# BUILD TRANSFORM
# ===============================================================

def build_transform(theta_keys):
    theta_keys = tuple(theta_keys)
    IDX = {
        k: i for i, k in enumerate(theta_keys)
    }
    # -----------------------------------------------------------
    # detect parametrizations OUTSIDE JIT
    # -----------------------------------------------------------

    distance_kind = detect_distance_kind(theta_keys)
    iota_kind = detect_iota_kind(theta_keys)
    mass_kind = detect_mass_kind(theta_keys)
    spin_kind = detect_spin_kind(theta_keys)
    #print(">> distance_kind:",distance_kind) ; quit()
    # -----------------------------------------------------------
    # build transforms
    # -----------------------------------------------------------

    # Found in Module_Distance/Inclination/Mass/Spins.py
    distance_transform = DISTANCE_BUILDERS[distance_kind](IDX)
    iota_transform = IOTA_BUILDERS[iota_kind](IDX)
    mass_transform = MASS_BUILDERS[mass_kind](IDX)
    spin_transform = SPIN_BUILDERS[spin_kind](IDX)

    # -----------------------------------------------------------
    # optional params
    # -----------------------------------------------------------

    ipsi = IDX.get("psi", None)
    iphi = IDX.get("phi_coal", None)

    iRA = IDX.get("RA", None)
    iDec = IDX.get("Dec", None)

    itcoal = IDX.get("t_coal", None)

    # -----------------------------------------------------------
    # differentiable transform
    # -----------------------------------------------------------

    @jit
    def transform(theta):
        dL = distance_transform(theta)
        iota = iota_transform(theta)
        Mc, eta = mass_transform(theta)
        (sx1,sy1,sz1,sx2,sy2,sz2) = spin_transform(theta)

        psi = theta[ipsi]
        phi_coal = theta[iphi]
        RA = theta[iRA]
        Dec = theta[iDec]
        t_coal = theta[itcoal]
 
        return [dL,iota,psi,phi_coal,RA,Dec,t_coal,Mc,eta,sx1,sy1,sz1,sx2,sy2,sz2]

    return transform
