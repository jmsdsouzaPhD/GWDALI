# ===============================================================
# MASS BUILDERS
# ===============================================================

# Mc, M, m1, m2, ln_Mc 
# eta, q, ln_eta, inv_eta, deltaM

from jax import jit
import jax.numpy as jnp

func_Mc = {}
func_Mc["Mc"] = jit(lambda x: x)
func_Mc["ln_Mc"] = jit(lambda x: jnp.exp(x))
Mc_keys = ["Mc","ln_Mc"]

func_eta = {}
func_eta['eta'] = jit( lambda x: x )
func_eta['q'] = jit( lambda x: x/(1+x)**2 )
func_eta['inv_eta'] = jit( lambda x: 1./x )
func_eta['ln_eta'] = jit( lambda x: jnp.exp(x) )
func_eta['deltaM'] = jit( lambda x: .25*(1.-x*x) )

eta_keys = func_eta.keys()

MASS_BUILDERS = {}

def register_mass(name):
    def wrapper(fn):
        MASS_BUILDERS[name] = fn
        return fn
    return wrapper

#---------------//-----------------------

@register_mass("m1_m2")
def build_m1m2(IDX):
    im1 = IDX["m1"]
    im2 = IDX["m2"]
    @jit
    def transform(theta):
        m1 = theta[im1]
        m2 = theta[im2]
        M = m1 + m2
        eta = (m1 * m2) / M**2
        Mc = M * eta**(3./5)
        return Mc, eta
    return transform

@register_mass("m1_M")
@register_mass("M_m1")
def build_Mm2(IDX):
    im1 = IDX["m1"]
    iM = IDX["M"]
    @jit
    def transform(theta):
        m1 = theta[im1]
        M = theta[iM]
        m2 = M-m1
        eta = (m1 * m2) / M**2
        Mc = M * eta**(3./5.)
        return Mc, eta
    return transform

@register_mass("m2_M")
@register_mass("M_m2")
def build_m1M(IDX):
    im2 = IDX["m2"]
    iM = IDX["M"]
    @jit
    def transform(theta):
        m2 = theta[im2]
        M = theta[iM]
        m1 = M-m2
        eta = (m1 * m2) / M**2
        Mc = M * eta**(3./5.)
        return Mc, eta
    return transform

for Mc_key in Mc_keys:
    for eta_key in eta_keys:
        @register_mass(f"{Mc_key}_{eta_key}")
        @register_mass(f"{eta_key}_{Mc_key}")
        def build_Mc_eta(IDX, Mc_key=Mc_key, eta_key=eta_key):
            iMc = IDX[Mc_key]
            ieta = IDX[eta_key]
            func1 = func_Mc[Mc_key]
            func2 = func_eta[eta_key]
            @jit
            def transform(theta):
                return (
                    func1( theta[iMc] ),
                    func2( theta[ieta] ),
                )
            return transform

for Mc_key in Mc_keys:
    @register_mass(f"{Mc_key}_M")
    @register_mass(f"M_{Mc_key}")
    def build_Mc_M(IDX, Mc_key=Mc_key):
        iM = IDX["M"]
        iX = IDX[Mc_key]
        func = func_Mc[Mc_key]
        @jit
        def transform(theta):
            M = theta[iM]
            Mc = func( theta[iX] )
            eta = (Mc/M)**(5./3)
            return Mc, eta
        return transform
#---------------//-----------------------

for eta_key in eta_keys:
    @register_mass(f"M_{eta_key}")
    @register_mass(f"{eta_key}_M")
    def build_M_eta(IDX, eta_key=eta_key):
        iM = IDX["M"]
        iX = IDX[eta_key]
        func = func_eta[eta_key]
        @jit
        def transform(theta):
            M = theta[iM]
            X = theta[iX]
            eta = func(X)
            #eta = q / (1.0 + q)**2
            Mc = M * eta**(3.0 / 5.0)
            return Mc, eta
        return transform

for eta_key in eta_keys:
    @register_mass(f"m1_{eta_key}")
    @register_mass(f"{eta_key}_m1")
    def build_m1_eta(IDX, eta_key=eta_key):
        im1 = IDX["m1"]
        iX = IDX[eta_key]
        func = func_eta[eta_key]
        @jit
        def transform(theta):
            m1 = theta[im1]
            eta = func( theta[iX] )
            q = 1./(2*eta) * ( (1.-2*eta) - jnp.sqrt(1.-4*eta) )
            m2 = q * m1
            M = m1 + m2
            Mc = M * eta**(3./5)
            return Mc, eta
        return transform

for eta_key in eta_keys:
    @register_mass(f"m2_{eta_key}")
    @register_mass(f"{eta_key}_m2")
    def build_m2_eta(IDX, eta_key=eta_key):
        im2 = IDX["m2"]
        iX = IDX[eta_key]
        func = func_eta[eta_key]
        @jit
        def transform(theta):
            m2 = theta[im2]
            eta = func( theta[iX] )
            q = 1./(2*eta) * ( (1.-2*eta) - jnp.sqrt(1.-4*eta) )
            m1 = m2 / q
            M = m1 + m2
            Mc = M * eta**(3./5)
            return Mc, eta
        return transform

@jit
def _solve_ratio(b):
    disc = b*b*(0.25 - b/27) + 0j
    sqrt_disc = jnp.sqrt(disc)

    u3 = b/2 + sqrt_disc
    v3 = b/2 - sqrt_disc

    u = u3**(1/3)
    v = v3**(1/3)
    return jnp.real(u + v)

for Mc_key in Mc_keys:
    @register_mass(f"{Mc_key}_m2")
    @register_mass(f"m2_{Mc_key}")
    def build_Mc_m2(IDX,Mc_key=Mc_key):
        iM2 = IDX["m2"]
        iX = IDX[Mc_key]
        func = func_Mc[Mc_key]
        @jit
        def transform(theta):
            m2 = theta[iM2]
            Mc = func( theta[iX] )
            b = (Mc/m2)**5
            r = _solve_ratio(b)
            m1 = m2 * r
            eta = m1*m2/(m1+m2)**2
            return Mc, eta
        return transform

for Mc_key in Mc_keys:
    @register_mass(f"{Mc_key}_m1")
    @register_mass(f"m1_{Mc_key}")
    def build_Mc_m1(IDX,Mc_key=Mc_key):
        iM1 = IDX["m1"]
        iX = IDX[Mc_key]
        func = func_Mc[Mc_key]
        @jit
        def transform(theta):
            m1 = theta[iM1]
            Mc = func( theta[iX] )
            b = (Mc/m1)**5
            q = _solve_ratio(b)
            m2 = m1 * q
            eta = m1*m2/(m1+m2)**2
            return Mc, eta
        return transform