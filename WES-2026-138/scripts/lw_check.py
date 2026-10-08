"""Check of the Limacher & Wood (2021) relation  int_0^inf (a^2 - ur^2) r dr = 0  in the rotor plane
for the vortex-cylinder field (planar rotor) and for the Madsen expression."""
import numpy as np, json
from vcfull import ur_cyl, uz_cyl, strengths
from scipy.special import ellipk, ellipe
from driver import run
RT, RROOT = 63.0, 0.3
xg, wg = np.polynomial.legendre.leggauss(120)
tg = 0.5*(xg+1); wt = 0.5*wg

def integ(f, lo, hi):
    """integral of f over [lo,hi] with (integrable) log singularities at both ends: smoothstep substitution"""
    s = tg*tg*(3-2*tg); ds = 6*tg*(1-tg)
    r = lo + (hi-lo)*s
    return np.sum(wt*f(r)*ds)*(hi-lo)

def field_integrals(Rk, gam, U0, rmax_fac=2000):
    Rk = np.asarray(Rk, float); gam = np.asarray(gam, float)
    def ur1(r, R, g):
        k2 = np.minimum(4*r*R/(R+r)**2, 1-1e-14); k = np.sqrt(k2)
        return -g/(2*np.pi)*np.sqrt(R/r)*((2-k2)/k*ellipk(k2) - 2/k*ellipe(k2))
    ur = lambda r: sum(ur1(r, R, g) for R, g in zip(Rk, gam) if abs(g) > 1e-14 and R > 1e-6)
    # in the plane z = 0 the axial induction is piecewise constant: a = -sum_{Rk>r} gam_k/(2 U0)
    edges = np.concatenate([[1e-6], Rk]) if Rk[0] > 1e-6 else Rk.copy()
    Ia = Ir_in = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        aj = -np.sum(gam[Rk > 0.5*(lo+hi)])/(2*U0)
        Ia += aj**2*(hi*hi-lo*lo)/2; Ir_in += integ(lambda r: (ur(r)/U0)**2*r, lo, hi)
    far = np.concatenate([[Rk[-1]], np.geomspace(Rk[-1]*1.05, Rk[-1]*rmax_fac, 60)])
    Ir_out = sum(integ(lambda r: (ur(r)/U0)**2*r, lo, hi) for lo, hi in zip(far[:-1], far[1:]))
    return Ia, Ir_in, Ir_out

out = {}
for a0 in (0.1, 0.3):
    Ia, Ii, Io = field_integrals([1.0], [-2*a0], 1.0)
    out[f'uniform_a{a0}'] = dict(Ia=Ia, Ir_in=Ii, Ir_out=Io, ratio=(Ii+Io)/Ia)
    print(f'uniform disc a={a0}: int a^2 r = {Ia:.6f}, int ur^2 r = {Ii+Io:.6f} (inside {Ii:.6f}, outside {Io:.6f}), ratio = {(Ii+Io)/Ia:.4f}')
for vw, rpm, pitch in ((11.0, 12.1, 0.0), (7.0, 8.0, 0.0), (15.0, 12.1, 10.45)):
  for n in ((35, 70, 140) if vw == 11.0 else (70,)):
    res = run(vw, rpm, pitch, nsec=n, radial=True, radmode=2, npform=1)
    dr = (RT-RROOT)/n; Rk = RROOT + dr*np.arange(n+1)
    for lab, aa in (('ainf', res['vc']['ainf']), ('aB', res['vc']['a'])):
        gam = strengths(aa, vw)
        Ia, Ii, Io = field_integrals(Rk, gam, vw)
        out[f'nrel_{vw}_{n}_{lab}'] = dict(Ia=Ia, Ir_in=Ii, Ir_out=Io, ratio=(Ii+Io)/Ia)
        print(f'NREL planar {vw} m/s Nsec={n} {lab}: int a^2 r = {Ia:.3f}, int ur^2 r = {Ii+Io:.3f} (in {Ii:.3f}, out {Io:.3f}), ratio = {(Ii+Io)/Ia:.4f}')
    if n == 70:
        rM = run(vw, rpm, pitch, nsec=n, radial=True, radmode=1, npform=1)
        r = rM['bem']['r']; urM = rM['bem']['ur']; a = res['vc']['ainf']
        IM = np.sum((urM/vw)**2*r)*dr; Iain = np.sum(a**2*r)*dr
        # Madsen expression continued outside the disc with the rotor-averaged cT (cT_av ~ CT R^2/r^2)
        CT = rM['cT']; rr = np.geomspace(RT*1.0005, RT*2000, 4000); x = rr/RT
        uo = vw/2.24*(CT/x**2)/(4*np.pi)*np.log((0.04**2+(x+1)**2)/(0.04**2+(x-1)**2))
        IMo = np.trapezoid((uo/vw)**2*rr, rr)
        out[f'madsen_{vw}'] = dict(Ia=Iain, Ir_in=IM, Ir_out=IMo, ratio_in=IM/Iain, ratio=(IM+IMo)/Iain)
        print(f'   Madsen {vw} m/s: int a^2 r = {Iain:.3f}, int ur^2 r inside disc = {IM:.3f} (ratio {IM/Iain:.4f}); with continuation outside {IM+IMo:.3f} (ratio {(IM+IMo)/Iain:.4f})')
json.dump(out, open('lw_check.json','w'), indent=1)
