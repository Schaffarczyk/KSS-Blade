#!/usr/bin/env python3
# Nsec convergence study of the vortex-cylinder u_r sweep (KSS V3, Sub3.f)
# Semi-analytical replica: fixed smooth a(r) (cubic spline of the converged
# Nsec=70 BEM solution, Bem_VC.out), scipy elliptic integrals (double prec.).
# Validated against Fortran VC.out at Nsec=70: max|delta u_r| = 2.3e-2 m/s.
# APS / Claude, 2026-07-31
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import ellipk, ellipe

U0, Rhub, Rtip = 11.0, 0.3, 63.0

bla = [l.split() for l in open('BlaDes.in').readlines()[1:] if l.strip()]
r_des = np.array([float(p[0]) for p in bla]); oop = np.array([float(p[4]) for p in bla])
oopspl = CubicSpline(r_des, oop, bc_type='natural')
rows = []
for line in open('Bem_VC.out').readlines()[1:]:
    p = line.split()
    if len(p) >= 27: rows.append((float(p[0]), float(p[1]), float(p[16])))
r70, a70, Ga70 = np.array(rows).T
aspl = CubicSpline(r70, a70, bc_type='natural'); Gaspl = CubicSpline(r70, Ga70, bc_type='natural')

def vc(nsec, clip):
    dr = (Rtip-Rhub)/nsec
    rm = Rhub + 0.5*dr + np.arange(nsec)*dr
    zm = oopspl(rm)
    Rk = Rhub + np.arange(nsec+1)*dr
    zk = oopspl(np.clip(Rk, r_des[0], r_des[-1]))
    am = aspl(np.clip(rm, r70[0], r70[-1]))
    ae = np.concatenate([[0.], am, [0.]])
    gamt = 2.*U0*(ae[1:]-ae[:-1])
    ur = np.zeros(nsec)
    for k in range(nsec+1):
        if abs(gamt[k]) < 1e-10 or Rk[k] < 1e-3: continue
        x = zm - zk[k]
        k2 = 4.*rm*Rk[k]/((Rk[k]+rm)**2 + x*x)
        k2 = np.clip(k2, 1e-14, 0.999999 if clip else 1-1e-14)
        sk = np.sqrt(k2)
        br = (2.-k2)/sk*ellipk(k2) - 2./sk*ellipe(k2)
        ur += -gamt[k]/(2*np.pi)*np.sqrt(Rk[k]/rm)*br
    return rm, ur, dr

Ns = [35, 70, 140, 280, 560, 1120, 2240]
print(f"{'N':>6} {'ur_last(clip)':>14} {'ur_last(free)':>14} {'ur@62.55':>9} {'Moment/kNm':>11}")
ylast = []
for n in Ns:
    _, urc, _ = vc(n, True)
    rmf, urf, dr = vc(n, False)
    Ga = Gaspl(np.clip(rmf, r70[0], r70[-1])); kap = np.arctan(oopspl(rmf, 1))
    M = 3*np.sum(1.225*Ga*urf*np.tan(kap)*dr*rmf)/1000
    ylast.append(urf[-1])
    print(f"{n:6d} {urc[-1]:14.3f} {urf[-1]:14.3f} {np.interp(62.55,rmf,urf):9.3f} {M:11.2f}")
A = np.polyfit(np.log(Ns), ylast, 1)
print(f"fit: ur_last = {A[0]:.3f} ln(N) + {A[1]:.3f};  |gam_tip|/(2 pi) = {2*U0*a70[-1]/(2*np.pi):.3f}")
