"""Power decomposition of the complete non-planar VC model at fixed circulation (cf. RC1, Eq. 11)."""
import numpy as np
from driver import run, shape, STATIONS
from vcfull import sweep, strengths
rs = np.array([s[0] for s in STATIONS]); R0, RT, RROOT = 2.8667, 63.0, 0.3
RHO, B = 1.225, 3

def zfun_of(kind, d):
    return lambda rr: shape(np.clip(rr, R0, RT), kind, d)

def decomp(vw=11.0, rpm=12.1, pitch=0.0, nsec=70, zfun=None, npform=1, use_aB=False, res=None):
    """returns dict with radial and axial power contributions [kW] at fixed Gamma of the KSS solution"""
    om = rpm*np.pi/30
    if res is None:
        res = run(vw, rpm, pitch, nsec=nsec, defl=zfun(rs), radial=True, radmode=2, npform=npform)
    r = res['vc']['r']; G = res['bem']['gam']; F = res['bem']['F']
    a = res['vc']['a'] if use_aB else res['vc']['ainf']
    dr = (RT-RROOT)/nsec
    Rk = RROOT + dr*np.arange(nsec+1)
    gam = strengths(a, vw)
    zm, zk = zfun(r), zfun(Rk)
    eps = 1e-3
    tank = (zfun(r+eps)-zfun(r-eps))/(2*eps)          # positive downwind
    ur, uz = sweep(r, zm, Rk, zk, gam)
    _, uz0 = sweep(r, 0*zm, Rk, 0*zk, gam)
    da = -(uz-uz0)/vw                                   # change of axial induction, non-planar minus planar
    dPr = -B*RHO*om*np.sum(G*r*ur*tank)*dr/1e3          # radial term (downwind-positive kappa)
    dPa = -B*RHO*om*np.sum(G*r*vw*da)*dr/1e3            # axial term
    dPaF = -B*RHO*om*np.sum(G*r*vw*da/np.maximum(F,1e-3))*dr/1e3
    return dict(r=r, G=G, F=F, a=a, ur=ur, da=da, tank=tank, dPr=dPr, dPa=dPa, dPaF=dPaF, res=res, om=om, dr=dr,
                a_planar_check=float(np.max(np.abs(-uz0/vw - a))))

if __name__ == '__main__':
    for lab, npf, aB in (('consistent, a_inf', 1, False), ('consistent, a_B', 1, True), ('V5 solution, a_inf',0,False)):
        d = decomp(zfun=zfun_of('linear', 6.0), npform=npf, use_aB=aB)
        print(f"{lab:22s} radial {d['dPr']:+8.2f} kW  axial {d['dPa']:+8.2f} kW (with 1/F: {d['dPaF']:+8.2f})  sum {d['dPr']+d['dPa']:+7.2f}   P={d['res']['power']:.1f} planar-check {d['a_planar_check']:.1e}")
