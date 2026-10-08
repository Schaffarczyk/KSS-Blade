"""Revision 2 (RC1/RC2): consistent non-planar formulation (KSS V6, NPform=1),
cylinder strengths from a_inf, complete VC model (radial + axial) at fixed circulation."""
import json, numpy as np
from scipy.interpolate import CubicSpline
from driver import run, shape, STATIONS
from vcfull import sweep, strengths, ur_cyl
from decomp import decomp, zfun_of, rs, R0, RT, RROOT, RHO, B
out = {}
U0, RPM = 11.0, 12.1; OM = RPM*np.pi/30
AC = 0.4

# ---------------- A. reference case: linear 6 m downwind deflection ----------------
zl = zfun_of('linear', 6.0)
pl  = run(U0, RPM, 0.0, radial=False, radmode=0)
A = dict(planar=dict(T=pl['thrust'], P=pl['power'], cP=pl['cP']))
for npf in (0, 1):
    for rm, lab in ((0, 'defl'), (1, 'madsen'), (2, 'vc')):
        r = run(U0, RPM, 0.0, defl=zl(rs), radial=True, radmode=rm, npform=npf)
        A[f'np{npf}_{lab}'] = dict(T=r['thrust'], P=r['power'], cP=r['cP'])
d = decomp(zfun=zl, npform=1)
A['vc_radial_kW'] = d['dPr']; A['vc_axial_kW'] = d['dPa']; A['vc_axial_overF_kW'] = d['dPaF']; A['vc_net_kW'] = d['dPr']+d['dPa']
dB = decomp(zfun=zl, npform=1, use_aB=True)
A['aB_radial_kW'] = dB['dPr']; A['aB_axial_kW'] = dB['dPa']
A['P_complete'] = A['np1_defl']['P'] + d['dPr'] + d['dPa']
A['madsen_radial_kW'] = A['np1_madsen']['P'] - A['np1_defl']['P']
Pref = 0.5*RHO*np.pi*RT**2*U0**3/1e3
A['cP_complete'] = A['P_complete']/Pref
out['A'] = A
print('A:', json.dumps({k: (v if not isinstance(v, dict) else {kk: round(vv, 4) for kk, vv in v.items()}) for k, v in A.items()}, indent=0))

# ---------------- B. local distributions at Nsec = 70 ----------------
rM = run(U0, RPM, 0.0, defl=zl(rs), radial=True, radmode=1, npform=1)
r = d['r']; G = d['G']
ur_aB = dB['ur']
fr  = -RHO*G*d['ur']*d['tank']                # in-plane force per unit radius and blade, radial term [N/m]
fa  = -RHO*G*U0*d['da']                       # axial term
frM = -RHO*G*rM['bem']['ur']*d['tank']
out['B'] = dict(r=r.tolist(), a=d['res']['vc']['a'].tolist(), ainf=d['a'].tolist(), F=d['F'].tolist(), gam=G.tolist(),
                ur_vc=d['ur'].tolist(), ur_vc_aB=ur_aB.tolist(), ur_m=rM['bem']['ur'].tolist(), da=d['da'].tolist(),
                f_r=fr.tolist(), f_a=fa.tolist(), f_rM=frM.tolist())
neg = r[d['ur'] < 0]
print('B: ur tip VC(ainf) %.2f  VC(aB) %.2f  Madsen %.2f ; negative ur between %.1f and %.1f m, min %.2f' %
      (d['ur'][-1], ur_aB[-1], rM['bem']['ur'][-1], neg.min(), neg.max(), d['ur'].min()))
for lo, hi in ((0.7, 0.92), (0.5, 0.9)):
    m = (r/RT > lo) & (r/RT < hi); q = d['ur'][m]/rM['bem']['ur'][m]
    print(f'   ratio VC/Madsen for {lo}<r/R<{hi}: {q.min():.2f} .. {q.max():.2f}')
dr70 = d['dr']
for lab, f in (('radial VC', fr), ('axial VC', fa), ('radial Madsen', frM)):
    print(f'   {lab}: total {B*OM*np.sum(f*r)*dr70/1e3:+.1f} kW, outer 10% {B*OM*np.sum((f*r)[r>0.9*RT])*dr70/1e3:+.1f} kW, inner 90% {B*OM*np.sum((f*r)[r<=0.9*RT])*dr70/1e3:+.1f} kW, peak |f| {np.abs(f).max():.1f} N/m at r={r[np.abs(f).argmax()]:.1f}')
fn = fr+fa
print(f'   net VC: total {B*OM*np.sum(fn*r)*dr70/1e3:+.2f} kW, max |f_net| {np.abs(fn).max():.1f} N/m at r={r[np.abs(fn).argmax()]:.1f}; outer10% {B*OM*np.sum((fn*r)[r>0.9*RT])*dr70/1e3:+.2f} inner90% {B*OM*np.sum((fn*r)[r<=0.9*RT])*dr70/1e3:+.2f}')

# ---------------- C. grid study of the outermost value ----------------
def continuous(res):
    """continuous a_B(r), a_inf(r), Gamma(r) from a KSS solution; tip loss F evaluated analytically"""
    rr = res['bem']['r']; rl = rr[-1]
    _a = CubicSpline(rr, res['vc']['a']); _p = CubicSpline(rr, np.radians(res['bem']['phi']/57.3*180/np.pi)); _g = CubicSpline(rr, res['bem']['gam'])
    aB = lambda x: _a(np.minimum(x, rl)); phi = lambda x: _p(np.minimum(x, rl))
    def F(x):
        xx = x/RT; s = np.sin(phi(x))
        ft = 2/np.pi*np.arccos(np.exp(-0.5*B*(1-xx)/(xx*s)))
        xr = RROOT/RT
        fr_ = 2/np.pi*np.arccos(np.exp(-0.5*B*np.maximum(xx-xr, 0)/(xx*s)))
        return ft*fr_
    def ainf(x):
        a = aB(x); ct = np.where(a <= AC, 4*F(x)*a*(1-a), 4*F(x)*(AC*AC+(1-2*AC)*a))
        return np.where(ct <= 4*AC*(1-AC), 0.5*(1-np.sqrt(np.maximum(0, 1-ct))), (0.25*ct-AC*AC)/(1-2*AC))
    gam = lambda x: np.where(x <= rl, _g(np.minimum(x, rl)), _g(rl)*np.clip((RT-x)/(RT-rl), 0, 1))
    return aB, ainf, gam, F
# frozen solution of the grid study: Nsec = 70 (as in the discussion paper and in Fig. 2a)
res70 = run(U0, RPM, 0.0, nsec=70, defl=zl(rs), radial=True, radmode=2, npform=1)
aBf, ainff, gamf, Ff = continuous(res70)
chk = np.max(np.abs(ainff(res70['bem']['r']) - res70['vc']['ainf']))
print('C: max |a_inf(analytic F) - a_inf(KSS)| at Nsec=70 midpoints: %.1e' % chk)
C = dict(nsec=[], ur_tip_ainf=[], ur_tip_aB=[], ur_fix_ainf=[], ur_fix_aB=[], ur_tip_ainf_guard=[], gam_tip_ainf=[], gam_tip_aB=[], prof={})
rfix = 61.0
def vc_ur(n, afun, eval_r=None, kcap=None):
    dr = (RT-RROOT)/n; rm = RROOT + dr*(np.arange(n)+0.5); Rk = RROOT + dr*np.arange(n+1)
    gam = strengths(afun(rm), U0); re = rm if eval_r is None else np.asarray(eval_r, float)
    ur = np.zeros_like(re)
    for R, g in zip(Rk, gam):
        if abs(g) < 1e-14 or R < 1e-6: continue
        z = zl(re)-zl(R); k2 = 4*re*R/((R+re)**2+z*z)
        if kcap: k2 = np.minimum(k2, kcap)
        from scipy.special import ellipk, ellipe
        k = np.sqrt(k2); ur += -g/(2*np.pi)*np.sqrt(R/re)*((2-k2)/k*ellipk(k2) - 2/k*ellipe(k2))
    return rm, ur, gam
for n in (35, 70, 140, 280, 560, 1120, 2240, 4480):
    rm, u1, g1 = vc_ur(n, ainff); _, u2, g2 = vc_ur(n, aBf)
    _, u1g, _ = vc_ur(n, ainff, kcap=1-1e-6)
    f1 = [np.interp(rfix, rm, u1)]; f2 = [np.interp(rfix, rm, u2)]   # fixed radius: interpolated between section midpoints
    C['nsec'].append(n); C['ur_tip_ainf'].append(float(u1[-1])); C['ur_tip_aB'].append(float(u2[-1])); C['ur_tip_ainf_guard'].append(float(u1g[-1]))
    C['ur_fix_ainf'].append(float(f1[0])); C['ur_fix_aB'].append(float(f2[0])); C['gam_tip_ainf'].append(float(g1[-1])); C['gam_tip_aB'].append(float(g2[-1]))
    if n in (35, 70, 280, 2240):
        m = rm > 0.8*RT; C['prof'][str(n)] = dict(r=rm[m].tolist(), ur_ainf=u1[m].tolist(), ur_aB=u2[m].tolist())
    print(f'   N={n:5d}  ur_tip a_inf {u1[-1]:6.3f} (guard {u1g[-1]:6.3f}) a_B {u2[-1]:6.3f} | at r={rfix}: a_inf {f1[0]:6.3f} a_B {f2[0]:6.3f} | gam_tip/2pi a_inf {abs(g1[-1])/2/np.pi:.3f} a_B {abs(g2[-1])/2/np.pi:.3f}')
# re-converged KSS runs
C['kss'] = []
for n in (35, 70, 140, 280):
    rr = run(U0, RPM, 0.0, nsec=n, defl=zl(rs), radial=True, radmode=2, npform=1)
    r0 = run(U0, RPM, 0.0, nsec=n, defl=zl(rs), radial=True, radmode=0, npform=1)
    dd = decomp(zfun=zl, npform=1, nsec=n, res=rr)
    C['kss'].append(dict(nsec=n, ur_tip=float(rr['vc']['ur'][-1]), P_defl=r0['power'], dPr=dd['dPr'], dPa=dd['dPa'], ainf_tip=float(rr['vc']['ainf'][-1])))
    print(f"   KSS N={n}: ur_tip {rr['vc']['ur'][-1]:.3f}  a_inf(tip) {rr['vc']['ainf'][-1]:.4f}  P_defl {r0['power']:.1f}  dPr {dd['dPr']:+.2f}  dPa {dd['dPa']:+.2f}  net {dd['dPr']+dd['dPa']:+.2f}")
lnN = np.log(C['nsec'][:7]); pf = np.polyfit(lnN, C['ur_tip_aB'][:7], 1); C['fit_aB'] = pf.tolist()
print('   log fit a_B strengths: slope %.3f const %.3f' % tuple(pf))
out['C'] = C
json.dump(out, open('study3.json', 'w'), indent=1)
