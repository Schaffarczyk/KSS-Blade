"""Local loads incl. curved tip, grid convergence of the local terms, drag influence (RC2)."""
import json, numpy as np
from driver import run
from decomp import decomp, zfun_of, rs, R0, RT, RROOT, RHO, B
from driver import shape
U0, RPM = 11.0, 12.1; OM = RPM*np.pi/30
out = {}
def zc(rr):   # linear 6 m, outer 5 % additionally inclined to kappa = 30 deg (downwind)
    rr = np.clip(np.asarray(rr, float), R0, RT); z = shape(rr, 'linear', 6.0); rb = 0.95*RT
    return z + np.where(rr > rb, (rr-rb)*(np.tan(np.radians(30))-6.0/(RT-R0)), 0.0)
for lab, zf in (('linear6', zfun_of('linear', 6.0)), ('bend6', zfun_of('bend', 6.0)), ('prebend-6', zfun_of('prebend', -6.0)), ('curvedtip', zc)):
    o = {}
    for n in (70, 140, 280):
        d = decomp(zfun=zf, npform=1, nsec=n)
        rM = run(U0, RPM, 0.0, nsec=n, defl=zf(rs), radial=True, radmode=1, npform=1)
        r = d['r']; G = d['G']
        fr = -RHO*G*d['ur']*d['tank']; fa = -RHO*G*U0*d['da']; fM = -RHO*G*rM['bem']['ur']*d['tank']
        m5 = r > 0.95*RT; m10 = r > 0.9*RT
        mom = lambda f, m: float(B*np.sum((f*r)[m])*d['dr']/1e3)       # kNm
        o[n] = dict(r=r.tolist(), fr=fr.tolist(), fa=fa.tolist(), fM=fM.tolist(), fta=d['res']['bem']['ct'].tolist(),
                    M5=dict(r=mom(fr, m5), a=mom(fa, m5), M=mom(fM, m5)), M10=dict(r=mom(fr, m10), a=mom(fa, m10), M=mom(fM, m10)),
                    Mtot=dict(r=mom(fr, r > 0), a=mom(fa, r > 0), M=mom(fM, r > 0)),
                    peak=dict(r=float(fr[np.abs(fr).argmax()]), M=float(fM[np.abs(fM).argmax()]), net=float((fr+fa)[np.abs(fr+fa).argmax()]),
                              r_at=float(r[np.abs(fr).argmax()]), net_at=float(r[np.abs(fr+fa).argmax()])))
        print(f"{lab:10s} N={n:3d} total kNm: radial VC {o[n]['Mtot']['r']:+7.2f} axial {o[n]['Mtot']['a']:+7.2f} Madsen {o[n]['Mtot']['M']:+7.2f} | outer5%: VC r {o[n]['M5']['r']:+6.2f} a {o[n]['M5']['a']:+6.2f} net {o[n]['M5']['r']+o[n]['M5']['a']:+6.2f} M {o[n]['M5']['M']:+6.2f} | outer10%: r {o[n]['M10']['r']:+6.2f} a {o[n]['M10']['a']:+6.2f} M {o[n]['M10']['M']:+6.2f} | peak f: VC r {o[n]['peak']['r']:+6.1f} @ {o[n]['peak']['r_at']:.1f}, M {o[n]['peak']['M']:+6.1f}, net {o[n]['peak']['net']:+6.1f} @ {o[n]['peak']['net_at']:.1f} N/m")
    out[lab] = o
# aerodynamic in-plane force per unit length for scale, and drag influence
res = run(U0, RPM, 0.0, nsec=70, defl=zfun_of('linear', 6.0)(rs), radial=True, radmode=2, npform=1)
b = res['bem']; phi = b['phi']/57.3
with np.errstate(all='ignore'):
    eps = np.where(np.abs(b['cl']) > 1e-3, b['cd']/b['cl'], 0.0)   # cylindrical root sections (c_L = 0): no reduction applied
ft = 0.5*RHO*b['chord']*b['w']**2*b['ct']
m = b['r'] > 0.25*RT
print('aerodynamic in-plane force per unit length: at r=60.8: %.0f N/m ; max %.0f' % (np.interp(60.8, b['r'], ft), ft.max()))
print('drag share of normal force eps*tan(phi): outer 75%%: %.4f .. %.4f ; eps*cot(phi) (Glauert reduction of wake circulation): %.3f .. %.3f, mean %.3f' %
      ((eps*np.tan(phi))[m].min(), (eps*np.tan(phi))[m].max(), (eps/np.tan(phi))[m].min(), (eps/np.tan(phi))[m].max(), (eps/np.tan(phi))[m].mean()))
# effect of lift-only strengths: a_inf from the lift part of the thrust coefficient
F = b['F']; a = res['vc']['a']; AC = 0.4
ct = np.where(a <= AC, 4*F*a*(1-a), 4*F*(AC*AC+(1-2*AC)*a))
ctL = ct*(b['cl']*np.cos(phi))/b['cn']
ainfL = np.where(ctL <= 0.96, 0.5*(1-np.sqrt(np.maximum(0, 1-ctL))), (0.25*ctL-AC*AC)/(1-2*AC))
from vcfull import sweep, strengths
dr = (RT-RROOT)/70; Rk = RROOT+dr*np.arange(71); zf = zfun_of('linear', 6.0)
tank = 6.0/(RT-R0)*(b['r'] > R0)
for lab2, aa in (('lift+drag', res['vc']['ainf']), ('lift only', ainfL)):
    ur, uz = sweep(b['r'], zf(b['r']), Rk, zf(Rk), strengths(aa, U0)); _, uz0 = sweep(b['r'], 0*b['r'], Rk, 0*Rk, strengths(aa, U0))
    for lab3, G in (('Gamma', b['gam']), ('Gamma(1-eps cot phi)', b['gam']*(1-np.clip(eps/np.tan(phi), 0, 1)))):
        dPr = -B*RHO*OM*np.sum(G*b['r']*ur*tank)*dr/1e3; dPa = B*RHO*OM*np.sum(G*b['r']*(uz-uz0))*dr/1e3
        print(f'   strengths from {lab2:9s}, force with {lab3:22s}: radial {dPr:+7.2f} axial {dPa:+7.2f} net {dPr+dPa:+6.2f} kW ; ur_tip {ur[-1]:.3f}  max|d ur| vs lift+drag')
    out['ur_'+lab2] = ur.tolist()
dur = np.abs(np.array(out['ur_lift only'])-np.array(out['ur_lift+drag']))
for lo in (0.0, 0.25, 0.5):
    mm = b['r'] > lo*RT; print('max |ur(lift only) - ur(lift+drag)| for r > %.2f R: %.3f m/s at r = %.1f m' % (lo, dur[mm].max(), b['r'][mm][dur[mm].argmax()]))
json.dump(out, open('study5.json', 'w'), indent=1)
