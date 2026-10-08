"""Parametric study (consistent formulation): radial-only (Madsen, VC) and complete VC model at fixed circulation."""
import json, numpy as np
from driver import run, shape, STATIONS
from decomp import decomp, zfun_of, rs, RT, RROOT, RHO, B
OPS = [(7.0, 8.0, 0.0), (9.0, 10.3, 0.0), (11.0, 12.1, 0.0), (13.0, 12.1, 6.6), (15.0, 12.1, 10.45)]
# deflection positive downwind; the prebend shape is applied UPWIND (negative)
CASES = [('linear', +1), ('bend', +1), ('prebend', -1)]
rows = []
for vw, rpm, pitch in OPS:
    base = run(vw, rpm, pitch, radial=False, radmode=0)
    for kind, sgn in CASES:
        for d in (3.0, 6.0, 9.0):
            zf = zfun_of(kind, sgn*d)
            r0 = run(vw, rpm, pitch, defl=zf(rs), radial=True, radmode=0, npform=1)
            r1 = run(vw, rpm, pitch, defl=zf(rs), radial=True, radmode=1, npform=1)
            dd = decomp(vw, rpm, pitch, zfun=zf, npform=1)
            row = dict(vw=vw, rpm=rpm, pitch=pitch, shape=kind, d=sgn*d, cT=base['cT'], P_planar=base['power'], T_planar=base['thrust'],
                       P_defl=r0['power'], T_defl=r0['thrust'], dP_M=r1['power']-r0['power'], dP_r=dd['dPr'], dP_a=dd['dPa'],
                       dP_kss_vc=dd['res']['power']-r0['power'])
            row['net'] = row['dP_r']+row['dP_a']
            for k in ('dP_M', 'dP_r', 'dP_a', 'net'): row['rel_'+k] = row[k]/row['P_defl']*100
            row['rel_geo'] = (row['P_defl']-row['P_planar'])/row['P_planar']*100
            row['diff_pct'] = (row['dP_r']-row['dP_M'])/abs(row['dP_M'])*100
            rows.append(row)
            print(f"{vw:4.0f} {kind:8s} {sgn*d:+3.0f} cT={base['cT']:.3f} Pdefl-Ppl={row['rel_geo']:+5.2f}%  radial M {row['rel_dP_M']:+5.2f}%  VC {row['rel_dP_r']:+5.2f}% (VC/M-1 = {row['diff_pct']:+5.1f}%)  axial {row['rel_dP_a']:+5.2f}%  net {row['rel_net']:+6.3f}%  [{row['net']:+6.2f} kW]  kss-check {row['dP_kss_vc']-row['dP_r']:+.2f}", flush=True)
json.dump(rows, open('study4.json', 'w'), indent=1)
