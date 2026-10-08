"""Figures of the second revision: fig1-local.pdf, fig2-grid-param.pdf"""
import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker
from scipy.interpolate import CubicSpline
from scipy.special import ellipk, ellipe
from driver import run
from vcfull import strengths
from decomp import zfun_of, rs, RT, RROOT, B
plt.rcParams.update({'axes.unicode_minus': False, 'font.size': 8.5, 'axes.linewidth': 0.6, 'axes.spines.top': False, 'axes.spines.right': False,
                     'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'legend.frameon': False, 'pdf.use14corefonts': True, 'pdf.compression': 9, 'font.family': 'sans-serif', 'font.sans-serif': ['Helvetica'], 'mathtext.fontset': 'custom', 'mathtext.rm': 'Helvetica', 'mathtext.it': 'Helvetica:italic', 'mathtext.bf': 'Helvetica:bold'})
CB, CO, CG, CK, CGR = '#0072B2', '#D55E00', '#009E73', '#222222', '#8a8a8a'
S3 = json.load(open('study3.json')); S4 = json.load(open('study4.json'))
Bd = S3['B']; r = np.array(Bd['r'])

# ---------------- Fig. 1 ----------------
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
a = ax[0]
a.axhline(0, color='#bbbbbb', lw=0.5)
a.plot(r, Bd['ur_m'], color=CB, lw=1.6, label='Madsen')
a.plot(r, Bd['ur_vc'], color=CO, lw=1.6, label='vortex cylinder, strengths from annulus induction')
a.plot(r, Bd['ur_vc_aB'], color=CGR, lw=1.1, ls='--', label='vortex cylinder, strengths from blade induction')
a.set_xlabel('r / m'); a.set_ylabel('radial induced velocity / (m/s)'); a.set_xlim(0, 63); a.set_ylim(-1.6, 7)
a.legend(loc='upper left', fontsize=7.2, handlelength=2.2)
a.text(0.02, 0.03, '(a)', transform=a.transAxes)
a = ax[1]
a.axhline(0, color='#bbbbbb', lw=0.5)
fr, fa, fM = np.array(Bd['f_r']), np.array(Bd['f_a']), np.array(Bd['f_rM'])
a.plot(r, fM, color=CB, lw=1.6, label='Madsen, radial term')
a.plot(r, fr, color=CO, lw=1.6, label='VC, radial term')
a.plot(r, fa, color=CG, lw=1.6, ls='--', label='VC, axial term')
a.plot(r, fr+fa, color=CK, lw=1.3, ls='-.', label='VC, radial + axial')
a.set_xlabel('r / m'); a.set_ylabel('additional in-plane force / (N/m)'); a.set_xlim(0, 63)
a.legend(loc='lower left', fontsize=7.2, handlelength=2.6)
a.text(0.93, 0.92, '(b)', transform=a.transAxes)
fig.tight_layout(w_pad=2.0); fig.savefig('fig1-local.pdf'); fig.savefig('fig1-local.png', dpi=170); plt.close(fig)

# ---------------- Fig. 2 ----------------
U0, RPM, AC = 11.0, 12.1, 0.4
zl = zfun_of('linear', 6.0)
res = run(U0, RPM, 0.0, nsec=70, defl=zl(rs), radial=True, radmode=2, npform=1)
rr = res['bem']['r']; rl = rr[-1]
_a = CubicSpline(rr, res['vc']['a']); _p = CubicSpline(rr, res['bem']['phi']/57.3)
aB = lambda x: _a(np.minimum(x, rl)); phi = lambda x: _p(np.minimum(x, rl))
def F(x):
    xx = x/RT; s = np.sin(phi(x))
    return 2/np.pi*np.arccos(np.exp(-0.5*B*(1-xx)/(xx*s))) * 2/np.pi*np.arccos(np.exp(-0.5*B*np.maximum(xx-RROOT/RT, 0)/(xx*s)))
def ainf(x):
    a_ = aB(x); ct = np.where(a_ <= AC, 4*F(x)*a_*(1-a_), 4*F(x)*(AC*AC+(1-2*AC)*a_))
    return np.where(ct <= 4*AC*(1-AC), 0.5*(1-np.sqrt(np.maximum(0, 1-ct))), (0.25*ct-AC*AC)/(1-2*AC))
def tip(n, afun, kcap=None):
    dr = (RT-RROOT)/n; rm = RROOT+dr*(np.arange(n)+0.5); Rk = RROOT+dr*np.arange(n+1); g = strengths(afun(rm), U0)
    re = rm[-1]; z = zl(re)-zl(Rk); k2 = 4*re*Rk/((Rk+re)**2+z*z)
    if kcap: k2 = np.minimum(k2, kcap)
    m = (np.abs(g) > 1e-14) & (Rk > 1e-6); k = np.sqrt(k2[m])
    return float(np.sum(-g[m]/(2*np.pi)*np.sqrt(Rk[m]/re)*((2-k2[m])/k*ellipk(k2[m]) - 2/k*ellipe(k2[m])))), abs(g[-1])/(2*np.pi)
NS = [35, 70, 140, 280, 560, 1120, 2240, 4480]
tA = [tip(n, ainf)[0] for n in NS]; tB = [tip(n, aB)[0] for n in NS]; tBg = [tip(n, aB, 1-1e-6)[0] for n in NS]
pf = np.polyfit(np.log(NS[:7]), tB[:7], 1)
G2 = dict(nsec=NS, tip_ainf=tA, tip_aB=tB, tip_aB_guard=tBg, fit_aB=pf.tolist(), gam_ainf=[tip(n, ainf)[1] for n in NS], gam_aB=[tip(n, aB)[1] for n in NS],
          kss=[(k['nsec'], k['ur_tip']) for k in S3['C']['kss']], madsen=Bd['ur_m'][-1])
json.dump(G2, open('fig2-data.json', 'w'), indent=1)
print('tip a_inf:', np.round(tA, 3)); print('tip a_B  :', np.round(tB, 3)); print('guard    :', np.round(tBg, 3)); print('fit a_B slope, const:', pf, ' gam_tip/2pi a_B:', G2['gam_aB'][1], 'a_inf', np.round(G2['gam_ainf'], 3))

fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw={'width_ratios': [1, 1.25]})
a = ax[0]
a.plot(NS, tB, 's', color=CGR, ms=4.5, label='strengths from blade induction')
xx = np.array([30, 5200]); a.plot(xx, pf[0]*np.log(xx)+pf[1], ':', color=CGR, lw=0.9)
a.plot(NS, tA, 'o-', color=CO, ms=4.5, lw=1.2, label='strengths from annulus induction')
kn, kv = zip(*G2['kss']); a.plot(kn, kv, 'D', mfc='none', mec=CK, ms=6.5, mew=0.9, label='KSS-Blade, re-converged')
a.axhline(G2['madsen'], color=CB, lw=1.4); a.text(36, G2['madsen']+0.15, 'Madsen', color=CK, fontsize=7.5)
a.set_xlabel('number of sections'); a.set_ylabel('radial velocity at outermost midpoint / (m/s)'); a.set_ylim(1.5, 12.0); a.set_xscale('log'); a.set_xlim(28, 5600); a.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([35, 70, 140, 280, 560, 1120, 2240, 4480])); a.xaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(['35', '70', '140', '280', '560', '1120', '2240', '4480'])); a.xaxis.set_minor_locator(matplotlib.ticker.NullLocator()); a.tick_params(axis='x', labelsize=7)
a.legend(loc='upper left', fontsize=7.2); a.text(0.92, 0.13, '(a)', transform=a.transAxes)
a = ax[1]
rows = {(x['shape']): x for x in S4 if x['vw'] == 11.0 and abs(x['d']) == 6.0}
labs = [('linear', 'linear\n6 m downwind'), ('bend', 'bending line\n6 m downwind'), ('prebend', 'prebend\n6 m upwind')]
keys = [('rel_dP_M', 'Madsen, radial term', CB, None), ('rel_dP_r', 'VC, radial term', CO, None), ('rel_dP_a', 'VC, axial term', CG, '//'), ('rel_net', 'VC, radial + axial', CK, None)]
w = 0.19
a.axhline(0, color='#888888', lw=0.6)
for j, (k, lab, c, h) in enumerate(keys):
    xs = np.arange(3) + (j-1.5)*(w+0.012); vals = [rows[s][k] for s, _ in labs]
    a.bar(xs, vals, w, color=c if h is None else 'white', edgecolor=c, hatch=h, lw=0.8, label=lab)
    for x_, v in zip(xs, vals):
        a.text(x_, v + (0.12 if v >= 0 else -0.12), f'{v:+.2f}', ha='center', va='bottom' if v >= 0 else 'top', fontsize=6.3, color=CK)
a.set_xticks(np.arange(3)); a.set_xticklabels([l for _, l in labs]); a.tick_params(axis='x', length=0)
a.set_ylabel('power change / % of deflected rotor'); a.set_ylim(-6.3, 8.6)
a.legend(loc='upper left', fontsize=7.2, ncol=2, columnspacing=1.0, handlelength=1.5); a.text(0.95, 0.04, '(b)', transform=a.transAxes)
fig.tight_layout(w_pad=2.0); fig.savefig('fig2-grid-param.pdf'); fig.savefig('fig2-grid-param.png', dpi=170); plt.close(fig)
