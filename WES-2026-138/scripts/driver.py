"""Driver for KSS V3 (RadMode 0/1/2) parametric runs."""
import os, re, shutil, subprocess, tempfile
import numpy as np

SRC = os.environ.get('KSS_SRC', '../KSS-Blade-V6')   # folder with the compiled binary 'kss' and the input files
BLADE = """r\tTwist   Chord  Pr-Name  oop-Deflection
"""
# NREL 5 MW stations as in BlaDes.in (r, twist, chord, profile)
STATIONS = [
 (2.8667,13.308,3.542,'CYL1'),(5.6000,13.308,3.854,'CYL1'),(8.3333,13.308,4.167,'CYL2'),
 (11.7500,13.308,4.557,'DU40'),(15.8500,11.480,4.652,'DU35'),(19.9500,10.162,4.458,'DU35'),
 (24.0500,9.011,4.249,'DU30'),(28.1500,7.795,4.007,'DU25'),(32.2500,6.544,3.748,'DU25'),
 (36.3500,5.361,3.502,'DU21'),(40.4500,4.188,3.256,'DU21'),(44.5500,3.125,3.010,'NA17'),
 (48.6500,2.319,2.764,'NA17'),(52.7500,1.526,2.518,'NA17'),(56.1667,0.863,2.313,'NA17'),
 (58.9000,0.370,2.086,'NA17'),(60.0000,0.200,2.000,'NA17'),(61.6333,0.106,1.419,'NA17'),
 (62.3170,0.530,0.714,'NA17'),(63.0000,0.000,0.010,'NA17')]
R0, RT = 2.8667, 63.0

def shape(r, kind, d):
    xi = (np.asarray(r) - R0) / (RT - R0)
    if kind == 'linear':   z = xi
    elif kind == 'bend':   z = xi**2 * (3 - xi) / 2          # cantilever bending line, tip slope 1.5 d/L
    elif kind == 'prebend': z = xi**3                        # outboard-concentrated prebend, tip slope 3 d/L
    elif kind == 'flat':   z = 0*xi
    else: raise ValueError(kind)
    return d * z

def machine(vw, rpm, pitch, nsec, radial, radmode, npform=1):
    return f"""DesMode    .F.
BladeNo    3.0
dens       1.225
Rhub       0.3
Rtip       63.0
windrage   {vw} {vw} 1
rpm        {rpm}
Pitchr     {pitch} {pitch} 1
Nsec        {nsec}
tipshape   P
Tiplen      3.0
twmax      20.
chmax       5.5
ImpChord   .F.
ImpTwist   .F. 7
ImpThick   .F.
TwistB     .F.
DesSchema  2
chroot     3.5
Radial     {'.T.' if radial else '.F.'}
RadMode    {radmode}
NPform     {npform}
"""

def run(vw=11.0, rpm=12.1, pitch=0.0, nsec=70, radial=True, radmode=2, defl=None, keep=None, npform=1):
    """defl: array of oop deflection at STATIONS (m). Returns dict."""
    if defl is None: defl = np.zeros(len(STATIONS))
    wd = tempfile.mkdtemp(prefix='kss_')
    for f in os.listdir(SRC):
        if f.endswith('.aer') or f in ('ProThick.in','ThickDis.in','kss'):
            shutil.copy(os.path.join(SRC,f), wd)
    with open(os.path.join(wd,'Machine.in'),'w') as f: f.write(machine(vw,rpm,pitch,nsec,radial,radmode,npform))
    with open(os.path.join(wd,'BlaDes.in'),'w') as f:
        f.write(BLADE)
        for (r,tw,ch,pr),z in zip(STATIONS,defl):
            f.write(f"{r:.4f}\t{tw:.3f}   {ch:.3f}     {pr}           {z:.4f}\n")
    out = subprocess.run(['./kss'], cwd=wd, capture_output=True, text=True, timeout=120).stdout
    g = lambda k: float(re.search(k+r'\s*=\s*([-0-9.Ee+]+)', out).group(1))
    res = dict(thrust=g('Thrust'), torque=g('Torque'), power=g('Power '), cP=g('cP'), cT=g('cT'))
    m = re.search(r'VC dTorque=\s*([-0-9.Ee+]+)', out)
    res['dtorvc'] = float(m.group(1)) if m else 0.0
    bem = np.genfromtxt(os.path.join(wd,'Bem.out'), skip_header=1, max_rows=nsec, delimiter=[8]*14+[9]*3+[8,9,9,5,3]+[9]*5, usecols=(0,1,2,3,4,5,16,24,25,26,7,9,10,12,13))
    res['bem'] = dict(r=bem[:,0], a=bem[:,1], ap=bem[:,2], F=bem[:,3], w=bem[:,4], chord=bem[:,5],
                      gam=bem[:,6], kappa=bem[:,7], ctloc=bem[:,8], ur=bem[:,9], phi=bem[:,10], cl=bem[:,11], cd=bem[:,12], cn=bem[:,13], ct=bem[:,14])
    vc = os.path.join(wd,'VC.out')
    if os.path.exists(vc) and radmode==2:
        v = np.genfromtxt(vc, skip_header=1)
        res['vc'] = dict(r=v[:,0], a=v[:,1], kappa=v[:,2], ur=v[:,3], dFz=v[:,4], ainf=v[:,5])
    if keep: shutil.copytree(wd, keep, dirs_exist_ok=True)
    shutil.rmtree(wd)
    return res

