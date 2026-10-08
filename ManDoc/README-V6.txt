================================================================
 KSS-Blade  V6                                      2026-10-08
 Consistent non-planar formulation (reply to RC1 on wes-2026-138)
================================================================
V6 corrects four points of the non-planar formulation of V5 that
were identified by Referee 1 of the WES discussion (source-level
analysis), plus one further inaccuracy found during the revision.

New line in Machine.in (after RadMode):
   NPform  1    consistent formulation (default if line is missing)
   NPform  0    formulation of V5 (reproduces the discussion paper)

Changes for NPform = 1
 1. Sub1.f  axial momentum balance of the annulus with cos(kappa)**2
            a/(1-a) = B c Cn cos^2(kappa) / (8 pi F r sin^2(phi))
            (fixed-point iteration, Glauert/Hansen branch, Newton)
 2. Sub1.f  torque integrated over the blade length ds = dr/cos(kappa)
 3. Sub1.f  sign of the radial-velocity force: the oop deflection in
            BlaDes.in is positive DOWNWIND, an outward u_r then
            REDUCES the driving force:
            dFzinc = - dens*Gamma*u_r*tan(kappa)*dr
            (Li et al. define kappa > 0 for UPWIND dihedral)
 4. Sub1.f/Sub3.f  vortex cylinder strengths from the ANNULUS
            induction a_inf (no tip loss) instead of the blade
            induction a_B:   4 a_inf (1-a_inf) = 4 F a_B (1-a_B)
            (Glauert branch as for a_B); a_inf is written to VC.out
 5. Sub1.f  kappa from a central difference about the section
            midpoint (V5: backward difference, lagging by dr/2)

RadMode 0 (deflected geometry, no u_r term), 1 (Madsen), 2 (VC)
as before. The axial induction change of the non-planar VC model is
NOT part of the Fortran code; it is evaluated in post-processing
(scripts/vcfull.py, scripts/decomp.py) at fixed circulation.

Reference case (NREL 5 MW modified, 11 m/s, 12.1 rpm, 6 m linear
downwind deflection, Nsec 70):
                      NPform 0 (=V5)     NPform 1 (V6)
   planar              699.77 kN  4718.40 kW   (same)
   RadMode 0           695.12     4650.46      697.92   4716.52
   RadMode 1 Madsen    695.12     4726.29      697.92   4640.08
   RadMode 2 VC        695.12     4728.71      697.92   4620.54
   complete VC (post-processing, radial -96.0 + axial +95.3 kW): 4715.85 kW

Build:
  gfortran -std=legacy -fno-automatic -O1 -o kss \
           mem.f KSS.f Sub1.f Sub2.f Sub3.f SubNum.f

Note: CYL1.aer and CYL2.aer (cylindrical root sections, cD = 1.0)
are required for the NREL-5MW case; they were missing in the
NREL-5MW folder of the GitHub repository (V5).
================================================================
