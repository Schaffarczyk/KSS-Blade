# KSS-Blade
Aerodynamic Design Code for Wind Turbine Blades

Korjahn (bewind) Schlipf (WETI), Schaffarczyk (KUAS)

V2: 2022 July 6

V3: 2022 Oct  4
two sample cases have been added:
a) NREL 5MW 
b) IEAwind/NREL 15 MW

V4: 2025 Oct 08
updated for Optimus Syria Project

V5: 2026 Mar 25
MSc Thesis Bhima Masare
implement "large deflections"

V6: 2026 Oct 08
consistent formulation for non-planar rotors (deflected or pre-bent blades),
following the review of the WES brief communication wes-2026-138:
- new line "NPform" in Machine.in: 1 = consistent formulation (default), 0 = formulation of V5
- axial momentum balance with cos(kappa)**2, torque integrated over the blade length
- corrected sign of the radial-velocity force (oop deflection is positive downwind)
- vortex cylinder strengths from the annulus induction a_inf (written to VC.out)
- kappa from a central difference about the section midpoint
- NREL-5MW: CYL1.aer and CYL2.aer (cylindrical root sections) added, they were missing
- folder WES-2026-138: Python scripts and results of the studies for the WES paper

Details and reference results: ManDoc/README-V6.txt
