echo 'compile KSS V6 (non-planar rotors, radial induction)'
gfortran -std=legacy -fno-automatic -O1 -fbounds-check -o KSS.exe mem.f KSS.f  Sub1.f  Sub2.f Sub3.f SubNum.f
echo 'ready'
