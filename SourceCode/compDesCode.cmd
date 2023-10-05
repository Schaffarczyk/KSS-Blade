echo 'compile module mem'
gfortran -c  mem.f
echo 'compile rest'
gfortran -c  -O3  -fno-automatic -fbounds-check KSS.f  Sub1.f  Sub2.f SubNum.f
echo 'link all together'
gfortran mem.o KSS.o  Sub1.o Sub2.o SubNum.f -o KSS.exe
