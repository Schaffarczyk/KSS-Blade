echo 'compile module mem'
gfortran -c  mem.f
echo 'compile rest'
gfortran -c  -O3  -fno-automatic -fbounds-check KSS.f  Sub.f SubNum.f
echo 'link all together'
gfortran mem.o KSS.o  Sub.o SubNum.f -o KSS.exe
