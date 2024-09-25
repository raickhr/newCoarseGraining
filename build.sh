FC=ifort

$FC -c configurationMod.F90
$FC -c main.F90 

$FC configurationMod.o main.o -o main.exe 

rm -f *.o *.mod