#!bin/bash
cd "/home/alessio/cpp/cpp-training-katz/ex7"
g++ ex7.cpp -g -o ex7.exe && "/home/alessio/cpp/cpp-training-katz/ex7/"ex7.exe
cd "/home/alessio/cpp/cpp-training-katz/ex7" 
rm -f "./CPLV.DAT"
gfortran ex7.f -o ex7_fortran.exe && "/home/alessio/cpp/cpp-training-katz/ex7/"ex7_fortran.exe

