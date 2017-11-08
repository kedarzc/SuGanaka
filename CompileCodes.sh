#!/bin/bash

# Delete all the previous files
rm *.mod
rm *.o

# Compile the code
gfortran -c MainModule.f90
gfortran CST_main.f90 -o CST_main MainModule.o

# Run the Code
./CST_main 
