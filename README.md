SuGanaka
--------
A Structural Finite Element Software. SuGanaka is a Sanskrit word
meaning 'Good Computer' (Su: Good, Ganaka:Someone who computes). The code 
is in its infancy but can simulate 2D structural problems using the 3-Node
Plane Stress Triangular elements (CST).

Here is a list of the files and the data stored in each of them.
1. Connectivity_ABQ.csv : Stores the nodal connectivity
2. MeshData_ABQ.csv : Stores the nodal information
3. CST_main.f90 : Main Source code for SuGanaka
4. MainModule.f90: A Module which stores various subroutines
5. Results.sdb: sbb->SuGanaka Data Base, stores nodal displacements and element stresses and strains
6. Results.vtk: A vtk file which can be used for visualization in Paraview
7. CompileCodes.sh: A shell script to compile the codes and run it

How to run the code
--------------------
Step1: chmod CompileCodes.sh
Step2: Run the code ./CompileCodes


