## Tau leaping simulations of Glycine moleuces on Cu110



[![MIT licensed](https://img.shields.io/github/license/mashape/apistatus.svg)](http://github.com/jbr36/TauLeapingCode/blob/master/license.md)

[![Release version](1.1)](https://github.com/jbr36/TauLeapingCode/releases/tag/1.1)
https://doi.org/10.5281/zenodo.159466

Software for sequential multi-scale simultions of the kinetics and dynamics of Glycine molecules on Cu(110). The DFT energies used in the lookup table will be published separately.

Below follow some general information about the code:
- Freely available under the MIT License.
- Fortran code, after compilation runs on Windows, Linux and Mac OS X.

#### Compilation
Use the Makefile to compile with gfotran.

#### Input preparation
The input parameters are specified in system.in.

 - ntriangle n       - n number of triangles representing each a Glycine molecule.
 - posNxy n1 n2      - number of copper atoms in x (n1) an y (n2) direction
 - temp t            - t the temperature for the rate calulations from DFT energy barriers  
 - hbondstr hbs      - strength of H-bond between neighboring O and N atoms in Gly
 - nkmcsteps n       - n number of tauleaping steps/ moves of a single triangle until algorithm stops
 - printlevel n      - n number of steps after which a full print out to files should be written

If isl-type1 and isl-type2 are set to 0 and their size is 0 the Gly molecules 
will be randomly distributed over the copper surface.

The triangles have a certain arrangement, as discussed in the paper. 

If type=1 heterochiral clusters will be arranged color 12 34, 
type=2 homochiral clussters are possible color 1 2 3 4 
colors denote: 1 green, 2 red, 3 yellow, 4 blue.

The colours are chosen with color-isl1 and color-isl2 according to the type. 
The implementation allows to have 1 or 2 islands on the surface.

isl1-size-xy and isl2-size-xy specify the size of the clusters in x y directions in terms of number of molecules.

isl1-pos-xy and isl2-pos-xy determine the left lower corner of the Gly clusters on the xy-plane of copper atoms.

isl1-defect-type and isl2-defect-type allow to take out 1 or 2 Gly molecules from the clusters at position
isl1-defects-xy x1 y1 x2 y2 and 
isl2-defects-xy x1 y1 x2 y2,
x1 y1 x2 y2 being the respective coordinates.

#### Defects
Additionally, there can be one or two 'defective' areas on the metal surface, 
where the Gly molecules are not allowed to move to.

mdefect1-type and mdefect2-type turn it with 1 on and with 0 off.

mdefect1-size-xy and mdefect2-size-xy x1 y1 give the number of atoms in x1 and y1 direction which cannot react.
mdefect1-pos-xy and mdefect2-pos-xy determine the left lower corner of the defective areas on the xy-plane of copper atoms.

All results are written to vtk files which can then be visualised in PARAVIEW from http://www.paraview.org/.
The final important feature is the possibility to restart from an input vtk file. 

#### Restart
Add
restartTrue restart.vtk

to system.in

The file restart.vtk copied from a triangle000*.vtk (resulting output) allows to restart from an intermediate output file.

Make sure the system.in is not changed apart from the restartTrue line being added.

More details can be found in the corresponding publication.

