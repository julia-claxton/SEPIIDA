Code for EPP-induced radiation forward propagation

Forked from grant's model. Readme TBD until code is finished.

need to edit `EDIT_THIS_FILE.mac` with your build directory. all user-editable simulation parameters also live there so take a look before running

How to build executable:
... TODO

```
cd "path/to/SEPIIDA_source"

cmake -DCMAKE_INSTALL_PREFIX="/path/to/geant4-install" -DGeant4_DIR="/path/to/geant4-install/lib" "path/to/SEPIIDA_source"

make
```

How to run executable:

`./SEPIIDA <number of particles> <particle> <energy> <pitch angle>`

particle:
* e- = electrons
* proton = protons
* gamma = photon
* For other particles, see: https://fismed.ciemat.es/GAMOS/GAMOS_doc/GAMOS.5.1.0/x11519.html

energy in kev

pitch angle in degrees

**TODO**
add columns of msis file

control brem splitting from command line
