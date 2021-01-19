
# ionizer

## simulation of ionization by charged particle tracks in silicon

based on  
Hans Bichsel: Straggling in thin silicon detectors  
Rev. Mod. Phys. 60 (1988) 663-99  
http://prola.aps.org/abstract/RMP/v60/i3/p663_1

and  
M. Brigida, M.N. Mazzotta et al:  
A new Monte Carlo code for full simulation of silicon strip detectors  
Nuclear Instruments and Methods in Physics Research A 533 (2004) 322–343
https://doi.org/10.1016/j.nima.2004.05.127

here translated to C++
- ionizer.cpp is without root, produces ion clusters in ionizer.dat
- ionizer.cc with ROOT, produces ionizer.root

example plots: 5 GeV electrons passing 150 um Si

![energy loss](teh.png "ionization in 150 um Si")

![energy loss](xy.png "x-y ionization map")

![energy loss](rz.png "R-z ionization map")
