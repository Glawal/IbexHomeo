# IbexHomeo
a specialized algorithm for finding homeostasis using interval methods.

To install IbexHomeo, you first to install Ibex (http://www.ibex-lib.org/download), then download theses files and install them, for examples in a sub-repository of examples, and compile them with the given makefile.

It take as input a minibex file (http://www.ibex-lib.org/doc/minibex.html) where your parameters are given in the Constants part of the file. The ones that vary will be added to the variables set by a preprocessing.
Options of IbexHomeo are the same as for IbexOpt (http://www.ibex-lib.org/doc/optim.html), plus an option --homeostasis. The option -t is for the loop timeout.
