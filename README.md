# IbexHomeo
a specialized algorithm for finding homeostasis using interval methods.

To install IbexHomeo, you need first to install Ibex (http://www.ibex-lib.org/download), then download theses files and install them, for examples in a sub-repository of examples.
To compile these files, it is recommanded to use gcc/7.5.0 on linux.
Then, you can follow the procedure given here : http://www.ibex-lib.org/doc/install.html#compiling-a-test-program applied to ibexoptbox4.

It take as input a minibex file (http://www.ibex-lib.org/doc/minibex.html) where your parameters are given in the Constants part of the file. The ones that vary will be added to the variables set by a preprocessing.
Options of IbexHomeo are the same as for IbexOpt (http://www.ibex-lib.org/doc/optim.html), plus an option --homeostasis. The option -t is for the loop timeout.
