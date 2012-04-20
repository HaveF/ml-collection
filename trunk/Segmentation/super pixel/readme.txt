Running Superpixel
==================
download superpixel http://www.cs.sfu.ca/~mori/research/superpixels
download boundary detector http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/segbench/code/segbench.tar.gz

untar superpixel code to $SUPERPIXEL
untar segbench.tar.gz to $SEGBENCH
cd $SUPERPIXEL/yu_imncut
compile all *.c files by command "mex *.c"
cd ..
All ./yu_imncut and all folders under $SEGBENCH to system path
run sp_demo

Running HOG
===========
Please download the code from 
http://pascal.inrialpes.fr/soft/olt
and follow the instructions.

Things and stuff
================
Please download the code from
http://ai.stanford.edu/~gaheitz/Research/TAS/tas.tgz

Untar the file and follow the instructions in the readme file.