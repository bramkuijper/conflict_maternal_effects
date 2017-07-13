#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# script to get mathematica output into the *.cpp file

import os, re


# regular expressions to replace variables 
regexps = [ 
        (r"γ\((\d)\)","gam\\1"), 
        (r"β\((\d)\)","beta\\1"), 
        (r"s\((\d)\)", "s\\1"), 
        (r"c\((\d)\)", "c\\1"), 
        (r"p\((\d)\)", "p\\1"), 
        (r"([u|v])\((\d)\)", "\\1\\2"), 
        (r"σ\((\d),(\d)\)", "sigma\\1\\2"), ("Power","pow"),("Sqrt","sqrt"), (r"λ","ev"), (r"List","") ]

# open base *.cpp file
f = open("resolution_numerical.cpp")
fl = f.read()
f.close()

# open eigenvectors file
f = open("eigenvectors.txt")
fl_eigenvec = f.read()
f.close()


# open eigenvectors file
f = open("eigenvectors_battle.txt")
fl_eigenvec_battle = f.read()
f.close()


# selection gradients 
f = open("selgrad.txt")
fl_selgrad = f.read()
f.close()

# selection gradients 
f = open("selgrad_offspring_battle.txt")
fl_selgrad_offspring_battle = f.read()
f.close()


# selection gradients 
f = open("selgrad_mother_battle.txt")
fl_selgrad_mother_battle = f.read()
f.close()

fl = re.sub("EIGENVECS_RESOLUTION",fl_eigenvec, fl) 
fl = re.sub("EIGENVECS_BATTLE",fl_eigenvec_battle, fl) 
fl = re.sub("SELGRADS_RESOLUTION",fl_selgrad, fl) 
fl = re.sub("SELGRADS_MOM",fl_selgrad_mother_battle, fl) 
fl = re.sub("SELGRADS_OFF",fl_selgrad_offspring_battle, fl) 

# change all indices etc to a c-able format
for regexp in regexps:
    fl = re.sub(regexp[0], regexp[1], fl)

f = open("resolution_numerical2.cpp","w")
f.write(fl)
f.close()
