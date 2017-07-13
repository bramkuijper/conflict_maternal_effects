#!/usr/bin/env bash
./xpo_conflict 0.01 0.01       0.2 0.9         0.25  0.1 0.5        0.25 0.25  1.0 0 &
./xpo_conflict 0.01 0.01       0.2 0.9         0.25  0.1 0.5        0.25 0.25  1.0 1 &
#    mu = atof(argv[1]); // mutation prob
#    sdmu = atof(argv[2]); // mutational distribution
#    c[0] = atof(argv[3]); 
#    c[1] = atof(argv[4]); // survival slope
#    p = atof(argv[5]); // initial freq envt 1
#    d = atof(argv[6]); // dispersal 
#    sigma[0] = atof(argv[7]); //prob envt 1 changes
#    sigma[1] = atof(argv[8]); // prob envt 2 changes
#    gam = atof(argv[9]); // prob envt 2 changes
#    maternal_expression = atoi(argv[10]); // dispersal at t=0
