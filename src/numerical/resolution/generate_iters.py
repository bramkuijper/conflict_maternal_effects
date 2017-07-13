#!/usr/bin/env python3

import numpy as np

sigma_combs = #[ [0.5, 0.5], [ 0.75, 0.75 ], [ 0.1, 0.5 ], [ 0.1, 0.75 ] ]
#sigma_combs = [ [ 0.25, 0.25 ] ]
d = [ 0.1 ]
n = [ 1.001 ] 
l = [ 0.5 ] 

step = 0.01
c1_parts = list(np.arange(0, 1 + step, step)) + [ 1 + step]

# list of gamma, beta coefficients in the following order:
# beta1, gamma1, beta2, gamma2
#gamma_beta_combis = [ [1, 1, 1, 1 ] ]
gamma_beta_combis = [ [ 1, 1, 1, 1], [1, 2, 1, 2], [ 1, 2, 2, 1]]
#gamma_beta_combis = [ [ 1,4,4,1] ]

ctr = 0

exe = "./xresolution_numerical"

for n_i in n:
    for l_i in l:
        for d_i in d:
            for sigma_comb_i in sigma_combs:
                sigma12_i = sigma_comb_i[0]
                sigma21_i = sigma_comb_i[1]

                for c1min_i in c1_parts[0:-1]:
                    for gamma_beta_i in gamma_beta_combis:

                        beta1_i = gamma_beta_i[0]
                        gamma1_i = gamma_beta_i[1]
                        beta2_i = gamma_beta_i[2]
                        gamma2_i = gamma_beta_i[3]

                        p1 = sigma21_i / (sigma12_i + sigma21_i)
                        print("echo " + str(ctr))
                        ctr+=1

                        print(exe + " " + str(p1) + " " + str(n_i) + " " + str(d_i) + " " + str(sigma12_i) + " " + str(sigma21_i) + " " + str(gamma1_i) + " " + str(gamma2_i) + " " + str(beta1_i) + " " + str(beta2_i) + " " + str(c1min_i) + " " + str(c1min_i + step) + " " + str(l_i))
