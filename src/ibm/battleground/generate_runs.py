#!/usr/bin/env python3

import numpy as np

exe_name = "./xpo_conflict"

mu_s = 0.01
mu_q = 0.01
sdmu = 0.01

step = 0.025
c1list = [ 0.83 ]
c2list = list(np.arange(0,1+step,step))

llist = [ 0.5 ]

d = 0.1
sigma12 = [ 0.25 ]
sigma21 = [ 0.25 ]

# list of gamma, beta coefficients in the following order:
# beta1, gamma1, beta2, gamma2
gamma_beta_combis = [[ 1, 2, 1, 2 ], [ 1, 1, 1, 1 ], [1, 2, 2, 1 ] ]
initqS = 0.9
initqNS = 0.1
inits1 = 0.9
inits2 = 0.1

ctr = 0 


for sigma12_i in sigma12:
    for sigma21_i in sigma21:
        for c1_i in c1list:
            for c2_i in c2list:
                for li in llist:
                    for gamma_beta_i in gamma_beta_combis:

                        beta1_i = gamma_beta_i[0]
                        gamma1_i = gamma_beta_i[1]
                        beta2_i = gamma_beta_i[2]
                        gamma2_i = gamma_beta_i[3]

                        for control_i in [ 0, 1]:

                            p = sigma21_i / (sigma12_i + sigma21_i)
                            print("echo " + str(ctr))
                            ctr+=1
                            print(exe_name + " " + str(mu_s) + " "  + str(sdmu) + " " + str(c1_i) + " " + str(c2_i) + " " + \
                                    str(p) + " " + str(d) + " " + str(li) + " " + str(sigma12_i) + " " + str(sigma21_i) + " " + str(gamma1_i) + " " + str(gamma2_i) + " " + str(beta1_i) + " " + str(beta2_i) + " "  + str(control_i))
