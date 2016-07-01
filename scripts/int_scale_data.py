import numpy as np

box = [150, 150, 250, 250, 250, 600]
integrals = np.zeros(len(box))
L = np.zeros(len(box))

run1_int = [.4035, .3632, .5301, .4605]
integrals[0] = np.average(run1_int)
run2_int = [.42335, .4177, .5192]
integrals[1] = np.average(run2_int)
run3_int = [.389, .32176]
integrals[2] = np.average(run3_int)
run4_int = [.377, .426]
integrals[3] = np.average(run4_int)
run5_int = [.3799, .48369]
integrals[4] = np.average(run5_int)
run7_int = [.4021, .376368]
integrals[5] = np.average(run7_int)

driving_scale_strat = [33.2, 36.591, 46.63, 50, 55.147, 108.031]
driving_scale_iso = [15.981, 17.552, 22.424, 23.991, 26.374, 51.688]

i = 0

while (i < len(box)):
        L[i] = box[i]*integrals[i]
        i = i+1

ratio_strat = L[:]/driving_scale_strat[:]
ratio_iso = L[:]/driving_scale_iso[:]

print(L)
print('---')
print(ratio_iso)
print(ratio_strat)
f_iso = np.average(ratio_iso)
f_strat = np.average(ratio_strat)
print('f_iso = {:}'.format(f_iso))
print('f_strat = {:}'.format(f_strat))
