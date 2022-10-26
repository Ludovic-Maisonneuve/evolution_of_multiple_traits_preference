import os
from functions import *

eps = 0.01

PT1_0 = 0.5
PT2_0 = 0.5
PM = 0.01
tPT1 = 0
tPT2 = 0
rho = eps
c = 0.1 * eps
cRI = 0.5 * eps
N = 10
pN = 10
s1 = 0 * eps
s2 = 0 * eps
a = np.exp(-1)

u_T1_0_1 = 0.3 * eps
u_T1_1_0 = 0.1 * eps
u_T2_0_1 = 0.2 * eps
u_T2_1_0 = 0.1 * eps

res = 250
L_tP = np.linspace(0, 1, res)

R = np.zeros((res, res))

for i, tPT1 in enumerate(L_tP):

    print((i / res * 100) // 1, '%')

    for j, tPT2 in enumerate(L_tP):

        if i > 79 * 2.5 and os.path.isfile('results/a/' + str(i) + '_' + str(j) + '.txt') == False:
            count, L_gamma_ES, L_nature = get_gamma_ES(PT1_0, PT2_0, PM, tPT1, tPT2, s1, s2, rho, c, cRI, a, u_T1_0_1,
                                                       u_T1_1_0, u_T2_0_1, u_T2_1_0,
                                                       N, pN)
            text_results = str(count) + '\n' + str(L_gamma_ES) + '\n' + str(L_nature)

            file = open('results/a/' + str(i) + '_' + str(j) + '.txt', 'a')
            file.write(text_results)
            file.close()
