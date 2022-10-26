import os
from functions import *

eps = 0.01

PT1_0 = 0.5
PT2_0 = 0.5
PM = 0.01
tPT1 = 1
tPT2 = 0.4
rho = eps
c = 0 * eps
cRI = 0.5 * eps
a = 1
N = 10
pN = 10
s1 = 2 * eps
s2 = 0

u_T1_0_1 = 0.3 * eps
u_T1_1_0 = 0.1 * eps
u_T2_0_1 = 0.2 * eps
u_T2_1_0 = 0.1 * eps

res = 250
L_log_a = np.linspace(-1, 1, res)
L_c = np.linspace(0, 0.5 * eps, res)

R = np.zeros((res, res))

for i, log_a in enumerate(L_log_a):

    print((i / res * 100) // 1, '%')
    a = np.exp(log_a)

    for j, c in enumerate(L_c):

        if i > 25 * 2.5 and os.path.isfile('results/c/' + str(i) + '_' + str(j) + '.txt') == False:
            count, L_gamma_ES, L_nature = get_gamma_ES(PT1_0, PT2_0, PM, tPT1, tPT2, s1, s2, rho, c, cRI, a, u_T1_0_1,
                                                       u_T1_1_0, u_T2_0_1, u_T2_1_0,
                                                       N, pN)
            text_results = str(count) + '\n' + str(L_gamma_ES) + '\n' + str(L_nature)

            file = open('results/c/' + str(i) + '_' + str(j) + '.txt', 'a')
            file.write(text_results)
            file.close()
