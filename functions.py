import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.polynomial.polynomial as nppol
from mpl_toolkits.axes_grid1 import AxesGrid


def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 457)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 228, endpoint=False),
        np.linspace(midpoint, 1.0, 229, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap


def inv(PT1_0, PT2_0, PM, tPT1, tPT2, s1, s2, rho, c, cRI, gamma_wt, gamma_m, a, u_T1_0_1, u_T1_1_0, u_T2_0_1, u_T2_1_0,
        N, pN):
    DeltaPM = DeltaPM_inv(PT1_0, PT2_0, PM, tPT1, tPT2, s1, s2, rho, c, cRI, gamma_wt, gamma_m, a, u_T1_0_1, u_T1_1_0,
                          u_T2_0_1, u_T2_1_0, N, pN)

    if DeltaPM > 0:
        return 1
    else:
        return 0


def DeltaPM_inv(PT1_0, PT2_0, PM, tPT1, tPT2, s1, s2, rho, c, cRI, gamma_wt, gamma_m, a, u_T1_0_1, u_T1_1_0, u_T2_0_1,
                u_T2_1_0, N, pN):
    PS = get_P_eq(PT1_0, s1, rho, (1 - gamma_wt) ** a, u_T1_0_1, u_T1_1_0)
    PN = get_P_eq(PT2_0, s2, rho, (gamma_wt) ** a, u_T2_0_1, u_T2_1_0)

    GS = PS * (1 - PS)
    GN = PN * (1 - PN)
    GM = PM * (1 - PM)

    DSM = GS * GM * (PS - 1 / 2) * rho * ((1 - gamma_m) ** a - (1 - gamma_wt) ** a)
    DNM = GN * GM * (PN - 1 / 2) * rho * (gamma_m ** a - gamma_wt ** a)

    DeltaPM = DSM * (s1 + rho * (PS - 1 / 2) * (PM * (1 - gamma_m) ** a + (1 - PM) * (1 - gamma_wt) ** a)) + \
              DNM * (s2 + rho * (PN - 1 / 2) * (PM * gamma_m ** a + (1 - PM) * gamma_wt ** a)) + \
              GM * (
                      - c * rho * (N + pN) / N * (
                      GS * ((1 - gamma_m) ** a - (1 - gamma_wt) ** a) + GN * (gamma_m ** a - gamma_wt ** a)) -
                      cRI * rho * pN / N * ((1 / 2 - PS) * ((1 - gamma_m) ** a - (1 - gamma_wt) ** a) * (PS - tPT1) + (
                      1 / 2 - PN) * (gamma_m ** a - gamma_wt ** a) * (PN - tPT2))
              )
    return DeltaPM


def get_gamma_ES(PT1_0, PT2_0, PM, tPT1, tPT2, s1, s2, rho, c, cRI, a, u_T1_0_1, u_T1_1_0, u_T2_0_1, u_T2_1_0, N, pN): # IMPORTANT !! the get_gamma_ES compute the evolutionary stable values of gamma. However this function does not always work so well. Then this function does not work I get the evolutionary stable values of gamma using a invasion graph (with the file invasion_graph.py).
    zero = 0

    count = 0
    L_gamma_ES = []
    L_nature = []
    res = 1000
    L_gamma = np.linspace(0, 1, res)

    L_inv = []

    if DeltaPM_inv(PT1_0, PT2_0, PM, tPT1, tPT2, s1, s2, rho, c, cRI, L_gamma[0], L_gamma[1], a, u_T1_0_1, u_T1_1_0,
                   u_T2_0_1, u_T2_1_0, N,
                   pN) <= zero:

        count += 1
        L_nature.append('x')
        L_gamma_ES.append(0)

    if DeltaPM_inv(PT1_0, PT2_0, PM, tPT1, tPT2, s1, s2, rho, c, cRI, L_gamma[res - 1], L_gamma[res - 2], a,
                       u_T1_0_1, u_T1_1_0,
                       u_T2_0_1, u_T2_1_0, N,
                       pN) <= zero:

        count += 1
        L_nature.append('x')
        L_gamma_ES.append(1)

    for i in range(1, res - 1):
        if DeltaPM_inv(PT1_0, PT2_0, PM, tPT1, tPT2, s1, s2, rho, c, cRI, L_gamma[i], L_gamma[i+1], a, u_T1_0_1, u_T1_1_0,
                   u_T2_0_1, u_T2_1_0, N,
                   pN) > zero:
            g_p = '+'
        else:
            g_p = '-'

        if DeltaPM_inv(PT1_0, PT2_0, PM, tPT1, tPT2, s1, s2, rho, c, cRI, L_gamma[i], L_gamma[i-1], a, u_T1_0_1, u_T1_1_0,
                   u_T2_0_1, u_T2_1_0, N,
                   pN) > zero:
            g_m = '+'
        else:
            g_m ='-'

        L_inv.append([g_m, g_p])

    # is_o = False
    # following_x = 0
    # following_oxxxo = 0
    # following_gamma = []

    L_sigle = []
    for i in range(0, res - 2):
        if L_inv[i][0] == '+' and L_inv[i][1] == '-':
            L_sigle.append(('<-'))
        elif L_inv[i][0] == '-' and L_inv[i][1] == '+':
            L_sigle.append(('->'))
        elif L_inv[i][0] == '-' and L_inv[i][1] == '-':
            L_sigle.append(('x'))
        else:
            L_sigle.append(('o'))

    #print(L_sigle)

    last_arrow = ' '
    is_last_x = False
    i_last_arrow = 0
    g_0 = False

    for i in range(0, res - 2):

        if last_arrow != ' ' and last_arrow != L_sigle[i] and (L_sigle[i] == '<-' or L_sigle[i] == '->'):
            count += 1
            L_gamma_ES.append(np.mean(L_gamma[i_last_arrow:i+1]))
            if last_arrow == '->':
                L_nature.append('x')
            else:
                L_nature.append('o')

        if (i == res - 3 or i == 0) and L_sigle[i] == 'o':
            count += 1
            L_gamma_ES.append(L_gamma[i+1])
            L_nature.append('o')

        if (i == res - 3 and 1 not in L_gamma_ES) and L_sigle[i] == 'x' and last_arrow == '->':
            count += 1
            L_gamma_ES.append(L_gamma[i+1])
            L_nature.append('x')

        if (i == res - 3 and 1 in L_gamma_ES) and L_sigle[i] == 'x' and last_arrow == '<-':
            count += 1
            L_gamma_ES.append(L_gamma[i + 1])
            L_nature.append('o')

        if (last_arrow == ' ' and 0 not in L_gamma_ES) and L_sigle[i] == 'x' and L_sigle[i+1] != 'x':
            count += 1
            L_gamma_ES.append(L_gamma[i+1])
            L_nature.append('x')

        if i == res - 3 and L_sigle[i] == '<-' and 1 in L_gamma_ES:
            count += 1
            L_gamma_ES.append(L_gamma[i+1])
            L_nature.append('o')

        if i == 0 and L_sigle[i] == '->' and 0 in L_gamma_ES:
            count += 1
            L_gamma_ES.append(L_gamma[i+1])
            L_nature.append('o')

        if L_sigle[i] == '->' and is_last_x == True and last_arrow == ' ':
             count += 1
             L_gamma_ES.append(L_gamma[i + 1])
             L_nature.append('o')

        if L_sigle[i] == 'x':
            is_last_x = True
        else:
            is_last_x = False

        if L_sigle[i] == '<-':

            last_arrow = '<-'
            i_last_arrow = i

        if L_sigle[i] == '->':

            last_arrow = '->'
            i_last_arrow = i

    return count, L_gamma_ES, L_nature


# def get_gamma_ES(PT1_0, PT2_0, PM, tPT1, tPT2, s1, s2, rho, c, cRI, a, u_T1_0_1, u_T1_1_0, u_T2_0_1, u_T2_1_0, N, pN): # IMPORTANT !! the get_gamma_ES compute the evolutionary stable values of gamma. However this function does not always work so well. Then this function does not work I get the evolutionary stable values of gamma using a invasion graph (with the file invasion_graph.py).
#     zero = 0
#
#     count = 0
#     L_gamma_ES = []
#     L_nature = []
#     res = 1000
#     L_gamma = np.linspace(0, 1, res)
#
#     if DeltaPM_inv(PT1_0, PT2_0, PM, tPT1, tPT2, s1, s2, rho, c, cRI, L_gamma[0], L_gamma[1], a, u_T1_0_1, u_T1_1_0,
#                    u_T2_0_1, u_T2_1_0, N,
#                    pN) <= zero:
#         count += 1
#         L_gamma_ES.append(0)
#         L_nature.append('x')
#
#     if DeltaPM_inv(PT1_0, PT2_0, PM, tPT1, tPT2, s1, s2, rho, c, cRI, L_gamma[-1], L_gamma[-2], a, u_T1_0_1, u_T1_1_0,
#                    u_T2_0_1, u_T2_1_0, N,
#                    pN) <= zero:
#         count += 1
#         L_gamma_ES.append(1)
#         L_nature.append('x')
#
#     for i in range(1, res - 1):
#         if DeltaPM_inv(PT1_0, PT2_0, PM, tPT1, tPT2, s1, s2, rho, c, cRI, L_gamma[i], L_gamma[i + 1], a, u_T1_0_1,
#                        u_T1_1_0, u_T2_0_1, u_T2_1_0,
#                        N,
#                        pN) <= zero and DeltaPM_inv(PT1_0, PT2_0, PM, tPT1, tPT2, s1, s2, rho, c, cRI, L_gamma[i],
#                                                    L_gamma[i - 1], a,
#                                                    u_T1_0_1, u_T1_1_0, u_T2_0_1, u_T2_1_0, N,
#                                                    pN) <= zero:
#             count += 1
#             L_gamma_ES.append(L_gamma[i])
#             L_nature.append('x')
#
#         elif DeltaPM_inv(PT1_0, PT2_0, PM, tPT1, tPT2, s1, s2, rho, c, cRI, L_gamma[i], L_gamma[i + 1], a, u_T1_0_1,
#                          u_T1_1_0, u_T2_0_1, u_T2_1_0,
#                          N,
#                          pN) >= zero and DeltaPM_inv(PT1_0, PT2_0, PM, tPT1, tPT2, s1, s2, rho, c, cRI, L_gamma[i],
#                                                      L_gamma[i - 1], a,
#                                                      u_T1_0_1, u_T1_1_0, u_T2_0_1, u_T2_1_0, N,
#                                                      pN) >= zero:
#             count += 1
#             L_gamma_ES.append(L_gamma[i])
#             L_nature.append('o')
#
#     return count, L_gamma_ES, L_nature


def get_P_eq(P_0, s, rho, h, u_0_1, u_1_0):
    a0 = u_0_1
    a1 = s - rho * h / 2 - u_0_1 - u_1_0
    a2 = 3 / 2 * rho * h - s
    a3 = - rho * h

    roots = list(nppol.polyroots([a0, a1, a2, a3]))

    r = []
    for i in roots:
        if i.imag == 0.0 and i.real <= 1 and i.real >= 0:
            r.append(i.real)
    r.sort()

    if len(r) == 1:
        return r[0]
    elif len(r) == 2:
        return r[1]
    elif len(r) == 3:
        print('dependence of the initial condition')
        print('<', r[1], '=>', r[0])
        print('>', r[1], '=>', r[2])
        if P_0 < r[1]:
            return r[0]
        elif P_0 > r[1]:
            return r[2]
        else:
            return r[1]
    else:
        if u_1_0 == 0  and u_0_1 > 0:
            return 1
        elif u_1_0 > 0  and u_0_1 == 0:
            return 0
        else:
            print('Error')

#
def ind(a, b):
    if a == b:
        return 1
    else:
        return 0


def phi(T1f, T2f, T1m, T2m, rho, a, gamma):
    return (1 - (1 - ind(T1f, T1m)) * rho * (1 - gamma) ** a) * (1 - (1 - ind(T2f, T2m)) * rho * (1 - gamma) ** a)


def T(T1f, T2f, PT1, PT2, N1, N2, rho, a, gamma):
    T = (1 - PT1) * (1 - PT2) * phi(T1f, T2f, 0, 0, rho, a, gamma) \
        + (1 - PT1) * PT2 * phi(T1f, T2f, 0, 1, rho, a, gamma) \
        + PT1 * (1 - PT2) * phi(T1f, T2f, 1, 0, rho, a, gamma) \
        + PT1 * PT2 * phi(T1f, T2f, 1, 1, rho, a, gamma)
    T = T * N1 / (N1 + N2)
    return T


def TRI(T1f, T2f, tPT1, tPT2, N1, N2, rho, a, gamma, cRI):
    T = (1 - tPT1) * (1 - tPT2) * phi(T1f, T2f, 0, 0, rho, a, gamma) \
        + (1 - tPT1) * tPT2 * phi(T1f, T2f, 0, 1, rho, a, gamma) \
        + tPT1 * (1 - tPT2) * phi(T1f, T2f, 1, 0, rho, a, gamma) \
        + tPT1 * tPT2 * phi(T1f, T2f, 1, 1, rho, a, gamma)
    T = T * N2 * cRI / (N1 + N2)
    return T


def iso_rep(T1f, T2f, PT1, PT2, tPT1, tPT2, N1, N2, rho, a, gamma, cRI, c):
    t = T(T1f, T2f, PT1, PT2, N1, N2, rho, a, gamma)
    tRI = TRI(T1f, T2f, tPT1, tPT2, N1, N2, rho, a, gamma, cRI)
    return tRI / (c + (1 - c) * (t + tRI))


def mean_iso_rep(PT1_0, PT2_0, tPT1, tPT2, s1, s2, rho, c, cRI, gamma, a, u_T1_0_1, u_T1_1_0, u_T2_0_1, u_T2_1_0, N1,
                 N2):
    PT1 = get_P_eq(PT1_0, s1, rho, (1 - gamma) ** a, u_T1_0_1, u_T1_1_0)
    PT2 = get_P_eq(PT2_0, s2, rho, (gamma) ** a, u_T2_0_1, u_T2_1_0)

    m_i_r = (1 - PT1) * (1 - PT2) * iso_rep(0, 0, PT1, PT2, tPT1, tPT2, N1, N2, rho, a, gamma, cRI, c) \
            + (1 - PT1) * PT2 * iso_rep(0, 1, PT1, PT2, tPT1, tPT2, N1, N2, rho, a, gamma, cRI, c) \
            + PT1 * (1 - PT2) * iso_rep(1, 0, PT1, PT2, tPT1, tPT2, N1, N2, rho, a, gamma, cRI, c) \
            + PT1 * PT2 * iso_rep(1, 1, PT1, PT2, tPT1, tPT2, N1, N2, rho, a, gamma, cRI, c)

    return m_i_r

# def coef(i_o, j_o, k_o, l_o, i_m, j_m, k_m, l_m, i_f, j_f, k_f, l_f, rT1P1, rP1T2, rT2P2):
#     coef = 0
#
#     if i_o == i_m and j_o == j_m and k_o == k_m and l_o == l_m:
#         coef += 1/2 * (1-rT1P1) * (1-rP1T2) * (1-rT2P2)
#     if i_o == i_f and j_o == j_f and k_o == k_f and l_o == l_f:
#         coef += 1/2 * (1-rT1P1) * (1-rP1T2) * (1-rT2P2)
#
#     if i_o == i_f and j_o == j_m and k_o == k_m and l_o == l_m:
#         coef += 1/2 * (rT1P1) * (1-rP1T2) * (1-rT2P2)
#     if i_o == i_m and j_o == j_f and k_o == k_f and l_o == l_f:
#         coef += 1/2 * (rT1P1) * (1-rP1T2) * (1-rT2P2)
#
#     if i_o == i_m and j_o == j_m and k_o == k_f and l_o == l_f:
#         coef += 1/2 * (1-rT1P1) * (rP1T2) * (1-rT2P2)
#     if i_o == i_f and j_o == j_f and k_o == k_m and l_o == l_m:
#         coef += 1/2 * (1-rT1P1) * (rP1T2) * (1-rT2P2)
#
#     if i_o == i_m and j_o == j_f and k_o == k_m and l_o == l_m:
#         coef += 1/2 * (rT1P1) * (rP1T2) * (1-rT2P2)
#     if i_o == i_f and j_o == j_m and k_o == k_f and l_o == l_f:
#         coef += 1/2 * (rT1P1) * (rP1T2) * (1-rT2P2)
#
#     if i_o == i_m and j_o == j_m and k_o == k_m and l_o == l_f:
#         coef += 1 / 2 * (1 - rT1P1) * (1 - rP1T2) * (rT2P2)
#     if i_o == i_f and j_o == j_f and k_o == k_f and l_o == l_m:
#         coef += 1 / 2 * (1 - rT1P1) * (1 - rP1T2) * (rT2P2)
#
#     if i_o == i_f and j_o == j_m and k_o == k_m and l_o == l_f:
#         coef += 1 / 2 * (rT1P1) * (1 - rP1T2) * (rT2P2)
#     if i_o == i_m and j_o == j_f and k_o == k_f and l_o == l_m:
#         coef += 1 / 2 * (rT1P1) * (1 - rP1T2) * (rT2P2)
#
#     if i_o == i_m and j_o == j_m and k_o == k_f and l_o == l_m:
#         coef += 1 / 2 * (1 - rT1P1) * (rP1T2) * (rT2P2)
#     if i_o == i_f and j_o == j_f and k_o == k_m and l_o == l_f:
#         coef += 1 / 2 * (1 - rT1P1) * (rP1T2) * (rT2P2)
#
#     if i_o == i_m and j_o == j_f and k_o == k_m and l_o == l_f:
#         coef += 1 / 2 * (rT1P1) * (rP1T2) * (rT2P2)
#     if i_o == i_f and j_o == j_m and k_o == k_f and l_o == l_m:
#         coef += 1 / 2 * (rT1P1) * (rP1T2) * (rT2P2)
#
#     return coef
#
#
# def get_P_eq_PrefTrait(PT10, PT20, PP10, PP20, tPT1, tPT2, s1, s2, rho, c, cRI, gamma_wt, a, uT1_01, uT1_10, uT2_01, uT2_10, uP1_01, uP1_10, uP2_01, uP2_10, N1, N2):
#
#     PT1p = PT10
#     PP1p = PP10
#     PT2p = PT20
#     PP2p = PP20
#
#     Delta = 1
#     n = 1
#     L_Delta = []
#
#     while Delta > 10 ** (-11) and n < 10000:
#
#         n+= 1
#
#         PT1 = PT1p
#         PP1 = PP1p
#         PT2 = PT2p
#         PP2 = PP2p
#
#         GT1 = PT1 * (1-PT1)
#         GP1 = PP1 * (1 - PP1)
#         GT2 = PT2 * (1 - PT2)
#         GP2 = PP2 * (1 - PP2)
#
#         DT1P1 = rho * GT1 * GP1 * (1 - gamma_wt) ** a
#         DT2P2 = rho * GT2 * GP2 * gamma_wt ** a
#
#         PT1p = PT1 + GT1 * (s1 + rho * (1 - gamma_wt) ** a * (PP1 - 1/2)) + (1-PT1p) * uT1_01 - PT1p * uT1_10
#         PP1p = PP1 + DT1P1 * (s1 + rho * (1 - gamma_wt) ** a * (PP1 - 1/2)) \
#                + GP1 * rho * (1 - gamma_wt) ** a * (c * (N1 + N2)/N1 * (PT1 - 1/2)
#                                                     + cRI * N1 / N2 * (PT1 - tPT1)) + (1-PP1p) * uP1_01 - PP1p * uP1_10
#         PT2p = PT2 + GT2 * (s2 + rho * gamma_wt ** a * (PP2 - 1 / 2)) + (1-PT2p) * uT2_01 - PT2p * uT2_10
#         PP2p = PP2 + DT2P2 * (s2 + rho * gamma_wt ** a * (PP2 - 1 / 2)) \
#                 + GP2 * rho * (1 - gamma_wt) ** a * (c * (N1 + N2) / N1 * (PT2 - 1 / 2)
#                                                      + cRI * N1 / N2 * (PT2 - tPT2)) + (1-PP2p) * uP2_01 - PP2p * uP2_10
#
#         Delta = np.sqrt((PT1 - PT1p) ** 2 + (PP1 - PP1p) ** 2 + (PT2 - PT2p) ** 2 + (PP2 - PP2p) ** 2) / 4
#         L_Delta.append(Delta)
#
#     return PT1p, PP1p, PT2p, PP2p
