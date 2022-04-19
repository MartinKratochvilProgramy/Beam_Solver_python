#function that numericaly solves initial beam
def beam_solver(point_loads, point_moments, linear_loads, L, Xa, Xb, J):
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches


    E = 2.1 * 10 ** 11  # Youngův modul PRO ŽELEZO
    # J = (0.2 * 0.5 ** 3) / 12  -zadává se v okně
    EJ = E * J


    N = 10000                                       # počet dělení po délce nosníku
    dx = L / N                                      # délka kroku
    X = np.arange(dx, L + dx, dx)                   # vektor pozic po délce nosníku
    fi = np.zeros(len(X))                           # vektor natočení
    v = np.zeros(len(X))                            # vektor průhybů

    reaction_forces = np.array([0.0, 0, 0])

    Shear = np.zeros(len(X))
    Moment = np.zeros(len(X))

    support_index_A = int(Xa / dx)
    support_index_B = int(Xb / dx) - 1
    dFi = 0.00001                                   #krok natočení pro numerické řešení průhybů

    def reactions_PL(n):
        xp = point_loads[n, 0]                       # pozice síly
        fx = point_loads[n, 1]                       # horizontální složka
        fy = point_loads[n, 2]                       # vertikální složka

        la_p = Xa - xp                              # rameno vůči A
        mp = fy * la_p                              # moment vůči A
        la_Rb = Xb - Xa                             # rameno B vůči A

        Rb = mp / la_Rb                             # reakce B
        Ra = -fy - Rb                               # reakce A
        Ha = -fx                                    # reakce v horizontálním směru

        return Ra, Rb, Ha

    def shear_moment_PL(n, Shear, Moment):
        xp = point_loads[n, 0]                       # pozice síly
        fy = point_loads[n, 2]                       # síla do osy y
        Ra = PL_record[n, 0]                        # reakce v A od jednotlivých sil
        Rb = PL_record[n, 2]                        # reakce v B -||-

        # projede nosník a spočte posouvací sílu a ohybový moment od síly v každém bodě

        for i, x in enumerate(X):
            shear = 0
            moment = 0
            if x > Xa:
                shear += Ra
                moment -= Ra * (x - Xa)

            if x > Xb:
                shear += Rb
                moment -= Rb * (x - Xb)

            if x > xp:
                shear += fy
                moment -= fy * (x - xp)

            Shear[i] += shear
            Moment[i] += moment

        return Shear, Moment

    def reactions_PM(n):
        m = point_moments[n, 1]
        la_vb = Xb - Xa

        Rb = m / la_vb
        Ra = -Rb

        return Ra, Rb

    def shear_moment_PM(n, Shear, Moment):
        xm = point_moments[n, 0]
        m = point_moments[n, 1]
        la_vb = Xb - Xa

        Rb = m / la_vb
        Ra = -Rb

        # projede nosník a spočte posouvací sílu a ohybový moment od momentu v každém bodě

        for i, x in enumerate(X):
            shear = 0
            moment = 0
            if x > Xa:
                shear += Ra
                moment -= Ra * (x - Xa)

            if x > Xb:
                shear += Rb
                moment -= Rb * (x - Xb)

            if x > xm:
                moment -= m

            Shear[i] += shear
            Moment[i] += moment

        return Shear, Moment

    def reactions_DL(n):
        xStart = linear_loads[n, 0]
        xEnd = linear_loads[n, 1]
        fy = linear_loads[n, 2]

        fy_Res = fy * (xEnd - xStart)
        x_res = xStart + 0.5 * (xEnd - xStart)

        la_p = Xa - x_res  # rameno vůči A
        mp = fy_Res * la_p  # moment vůči A
        la_Rb = Xb - Xa  # rameno B vůči A

        Rb = mp / la_Rb  # reakce B
        Ra = -fy_Res - Rb  # reakce A

        return Ra, Rb

    def shear_moment_DL(n, Shear, Moment):
        xStart = linear_loads[n, 0]
        xEnd = linear_loads[n, 1]
        fy = linear_loads[n, 2]
        Ra = DL_record[n, 0]
        Rb = DL_record[n, 1]

        for i, x in enumerate(X):
            shear = 0
            moment = 0
            if x > Xa:
                shear += Ra
                moment -= Ra * (x - Xa)

            if x > Xb:
                shear += Rb
                moment -= Rb * (x - Xb)

            if x > xStart and x <= xEnd:
                shear += fy * (x - xStart)
                moment -= fy * (x - xStart) * 0.5 * (x - xStart)
            elif (x > xEnd):
                shear += fy * (xEnd - xStart)
                moment -= fy * (xEnd - xStart) * (x - xStart - 0.5 * (xEnd - xStart))

            Shear[i] += shear
            Moment[i] += moment

        return Shear, Moment

    def plot_out():
        fig, ax = plt.subplots(4, sharex=True, figsize=(12, 8))

        if len(point_loads) != 0:
            F_max = max(abs(point_loads[:, 2]))
        else:
            F_max = 0

        if len(point_moments) != 0:
            M_max = max(abs(point_moments[:, 1]))
        else:
            M_max = 0

        if len(linear_loads) != 0:
            q_max = max(abs(linear_loads[:, 2]))
        else:
            q_max = 0

        y_max = max(F_max, M_max, q_max) * 1.2
        y_min = -max(F_max, M_max, q_max) * 1.2
        y_abs = max(abs(y_max), abs(y_min))

        # VYKRESLIT SÍLY
        for n, f in enumerate(point_loads):
            if point_loads[n, 2] < 0:
                ax[0].arrow(point_loads[n, 0], -point_loads[n, 2], 0, -0.2*point_loads[n, 2] + point_loads[n, 2], head_width=L*0.05,
                            head_length=-0.2*point_loads[n, 2])
            else:
                ax[0].arrow(point_loads[n, 0], -point_loads[n, 2], 0, -0.2*y_abs + point_loads[n, 2], head_width=L*0.05,
                            head_length=0.2*y_abs)
        # VYKRESLIT MOMENTY
        for n, f in enumerate(point_moments):
            ax[0].plot([point_moments[n, 0], point_moments[n, 0]], [point_moments[n, 1], -point_moments[n, 1]], color='red')
            ax[0].plot([point_moments[n, 0], point_moments[n, 0] + abs(point_moments[n, 1]) * 0.005], [point_moments[n, 1], point_moments[n, 1]],
                       color='red')
            ax[0].plot([point_moments[n, 0], point_moments[n, 0] - abs(point_moments[n, 1]) * 0.005], [-point_moments[n, 1], -point_moments[n, 1]],
                       color='red')
        # VYKRESLIT SPOJITÉ ZATÍŽENÍ
        for n, f in enumerate(linear_loads):
            rect = mpatches.Rectangle((linear_loads[n, 0], 0), linear_loads[n, 1] - linear_loads[n, 0], -linear_loads[n, 2],
                                      # fill=False,
                                      alpha=0.1,
                                      facecolor="blue")
            ax[0].add_patch(rect)

        ax[0].set_ylim(y_min, y_max)
        ax[0].set_xlim([-L*0.1, L*1.1])
        ax[0].arrow(Xa, -0.4*y_abs, 0, 0.2*y_abs, head_width=0.05*L, head_length=0.2*y_abs, color='red')
        ax[0].arrow(Xb, -0.4*y_abs, 0, 0.2*y_abs, head_width=0.05*L, head_length=0.2*y_abs, color='red')
        ax[0].plot([L, 0], [0, 0], color='black', linewidth='2')


        ax[1].set_title('Posouvací síla [kN]')
        ax[1].plot(X, Shear, color='blue')

        ax[2].set_title('Ohybový moment [kNm]')
        ax[2].plot(X, Moment, color='red')


        ax[3].plot(X, -v)
        ax[3].set_title('Průhyb [mm]')

        plt.show()

    def solve_displacements(Moment, EJ, fi_guess):

        # INICIALIZACE ODHADU ŘEŠENÍ
        fi[support_index_A] = fi_guess

        # PROJEDE NOSNÍK OD PODPORY A A SPOČTE NATOČENÍ A PRŮHYB
        for i, m in enumerate(Moment[support_index_A::]):
            index = i + support_index_A
            if index < len(Moment) - 1:
                fi[index + 1] = fi[index] + dx / (2 * EJ) * (Moment[index + 1] + Moment[index]) * 1000
                v[index + 1] = v[index] + dx / 2 * (fi[index + 1] + fi[index])

        err = 0 - v[support_index_B]

        #POKUD PŘED PODPOROU A JE ČÁST NOSNÍKU, PROJEDE TUTO ČÁST A SPOČTE PRŮHYBY
        if Xa != 0 and abs(err) < 0.001:
            fi[support_index_A-1] = fi[support_index_A]
            for i in range(support_index_A - 1):
                index = support_index_A - i - 1
                fi[index - 1] = fi[index] - dx / (2 * EJ) * (Moment[index] + Moment[index - 1]) * 1000
                v[index - 1] = v[index] - dx / 2 * (fi[index] + fi[index - 1])

        return fi, v, err

    # PROJEDE VŠECHNA ZATÍŽENÍ A SPOČTE REAKCE V PODPORÁCH OD SÍLY
    PL_record = np.empty([0, 3])
    if (len(point_loads) > 0):
        for n, p in enumerate(point_loads):
            Ra, Rb, Ha = reactions_PL(n)
            PL_record = np.append(PL_record, [np.array([Ra, Ha, Rb])], axis=0)

            reaction_forces[0] += Ra
            reaction_forces[1] += Ha
            reaction_forces[2] += Rb

    # PROJEDE OHYBOVÉ MOMENTY A SPOČTE REAKCE V PODPORÁCH OD MOMENTU
    PM_record = np.empty([0, 2])
    if (len(point_moments) > 0):
        for n, p in enumerate(point_moments):
            Ra, Rb = reactions_PM(n)
            PM_record = np.append(PM_record, [np.array([Ra, Rb])], axis=0)

            reaction_forces[0] += Ra
            reaction_forces[2] += Rb

    # PROJEDE OHYBOVÉ MOMENTY A SPOČTE REAKCE V PODPORÁCH OD SPOJITÉHO ZATÍŽENÍ
    DL_record = np.empty([0, 2])
    if (len(linear_loads) > 0):
        for n, p in enumerate(linear_loads):
            Ra, Rb = reactions_DL(n)
            DL_record = np.append(DL_record, [np.array([Ra, Rb])], axis=0)

            reaction_forces[0] += Ra
            reaction_forces[2] += Rb

    # SPOČTE POSOUVACÍ SÍLU A MOMENT PO DÉLCE NOSNÍKU OD SÍLY
    if (len(point_loads) > 0):
        for n, p in enumerate(point_loads):
            Shear, Moment = shear_moment_PL(n, Shear, Moment)

    # SPOČTE POSOUVACÍ SÍLU A MOMENT PO DÉLCE NOSNÍKU OD MOMENTU
    if (len(point_moments) > 0):
        for n, p in enumerate(point_moments):
            Shear, Moment = shear_moment_PM(n, Shear, Moment)

    # SPOČTE POSOUVACÍ SÍLU A MOMENT PO DÉLCE NOSNÍKU OD SPOJITÉHO ZATÍŽENÍ
    if (len(linear_loads) > 0):
        for n, p in enumerate(linear_loads):
            Shear, Moment = shear_moment_DL(n, Shear, Moment)

    fi_guess = 0.001    #odhad počátečního natočení

    # SPOČTE NATOČENÍ NUMERICKY
    print("Solving...")
    while True:
        fi, v, err = solve_displacements(Moment, EJ, fi_guess)
        if err < 0:
            fi_guess -= dFi
        elif err > 0:
            fi_guess += dFi

        if abs(err) < 0.0001:
            break

    print(
        f"Ra = {round(reaction_forces[0], 2)} kN \t Rb = {round(reaction_forces[2], 2)} kN")
    print(f"Mo_max = {round(max(Moment), 3)} kNm, Mo_min = {round(min(Moment), 3)} kNm")

    plot_out()


