# Método de Baristow 


import numpy as np
import pandas as pd


def bairstow(coeficientes, r=0.1, s=0.1,
             tol=1e-8, intmax=100):
    """
    Encontra raízes de um polinômio pelo método de Bairstow.

    coeficientes:
        Lista dos coeficientes em ordem decrescente.
        Exemplo: x^3 - 6x^2 + 11x - 6
        -> [1, -6, 11, -6]

    r e s:
        Estimativas iniciais do fator x² - r*x - s.

    tol:
        Tolerância numérica.

    intmax:
        Máximo de iterações por fator quadrático.

    Retorna:
        raízes, histórico e status.
    """

    a = np.array(coeficientes, dtype=float)

    # Remove coeficientes iniciais iguais a zero
    a = np.trim_zeros(a, trim="f")

    if len(a) < 2:
        raise ValueError("Informe um polinômio não constante.")

    if a[0] == 0:
        raise ValueError("O primeiro coeficiente deve ser diferente de zero.")

    # Normaliza o polinômio para facilitar os cálculos
    a = a / a[0]

    raizes = []
    historico = []

    # Divide o polinômio até restarem grau 1 ou 2
    while len(a) > 3:

        convergiu = False

        for it in range(1, intmax + 1):

            n = len(a) - 1

            # Divisão sintética pelo fator x² - r*x - s
            b = np.zeros(n + 1)
            c = np.zeros(n + 1)

            b[0] = a[0]
            b[1] = a[1] + r*b[0]

            for i in range(2, n + 1):
                b[i] = a[i] + r*b[i-1] + s*b[i-2]

            # Recorrência auxiliar para calcular as correções
            c[0] = b[0]
            c[1] = b[1] + r*c[0]

            for i in range(2, n + 1):
                c[i] = b[i] + r*c[i-1] + s*c[i-2]

            # Resíduo da divisão quadrática
            R = b[n-1]
            S = b[n]

            # Sistema linear para corrigir r e s
            J = np.array([
                [c[n-2], c[n-3]],
                [c[n-1], c[n-2]]
            ], dtype=float)

            residuo = np.array([R, S], dtype=float)

            try:
                dr, ds = np.linalg.solve(J, -residuo)
            except np.linalg.LinAlgError:
                return raizes, pd.DataFrame(historico), (
                    "Falha: sistema singular. Tente outros valores iniciais."
                )

            r += dr
            s += ds

            historico.append({
                "grau_atual": n,
                "iteracao": it,
                "r": r,
                "s": s,
                "residuo_R": R,
                "residuo_S": S,
                "correcao_r": dr,
                "correcao_s": ds
            })

            if max(abs(dr), abs(ds)) < tol:
                convergiu = True
                break

        if not convergiu:
            return raizes, pd.DataFrame(historico), (
                "Falha de convergência. Tente outros valores iniciais."
            )

        # Raízes do fator x² - r*x - s = 0
        delta = r**2 + 4*s

        raiz1 = (r + np.lib.scimath.sqrt(delta)) / 2
        raiz2 = (r - np.lib.scimath.sqrt(delta)) / 2

        raizes.extend([raiz1, raiz2])

        # O quociente da divisão é b[0:n-1]
        a = b[:n-1].copy()

        # Reinicia as estimativas para o próximo fator
        r, s = 0.1, 0.1

    # Resolve diretamente o polinômio restante
    raizes_restantes = np.roots(a)
    raizes.extend(raizes_restantes)

    return raizes, pd.DataFrame(historico), "Concluído."


# Exemplo: x^4 - 5x^2 + 4 = 0
coef = [1, 0, -5, 0, 4]

raizes, dados, status = bairstow(coef)

print("Status:", status)
print("Raízes:", raizes)
print(dados)

dados.to_csv("resultados_Bairstow.csv", index=False)
