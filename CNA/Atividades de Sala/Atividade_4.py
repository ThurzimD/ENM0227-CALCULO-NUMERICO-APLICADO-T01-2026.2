
import numpy as np
import pandas as pd


def muller(f, x0, x1, x2, tol=1e-6, intmax=100):
    # Müller pode trabalhar com aproximacoes complexas
    x0, x1, x2 = complex(x0), complex(x1), complex(x2)

    historico = []
    fx_historico = []
    erros = []

    for inte in range(1, intmax + 1):
        f0, f1, f2 = f(x0), f(x1), f(x2)

        # Distancias entre os tres pontos
        h0 = x1 - x0
        h1 = x2 - x1

        if abs(h0) < 1e-14 or abs(h1) < 1e-14:
            raise ValueError("Aproximacoes iniciais muito proximas.")

        # Diferencas divididas
        d0 = (f1 - f0) / h0
        d1 = (f2 - f1) / h1

        # Coeficientes da parabola interpoladora
        a = (d1 - d0) / (h1 + h0)
        b = a * h1 + d1
        c = f2

        discriminante = b**2 - 4*a*c
        raiz_disc = np.sqrt(discriminante)

        # Escolhe o denominador de maior modulo
        den1 = b + raiz_disc
        den2 = b - raiz_disc
        denominador = den1 if abs(den1) >= abs(den2) else den2

        if abs(denominador) < 1e-14:
            raise ValueError("Denominador muito proximo de zero.")

        # Nova aproximacao
        x3 = x2 - (2*c) / denominador

        erro = abs(x3 - x2)

        historico.append(x3)
        fx_historico.append(f(x3))
        erros.append(erro)

        # Atualiza as tres aproximacoes
        x0, x1, x2 = x1, x2, x3

        if abs(f(x3)) < tol:
            break

    dados = pd.DataFrame({
        "Iteracao": range(1, len(historico) + 1),
        "x": historico,
        "f(x)": fx_historico,
        "Erro_abs": erros
    })

    status = (
        "Convergencia atingida."
        if abs(f(x2)) < tol
        else "Maximo de iteracoes atingido."
    )

    return x2, dados, status
