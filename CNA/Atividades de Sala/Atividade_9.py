
import math
import pandas as pd


def razao_aurea(f, a, b, tol=1e-6, intmax=100, objetivo="min"):
    """
    f: funcao a minimizar ou maximizar.
    a, b: extremos do intervalo.
    objetivo: "min" para minimo ou "max" para maximo.
    """

    if a >= b:
        raise ValueError("O extremo a deve ser menor que b.")

    if objetivo not in ("min", "max"):
        raise ValueError("objetivo deve ser 'min' ou 'max'.")

    # Transforma maximizacao em minimizacao de -f
    def g(x):
        return f(x) if objetivo == "min" else -f(x)

    r = (math.sqrt(5) - 1) / 2

    # Dois pontos internos do intervalo
    x1 = b - r * (b - a)
    x2 = a + r * (b - a)

    f1 = g(x1)
    f2 = g(x2)

    historico = []

    for inte in range(1, intmax + 1):
        largura = b - a

        # Guarda o estado do intervalo antes da reducao
        historico.append({
            "Iteracao": inte,
            "a": a,
            "b": b,
            "x1": x1,
            "x2": x2,
            "f(x1)": f(x1),
            "f(x2)": f(x2),
            "Largura_intervalo": largura
        })

        # Criterio de parada baseado no tamanho do intervalo
        if largura < tol:
            break

        if f1 > f2:
            # O otimo esta na parte direita restante
            a = x1
            x1 = x2
            f1 = f2

            x2 = a + r * (b - a)
            f2 = g(x2)

        else:
            # O otimo esta na parte esquerda restante
            b = x2
            x2 = x1
            f2 = f1

            x1 = b - r * (b - a)
            f1 = g(x1)

    # Melhor ponto entre as duas aproximacoes finais
    x_otimo = x1 if g(x1) < g(x2) else x2

    dados = pd.DataFrame(historico)

    status = (
        "Convergencia atingida."
        if (b - a) < tol
        else "Maximo de iteracoes atingido."
    )

    return x_otimo, f(x_otimo), dados, status
