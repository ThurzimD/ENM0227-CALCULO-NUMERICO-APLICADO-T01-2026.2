
import numpy as np
import pandas as pd


def gauss_seidel(A, b, x0=None, tol=1e-6, intmax=100):
    A = np.array(A, dtype=float, copy=True)
    b = np.array(b, dtype=float, copy=True)

    n = len(b)

    if A.shape != (n, n):
        raise ValueError("A matriz A deve ser quadrada.")

    if not np.all(np.isfinite(A)) or not np.all(np.isfinite(b)):
        raise ValueError("A matriz e o vetor devem ter valores finitos.")

    # Chute inicial: vetor de zeros, se nao for fornecido
    if x0 is None:
        x = np.zeros(n)
    else:
        x = np.array(x0, dtype=float, copy=True)

    if x.shape != (n,):
        raise ValueError("O chute inicial deve ter n elementos.")

    if np.any(np.abs(np.diag(A)) < 1e-14):
        raise ValueError("A matriz possui elemento diagonal nulo.")

    historico = []
    erros = []
    residuos = []

    for inte in range(1, intmax + 1):
        x_anterior = x.copy()

        # Atualiza as incognitas uma por uma
        for i in range(n):
            soma_antes = np.dot(A[i, :i], x[:i])
            soma_depois = np.dot(A[i, i + 1:], x_anterior[i + 1:])

            x[i] = (
                b[i] - soma_antes - soma_depois
            ) / A[i, i]

        # Erro entre duas aproximacoes sucessivas
        erro = np.linalg.norm(x - x_anterior)

        # Residuo: mede quanto A*x difere de b
        residuo = np.linalg.norm(A @ x - b)

        historico.append(x.copy())
        erros.append(erro)
        residuos.append(residuo)

        if erro < tol and residuo < tol:
            break

    dados = pd.DataFrame({
        "Iteracao": range(1, len(historico) + 1),
        "Aproximacao": [v.tolist() for v in historico],
        "Erro": erros,
        "Residuo": residuos
    })

    if erros[-1] < tol and residuos[-1] < tol:
        status = "Convergencia atingida."
    else:
        status = "Maximo de iteracoes atingido."

    return x, dados, status
