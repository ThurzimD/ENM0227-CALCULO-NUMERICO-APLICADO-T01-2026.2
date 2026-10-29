
import numpy as np
import pandas as pd


def gauss_ingenua(A, b, tol=1e-14):
    # Cria copias para nao modificar os dados originais
    A = np.array(A, dtype=float, copy=True)
    b = np.array(b, dtype=float, copy=True)

    n = len(b)

    if A.shape != (n, n):
        raise ValueError("A matriz A deve ser quadrada.")

    if not np.all(np.isfinite(A)) or not np.all(np.isfinite(b)):
        raise ValueError("A matriz e o vetor devem ter valores finitos.")

    historico = []

    # Eliminacao progressiva
    for k in range(n - 1):

        # Pivo da etapa atual
        pivo = A[k, k]

        if abs(pivo) < tol:
            raise ValueError(
                "Pivo nulo ou muito pequeno. Use pivoteamento parcial."
            )

        # Zera os elementos abaixo do pivo
        for i in range(k + 1, n):
            multiplicador = A[i, k] / pivo

            A[i, k:] -= multiplicador * A[k, k:]
            b[i] -= multiplicador * b[k]

        historico.append({
            "Etapa": k + 1,
            "Matriz": A.copy(),
            "Vetor_b": b.copy()
        })

    if abs(A[n - 1, n - 1]) < tol:
        raise ValueError("Sistema singular ou pivo final muito pequeno.")

    # Substituicao retroativa
    x = np.zeros(n)

    for i in range(n - 1, -1, -1):
        soma = np.dot(A[i, i + 1:], x[i + 1:])
        x[i] = (b[i] - soma) / A[i, i]

    dados = pd.DataFrame([
        {
            "Etapa": item["Etapa"],
            "Matriz": item["Matriz"].tolist(),
            "Vetor_b": item["Vetor_b"].tolist()
        }
        for item in historico
    ])

    return x, dados
