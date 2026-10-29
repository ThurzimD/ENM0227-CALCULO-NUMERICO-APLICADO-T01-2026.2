
import numpy as np
import pandas as pd


def inversa_lu(A, tol=1e-14):
    A = np.array(A, dtype=float, copy=True)

    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("A matriz deve ser quadrada.")

    if not np.all(np.isfinite(A)):
        raise ValueError("A matriz deve ter valores finitos.")

    n = A.shape[0]

    # Inicializa L como identidade e U como copia de A
    L = np.eye(n)
    U = A.copy()

    # P representa as permutacoes de linhas
    P = np.eye(n)

    # Decomposicao PA = LU com pivoteamento parcial
    for k in range(n):

        linha_pivo = k + np.argmax(np.abs(U[k:, k]))

        if abs(U[linha_pivo, k]) < tol:
            raise ValueError("Matriz singular: inversa inexistente.")

        if linha_pivo != k:
            U[[k, linha_pivo]] = U[[linha_pivo, k]]
            P[[k, linha_pivo]] = P[[linha_pivo, k]]

            # Troca somente a parte ja calculada de L
            if k > 0:
                L[[k, linha_pivo], :k] = L[[linha_pivo, k], :k]

        for i in range(k + 1, n):
            L[i, k] = U[i, k] / U[k, k]
            U[i, k:] -= L[i, k] * U[k, k:]

    # Resolve A*x = e_j para cada coluna da identidade
    identidade = np.eye(n)
    A_inv = np.zeros((n, n))
    historico = []

    for j in range(n):
        # PA*x = P*e_j, ou seja, L*U*x = P*e_j
        rhs = P @ identidade[:, j]

        # Substituicao progressiva: L*y = P*e_j
        y = np.zeros(n)

        for i in range(n):
            y[i] = rhs[i] - np.dot(L[i, :i], y[:i])

        # Substituicao retroativa: U*x = y
        x = np.zeros(n)

        for i in range(n - 1, -1, -1):
            soma = np.dot(U[i, i + 1:], x[i + 1:])
            x[i] = (y[i] - soma) / U[i, i]

        A_inv[:, j] = x

        historico.append({
            "Coluna_identidade": j + 1,
            "Solucao": x.copy()
        })

    dados = pd.DataFrame({
        "Coluna_identidade": [
            item["Coluna_identidade"] for item in historico
        ],
        "Solucao": [
            item["Solucao"].tolist() for item in historico
        ]
    })

    return A_inv, dados, L, U, P
