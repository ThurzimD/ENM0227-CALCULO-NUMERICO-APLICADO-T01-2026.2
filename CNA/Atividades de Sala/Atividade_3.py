
import numpy as np
import pandas as pd


def newton_sistemas(F, J, x0, tol=1e-6, intmax=100):
    """
    F: funcao que retorna o vetor das equacoes.
    J: funcao que retorna a matriz Jacobiana.
    x0: vetor com os chutes iniciais.
    """

    x = np.array(x0, dtype=float)

    # Listas para armazenar os resultados
    historico = []
    norma_residuo = []
    erros = []

    for inte in range(1, intmax + 1):
        fx = np.asarray(F(x), dtype=float)

        # Norma do vetor de residuos
        residuo = np.linalg.norm(fx)

        historico.append(x.copy())
        norma_residuo.append(residuo)

        # Criterio de parada
        if residuo < tol:
            break

        jacobiana = np.asarray(J(x), dtype=float)

        try:
            # Resolve J * delta = -F.
            # Evita calcular a inversa da Jacobiana.
            delta = np.linalg.solve(jacobiana, -fx)
        except np.linalg.LinAlgError:
            raise ValueError(
                "Jacobiana singular: nao foi possivel continuar."
            )

        erro = np.linalg.norm(delta)
        erros.append(erro)

        # Atualiza todas as variaveis simultaneamente
        x = x + delta

    # Organiza os resultados em uma tabela
    dados = pd.DataFrame({
        "Iteracao": range(1, len(historico) + 1),
        "Aproximacao": [v.tolist() for v in historico],
        "Norma_F": norma_residuo
    })

    if np.linalg.norm(F(x)) < tol:
        status = "Convergencia atingida."
    else:
        status = "Maximo de iteracoes atingido."

    return x, dados, status
