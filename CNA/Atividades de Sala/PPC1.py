# Ruger Kutta 4º Ordem


import pandas as pd


def runge_kutta_4(f, t0, y0, tf, h):
    """
    Resolve dy/dt = f(t, y) pelo método RK4.

    Parâmetros:
        f  : função da EDO
        t0 : tempo inicial
        y0 : condição inicial
        tf : tempo final
        h  : tamanho do passo

    Retorna:
        tabela com os valores calculados.
    """

    t = t0
    y = y0

    historico = []

    # Registra a condição inicial
    historico.append({
        "t": t,
        "y": y
    })

    # Avança enquanto não atingir o tempo final
    while t < tf:
        # Evita ultrapassar o tempo final
        passo = min(h, tf - t)

        # Quatro estimativas de inclinação
        k1 = f(t, y)
        k2 = f(t + passo / 2, y + passo * k1 / 2)
        k3 = f(t + passo / 2, y + passo * k2 / 2)
        k4 = f(t + passo, y + passo * k3)

        # Atualização da solução
        y = y + (passo / 6) * (
            k1 + 2*k2 + 2*k3 + k4
        )

        # Atualiza a variável independente
        t = t + passo

        # Armazena o resultado do passo
        historico.append({
            "t": t,
            "y": y
        })

    return pd.DataFrame(historico)


# Exemplo: dy/dt = -2y, com y(0) = 1
def f(t, y):
    return -2*y


dados = runge_kutta_4(
    f=f,
    t0=0,
    y0=1,
    tf=2,
    h=0.1
)

print(dados)

# Salva os resultados
dados.to_csv("resultados_RK4.csv", index=False)
