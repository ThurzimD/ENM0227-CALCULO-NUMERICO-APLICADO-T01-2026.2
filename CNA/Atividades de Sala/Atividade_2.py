# Atividade 2 
#Newton-Raphson e Secante

# Os métodos de Newton-Raphson e da Secante servem para encontrar numericamente uma raiz de uma equação, isto é, um valor \(x\) para o qual f(x) = 0

#Formulas Newton-Raphson
# x_n+1 = x_n - f(x_n)/f'(x_n)

#Formulas Secante
# x_n+1 = x_n - f(x_n) x ((x_n - x_n-1) / (f(x_n) - f(x_n-1)) )

# erro Relativo 
# Er = abs((x_n - x_n-1) / x_n) x 100%


# Bibliotecas

import pandas as pd
import math


# =========================================================
# 1. METODO DE NEWTON-RAPHSON
# =========================================================

def newton_raphson(f, df, x0, tol=1e-6, intmax=100):

    iteracoes = []
    x_valores = []
    fx_valores = []
    erro_valores = []
    erro_relativo_valores = []

    x = x0

    for inte in range(1, intmax + 1):

        fx = f(x) #Precisa definir 

        # Verifica se a raiz foi encontrada
        if abs(fx) < tol:
            iteracoes.append(inte)
            x_valores.append(x)
            fx_valores.append(fx)
            erro_valores.append(None)
            erro_relativo_valores.append(None)

            dados = pd.DataFrame({
                "Iteracao": iteracoes,
                "x": x_valores,
                "f(x)": fx_valores,
                "Erro_abs": erro_valores,
                "Erro_relativo_%": erro_relativo_valores
            })

            return x, dados, "Convergencia atingida."

        dfx = df(x)

        # Evita divisao por zero
        if abs(dfx) < 1e-14:
            dados = pd.DataFrame({
                "Iteracao": iteracoes,
                "x": x_valores,
                "f(x)": fx_valores,
                "Erro_abs": erro_valores,
                "Erro_relativo_%": erro_relativo_valores
            })

            return x, dados, "Derivada nula ou muito proxima de zero."

        # Formula de Newton-Raphson
        x_novo = x - fx / dfx

        erro = abs(x_novo - x)

        if x_novo != 0:
            erro_relativo = erro / abs(x_novo) * 100
        else:
            erro_relativo = None

        iteracoes.append(inte)
        x_valores.append(x)
        fx_valores.append(fx)
        erro_valores.append(erro)
        erro_relativo_valores.append(erro_relativo)

        # Atualiza a aproximacao
        x = x_novo

    dados = pd.DataFrame({
        "Iteracao": iteracoes,
        "x": x_valores,
        "f(x)": fx_valores,
        "Erro_abs": erro_valores,
        "Erro_relativo_%": erro_relativo_valores
    })

    return x, dados, "Maximo de iteracoes atingido."


# =========================================================
# 2. METODO DA SECANTE
# =========================================================

def secante(f, x0, x1, tol=1e-6, intmax=100):

    iteracoes = []
    x_valores = []
    fx_valores = []
    erro_valores = []
    erro_relativo_valores = []

    for inte in range(1, intmax + 1):

        f0 = f(x0)
        f1 = f(x1)

        # Verifica se uma das aproximacoes e raiz
        if abs(f1) < tol:
            iteracoes.append(inte)
            x_valores.append(x1)
            fx_valores.append(f1)
            erro_valores.append(None)
            erro_relativo_valores.append(None)
            break

        denominador = f1 - f0

        # Evita divisao por zero
        if abs(denominador) < 1e-14:
            break

        # Formula da Secante
        x_novo = x1 - f1 * (x1 - x0) / denominador

        erro = abs(x_novo - x1)

        if x_novo != 0:
            erro_relativo = erro / abs(x_novo) * 100
        else:
            erro_relativo = None

        iteracoes.append(inte)
        x_valores.append(x1)
        fx_valores.append(f1)
        erro_valores.append(erro)
        erro_relativo_valores.append(erro_relativo)

        # Atualiza os dois pontos
        x0 = x1
        x1 = x_novo

    dados = pd.DataFrame({
        "Iteracao": iteracoes,
        "x": x_valores,
        "f(x)": fx_valores,
        "Erro_abs": erro_valores,
        "Erro_relativo_%": erro_relativo_valores
    })

    if len(iteracoes) == 0:
        return x1, dados, "Nao foi possivel continuar: denominador nulo."

    if abs(f(x1)) < tol:
        mensagem = "Convergencia atingida."
    elif abs(f(x1)) < 1e-14:
        mensagem = "Convergencia atingida."
    else:
        mensagem = "Metodo interrompido ou maximo de iteracoes atingido."

    return x1, dados, mensagem
