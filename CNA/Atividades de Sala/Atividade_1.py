### Atividade de Sala 01 ###

#Atividade Prática de Programação Métodos da Bisseção e da Falsa Posição


# Aluno: Arthur D. / 222005976
# Atividade 1: Escreva um programa de computador que calcule, utilizando o método da Bissecção.
# a raiz da seguinte equação algébrica: 
# 
# f(H) = ((0,471405*(20*H)^(5/3)) / (20+2H)^(2/3))-5 = 0 
#-------------------------------------------------------------
# Resposta esperada: H = 0,7023 (m)
# Faça um plot: Ed*Niter
# Ed = ((xr-xv)/xv) x 100
#------------------------------------------------------------

#------------ Bissecção-----
# x0+x1/2 
#---------------------------

#------------Falsa Posição---------
# xr = (x0 x f(x1) - x1 x f(x0)) / (f(x1) - f(x0))


import pandas as pd 
import math as math
import csv as csv


def function(H):
    fun=0.471405 * ((20 * H)**(5/3) / (20 + (2 * H))**(2/3)) - 5
    return fun

# Definir um intervalo Inicial

x0 = 0.5 #Hl
x1 = 1.0 #Hu

def verificacao (x0, x1):
    verf = function(x0)*function(x1)
    print("x0 =", x0)
    print("x1 =", x1)
    print("Produto de f(x0) x f(x1)=", verf)
    if verf < 0:
       a = print(f"Resultado é {verf}, Intervalo valido")
       return True
    else:
       b = print(f"Resultado é {verf}, Intervalo não valido")
    return False


def bissecao(x0, x1):
    intmax = 100
    tol = 1e-6
    xr = x0

    # condições 

    if function(x0) ==0:
        return x0, 0, "Raiz exata no extremo x0."
    
    if function(x1) == 0:
        return x1, 0, "Raiz exata no extremo x1."

    if not verificacao(x0, x1):
        return None, 0, "Intervalo inicial inválido."

    #Lista para Armazenamento 
    xr_valores_B = []
    interacao_B = []
    fxr_valores_B = []

    #Laço
    
    for inte in range(1, intmax + 1):
        xr = (x0 + x1) / 2 

        #Armazenando
        xr_valores_B.append(xr)
        interacao_B.append(inte)
        fxr_valores_B.append(function(xr))


        if function(xr) ==0:
            return xr, inte, "Raiz exata."

        if abs(function(xr)) < tol:
            return xr, intmax, xr_valores_B, interacao_B, fxr_valores_B
        
        if function(x0) * function(xr) < 0:
            x1 = xr
        else:
            x0 = xr

    return xr, intmax, xr_valores_B, interacao_B, fxr_valores_B

raiz, niter, xr_lista, iter_lista, fxr_lista = bissecao(x0, x1)

print("BISSECÇÃO:")
print(round(raiz, 4))


# Método da Falsa Posição
def falsa_posicao(x0, x1):
    intmax = 150
    tol = 1e-6

    if function(x0) == 0:
        return x0, 0, "Raiz exata no extremo x0."

    if function(x1) == 0:
        return x1, 0, "Raiz exata no extremo x1."

    if not verificacao(x0, x1):
        return None, 0, "Intervalo inicial inválido."

    xr = x0

    for inte in range(1, intmax + 1):
        f0 = function(x0)
        f1 = function(x1)

        denominador = f1 - f0

        if denominador == 0:
            return xr, inte - 1, "Denominador nulo."

        # Fórmula da Falsa Posição
        xr = (x0 * f1 - x1 * f0) / denominador

        fxr = function(xr)

        if abs(fxr) < tol:
            return xr, inte, "Convergência atingida."

        if f0 * fxr < 0:
            x1 = xr
        else:
            x0 = xr

    return xr, intmax, "Máximo de iterações atingido."

print("\nFALSA POSICAO:")
print(falsa_posicao(x0, x1))



#Salvando os Dados:
dados_Bisseccao = pd.DataFrame({
    "Iteração": iter_lista,
    "xr": xr_lista,
    "f(xr)": fxr_lista
})

print(dados_Bisseccao)

dados_Bisseccao.to_csv("resultados_bissecao.txt", index=False)
