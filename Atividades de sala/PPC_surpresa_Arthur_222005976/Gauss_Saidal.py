#Bibliotecas
import math as math
import pandas as pd



#Gauss Seidel Problema da Placa
# Algoritimo

T=[]

# Definindo malhas e as condições de contorno

n=20
T= [[0 for _ in range(n+2)] for _ in range(n+2)]
for i in range(n+2):
    for k in range(n+2):
        T[i][k] = 20
        T[n+1][k] = 100


# Implementação da atualização interativa
inter = 0
max_inter=100000
tol=10e-6 #°C
erro = tol*2
interacoes=[]


while (erro > tol and inter < max_inter):
    erro = 0.0
    for j in range(1, n+1):
        for i in range(1,n+1):
            antigo = T[i][j]
            T[i][j] = ((T[i-1][j] + T[i+1][j] + T[i][j-1] + T[i][j+1])) / 4
            diferenca = abs(T[i][j]-antigo)
            erro = max(erro, diferenca)
    inter = inter + 1

r=[[0 for _ in range(n+2)] for _ in range(n+2)]
maxr = 0
for j in range(n+1):
    for i in range(n+1):
        r[i][j] = (4 * (T[i][j])) - (T[i-1][j] + T[i+1][j] + T[i][j-1] + T[i][j+1])
        maxr = max(abs(r[i][j]), abs(maxr))
print(maxr)



#Armazenamento de dados
dados = {
    'Nº de Interação total': [inter],
    'Nº do erro final': [erro],
    f'Temperatura do nó central para n={n}': [T[int((n+1)/2)][int((n+1)/2)]],
    'Valor de Rij': [maxr],
}

df1 = pd.DataFrame(dados)
df1.to_csv(f'Dados_{n+2}x{n+2}.txt', sep='\t', index=False, encoding = 'utf-8')

df2 = pd.DataFrame(T)
df2.to_csv(f'Matriz_{n+2}x{n+2}_Temperatura_.txt', sep='\t', index=False, encoding = 'utf-8')



