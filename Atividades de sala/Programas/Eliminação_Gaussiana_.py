#Eliminação Gaussiana por Pivoteamento 

# Ampliação para um sistema 5 × 5

#AF = b
import math as math


d1=15
d2=65
d3= 45


Matriz_A=[
    [math.sin(d1), math.sin(d2), 1, 0, 0],
    [ math.cos(d1), math.cos(d2), 0, 0, 0],
     [0, 0, -1, math.cos(d3), 0],
     [0, 0, 0, -(math.sin(d3)), -1],
     [0, 0, 0, 0, -1]
]

Matriz_F=[
        [100],
        [-80],
        [20],
        [35],
        [15]

]


Matriz_b=[
    [0.0],
    ]*5

for i in range(5):
    Matriz_b[i] = Matriz_A[i][i]*Matriz_F[i][0]

n=5
it=0

for k in range(n):
    for i in range(k+1, n):
        fator=(Matriz_A[i][k])/(Matriz_A[k][k])
        for (j) in range(k+1,n):
            Matriz_A[i][j]=Matriz_A[i][j] - fator*(Matriz_A[k][k])
        Matriz_b[i]=(Matriz_b[i])- fator*Matriz_b[k]
        Matriz_A[i][k]=0.0


