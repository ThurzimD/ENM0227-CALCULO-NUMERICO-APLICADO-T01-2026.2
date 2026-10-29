#-------------------Eliminação Gaussiana--------------
import math as math 
# Estabelencdo os valores 

Q01=5
c01=10
Q03=8
c03=20
Q12=3
Q15=3
Q23=1
Q24=1
Q25=1
Q31=1
Q34=8
Q44=11
Q54=2
Q55=2

MatrizA=[
    [(Q12 + Q15), 0, -Q31, 0, 0],
    [-Q12, (Q23 + Q24 + Q25), 0, 0, 0],
    [0, -Q23, (Q31 + Q34), 0, 0],
    [0, -Q24, -Q34, Q44, -Q54],
    [-Q15, -Q25, 0, 0, (Q54 + Q55)]
]

MatrizB=[
    (Q01*c01),
    0,
    (Q03*c03),
    0,
    0
]

MatrizC=[
    0,
    0,
    0,
    0,
    0
]

n=5
it=0
for k in range(n-1):
    it+=1
    for i in range(k+1, n):
        fator=(MatrizA[i][k])/(MatrizA[k][k])
        for (j) in range(k+1,n):
            MatrizA[i][j]=MatrizA[i][j] - fator*(MatrizA[k][j])
        MatrizB[i]=(MatrizB[i])- fator*MatrizB[k]
        MatrizA[i][k]=0.0
    #print("Matriz a cada Pivo:","Pivo",it, MatrizA)


MatrizC[-1]=MatrizB[-1]/MatrizA[-1][-1]


for i in range(n-1, -1, -1):
    soma=MatrizB[i]
    for j in range(i+1, n):
        soma=soma-(MatrizA[i][j]*MatrizC[j])
    MatrizC[i]=soma/MatrizA[i][i]
print(MatrizC)