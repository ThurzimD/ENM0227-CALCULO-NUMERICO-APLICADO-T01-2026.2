#Atividade 4 - metodod de Miller
import math as math
# dado a equação: 
#  (Q**2*T)/(g*A**3) = 1 
# 
#-Armazeno de dados-
lis_it=[]



#---------Variaveis---------
b=4 #m
z=1.5 
g=9.81 #m/s²

a0=-b
a1=-2*z
a3=g*(b**3)
a4=3*g*(b**2)*z
a5=3*g*b*(z**2)
a6=g*(z**3)

def f_polinomial(H, Q):
    return (a6*(H**6))+(a5*(H**5))+(a4*(H**4))+(a3*(H**3))+(a1*(Q**2)*(H))+a0*(Q**2)

################ MILLER #############
#--Codição inicial---
H0=0.5
H1=1
H2=1.5
Q=20
it= 0

h0=H1-H0
h1=H2-H1

f0=f_polinomial(H0, Q)
f1=f_polinomial(H1, Q)
f2=f_polinomial(H2, Q)

delta0= (f1-f0)/h0
delta1= (f2-f1)/h1

a=(delta1-delta0)/(h1+h0)
bm= a*h1 + delta1

c=f2

### Discriminante ###
D= math.sqrt((bm**2)-4*a*c)

if (abs(bm+D)> abs(bm-D)):
    E=bm+D
else:
    E=bm-D

delta_H= -(2*c/E)

H3=H2+delta_H

#### Interações ####
Er=1
while Er > 10**-8:
    H0 = H1
    H1 = H2
    H2 = H3
    h0=H1-H0
    h1=H2-H1
    f0=f_polinomial(H0, Q)
    f1=f_polinomial(H1, Q)
    f2=f_polinomial(H2, Q)
    f3=f_polinomial(H3, Q)

    delta0= (f1-f0)/h0
    delta1= (f2-f1)/h1

    a=(delta1-delta0)/(h1-h0)
    bm= a*h1 + delta1

    c=f2

    ### Discriminante ###
    D= math.sqrt((bm**2)-4*a*c)

    if (abs(bm+D)> abs(bm-D)):
        E=bm+D
    else:
        E=bm-D

    delta_H= -(2*c/E)

    H3=H2+delta_H
    Er=abs(H3-H2)
    it+=1

    print(it)
    print(H3)
    print(f3)
    print(Er)

