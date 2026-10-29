# ------- Biblioteca --------
import math as math

#----- Definições de Var -------
R1= 6197.01
R2= 29250.23
# ----- Inicio ----
Q10= 0.03 # m³/s 
Q20= 0.02 #m³/s
#correções


# definindo Funções
def f1(Q10, Q20):
    return (Q10 + Q20 - 0.05)

def f2(Q10, Q20):
    return (R1*Q10**2 - R2*Q20**2)

# ---- Jacobiana ----
J11=1
J12=1
J21=2*R1*Q10
J22=-2*R2*Q20

D = J11*J22 - J12*J21

r1= f1(Q10, Q20)
r2= f2(Q10,Q20)

deltaQ1= (r1*J22 + J12*r2)/D
deltaQ2= (J11*r2 - r1*J21)/D