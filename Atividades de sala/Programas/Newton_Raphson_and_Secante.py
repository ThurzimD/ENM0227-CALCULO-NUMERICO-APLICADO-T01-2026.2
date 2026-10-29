#----Bibliotecas----
import math as math

#----Deifinição das Variaveis----
A = 0.01 #Area (m²)
h = 10 #Coeficiente de transferencia de calor (W/(m²K))
emi = 0.8 #emissividade da superficie (admensional)
sigma = 5.67*(10**-8) #Constante de Stefan Boltzmann (W/(m²K))
Tarr = 300 # Temperatura ambiente (K)
P = 20 #potencia termica (W)


#-----Chute inicial------
Ti = 350 # (K)

#---Lista para armazenar dados----
lista_it=[]
lista_Ts=[]
lista_f=[]
lista_erro=[]

#----Definição de Funções----
def f(Ti):
    return (h*A*(Ti-Tarr))+(emi*sigma*A*((Ti**4)-(Tarr**4)))-P

def F(Ti):
    return (4*emi*sigma*A*(Ti**3))+(h*A)

#------ Metodo de Newton Raphson ------

Erro = 1
it= 0
maxIt = 100

while (Erro > (1e-6) and it<=maxIt):
    fun=f(Ti)
    fun_d = F(Ti)
    Ts = Ti - (fun/fun_d)
    Erro = abs(Ts-Ti)
    it +=1
    Ti=Ts
    lista_it.append(it)
    lista_Ts.append(Ts)
    lista_f.append(fun)
    lista_erro.append(Erro)
else:
    print(lista_it)
    print(lista_Ts)
    print(lista_f)
    print(lista_erro)

    print("Temperatura de Equilibrio encontrada:", Ts, "K")
    print("Numero de interações necessarias:", it)
    print("Finalizado")


#----- Armazenamento de dados------



#----- Metodo da Secante-----

Erro = 1
it= 0
maxIt= 100
Ts=0
T0 = 400 #K
T1 = 425 #K

while (Erro > (1e-6) and it<=maxIt):
    Ts = T1 - (f(T1)*((T1 - T0)/(f(T1)-f(T0))))
    Erro = abs(Ts-T1)
    it +=1
    T0=T1
    T1=Ts
    print(Ts)
else:
    print("Finalizado Secante")