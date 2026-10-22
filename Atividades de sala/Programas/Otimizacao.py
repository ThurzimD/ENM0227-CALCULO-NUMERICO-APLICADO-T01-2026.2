#Atividade Prática de Programação Otimização Unidimensional da Espessura de Isolamento Térmico

#Bibliotecas
import math as math


#Parâmetros do caso

r1 = 0.04 #raio externo da tubulação sem isolamento (m)
L = 20 # Comprimento da tubulação (m)
Ti = 180 # Temperatura na interface tubo–isolante (°C)
T_inf = 40 # Temperatura ambiente (°C)
k = 0.045 # Condutividade térmica do isolante ( W/(m*K))
h = 12 # Coeficiente de convecção externa (W/(m^2 * K))
t = 30000 # Tempo acumulado de operação (h)
Pe = 0.70 # Preço da energia térmica reposta (R$/kWh)
ci = 9000 # Custo instalado do isolante por volume (R$/m^3)
#e = (0 , 0.12) # Intervalo admissível (m)

def Tx_de_perda(e):
    Rcond = (math.log((r1 + e) / r1)) / (2*k*(math.pi))
    Rconv = 1 / (2*math.pi*h*(r1 + e))
    q = (Ti - T_inf) / (Rcond + Rconv) # (W/m)
    return q

#fluxo_de_calor = Tx_de_perda(e[1])

#Q_potencia = L*(Tx_de_perda(e[1]))

def custo(e):
    C_e = Pe * ((L * Q_potencia * t)/1000 ) #Custo de energia acumulado no período

    C_i = ci * math.pi * L * ((r1 + e)**2 - r1**2) # Custo Instalado

    C = C_e + C_i
    return C, C_i, C_e

e = 0.0
while e <= 0.12:
    fluxo_de_calor = Tx_de_perda(e)
    Q_potencia = L*(Tx_de_perda(e))
    Custo_Total, Custo_i, Custo_e=custo(e)
    print(f"Custo total para espessura {e} é: {Custo_Total:.2f}!")
    print(f"Custo de Energia: {Custo_e}!")
    print(f"Custo Instalado: {Custo_i}!")
    e = round(e + 0.02, 2)

#Implementação do método da razão áurea 

#   Preparar o intervalo e a tolerância

tol = 10e-6 
int_max = 200
a = 0.0
b= 0.12
tensao = (math.sqrt(5) - 1) / 2 

# intervalos
x1 = b - (tensao * (b - a))
x2 = a + (tensao * (b - a))

# Programar a redução sucessiva do intervalo

inte = 0
#while abs(b-a) > tol and inte < int_max:

