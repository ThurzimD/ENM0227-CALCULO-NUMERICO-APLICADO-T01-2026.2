# Codigo para teste do metodo numerico de Runger Kutter 4 (RK4)

#pip install "biblioteca" (Lembrete)

#----------Bibliotecas-----
import math 
import matplotlib.pyplot as plt

# ---------- Definição da Função da EDO ----------
def f(x, yi, Res, St):
    return (1 - yi - (3/8)*Res*(yi**2)) / St # definição da EDO dy/dx = f(x, y), equação adimensionalizada.

#x   -> tempo adimensional t* (não é usado aqui pq a equação é autônoma,mas mantemos no argumento pra seguir o padrão f(x,y) do RK4);
# yi  -> velocidade adimensional v_z* no instante atual;
# Res -> número de Reynolds de partícula; 
# Res=0 -> regime de Stokes puro;
# St  -> número de Stokes;

# ---------- Método de RK4 ----------
def frk4(x, yi, h, Res, St):
    k1 = f(x, yi, Res, St) #Primeiro termo (ínico do passo);
    k2 = f(x + h/2, yi + (h/2)*k1, Res, St) #Segundo termo;
    k3 = f(x + h/2, yi + (h/2)*k2, Res, St) #Terceiro termo;
    k4 = f(x + h, yi + h*k3, Res, St) #Quarto termo (fim do passo);
    return yi + (h/6)*(k1 + 2*k2 + 2*k3 + k4)

# ---------- Roda a simulação completa ---------- LOOP
def simula(St, Res, h, t_final):
    x, yi = 0.0, 0.0
    resultados = [(x, yi)]
    while x < t_final:
        yi = frk4(x, yi, h, Res, St)
        x += h
        resultados.append((x, yi))
    return resultados

# ---------- Solução analítica (Re -> 0, Stokes puro) ----------
def analitica_stokes(t, St):
    return 1 - math.exp(-t/St)


# =====================  1: Re -> 0, comparar vários St ================
for St in [0.5, 1.0, 2.0]:
    num = simula(St=St, Res=0.0, h=0.01, t_final=10.0)
    x_last, y_last = num[-1]
    y_exato = analitica_stokes(x_last, St)
    print(f"St={St}: y_num={y_last:.5f}, y_exato={y_exato:.5f}, erro={abs(y_last-y_exato):.2e}")


# ===================== 2: variar o passo h =====================
St_fixo = 1.0
for h in [0.1, 0.05, 0.01, 0.001]:
    num = simula(St=St_fixo, Res=0.0, h=h, t_final=5.0)
    x_last, y_last = num[-1]
    y_exato = analitica_stokes(x_last, St_fixo)
    print(f"h={h}: erro={abs(y_last-y_exato):.2e}")


# ===================== 3: caso com Re != 0 =====================
resultado_oseen = simula(St=1.0, Res=1.0, h=0.01, t_final=10.0)


# ==================== 4: validar com solução exata de Riccati  =====================
def analitica_riccati(t, St, Res):
    v_star = (-1 + math.sqrt(1 + 1.5*Res)) / (0.75*Res)
    Q = 3/(8*St*Res)
    P = -(3/(4*St*Res))*v_star - 1/St
    return v_star + Q/P + (-1/v_star - Q/P) * math.exp(-P*t) ** -1

for St in [1.0]:
    for Res in [0.5, 1.0, 2.0]:
        num = simula(St=St, Res=Res, h=0.01, t_final=10.0)
        x_last, y_last = num[-1]
        print(f"St={St}, Res={Res}: y_num={y_last:.5f}")


# ===================== ITEM 5: plotar para vários Res =====================


plt.figure()
for Res in [0.0, 0.5, 1.0, 2.0, 5.0]:
    dados = simula(St=1.0, Res=Res, h=0.01, t_final=10.0)
    xs = [p[0] for p in dados]
    ys = [p[1] for p in dados]
    plt.plot(xs, ys, label=f"Res={Res}")

plt.xlabel("t*")
plt.ylabel("v_z*")
plt.legend()
plt.title("Efeito do número de Reynolds na sedimentação")
plt.show()
