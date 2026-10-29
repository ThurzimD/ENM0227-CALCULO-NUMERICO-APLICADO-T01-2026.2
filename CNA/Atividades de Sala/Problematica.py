import math
import csv
import pandas as pd

# =============================================================================
# TOOLKIT MESTRE DE CÁLCULO NUMÉRICO APLICADO (CNA) - PYTHON PURO
# Feito para a Prova: Altere apenas as funções no final do arquivo!
# =============================================================================

# -----------------------------------------------------------------------------
# MÓDULO 1: ZEROS DE FUNÇÕES
# -----------------------------------------------------------------------------
def bissecao(f, a, b, tol=1e-6, max_iter=100):
    historico = []
    a_curr, b_curr = float(a), float(b)
    
    for i in range(1, max_iter + 1):
        c = (a_curr + b_curr) / 2.0
        fc = f(c)
        
        # Salvando dados para o Pandas/CSV
        historico.append([i, a_curr, b_curr, c, fc, (b_curr - a_curr) / 2.0])
        
        if abs(fc) < tol or (b_curr - a_curr) / 2.0 < tol:
            break
            
        if f(a_curr) * fc < 0:
            b_curr = c
        else:
            a_curr = c
            
    df = pd.DataFrame(historico, columns=['Iter', 'a', 'b', 'Raiz(c)', 'f(c)', 'Erro'])
    return c, df

def newton_raphson(f, df_dx, x0, tol=1e-6, max_iter=100):
    historico = []
    x = float(x0)
    
    for i in range(1, max_iter + 1):
        fx = f(x)
        derivada = df_dx(x)
        
        if abs(derivada) < 1e-12:
            print("Erro: Derivada nula.")
            break
            
        dx = fx / derivada
        x_novo = x - dx
        
        historico.append([i, x, fx, derivada, x_novo, abs(dx)])
        x = x_novo
        
        if abs(dx) < tol:
            break
            
    df = pd.DataFrame(historico, columns=['Iter', 'x_antigo', 'f(x)', 'f\'(x)', 'x_novo', 'Erro'])
    return x, df

# -----------------------------------------------------------------------------
# MÓDULO 2: SISTEMAS LINEARES E MATRIZES (USANDO LISTAS NATIVAS)
# -----------------------------------------------------------------------------
def decomposicao_lu(A_orig, b_orig):
    """ Fatoração LU de Doolittle sem numpy """
    n = len(b_orig)
    # Criando matrizes vazias de zeros
    L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    U = [[0.0 for j in range(n)] for i in range(n)]
    
    # 1. Fatoração
    for k in range(n):
        for j in range(k, n):
            soma = sum(L[k][i] * U[i][j] for i in range(k))
            U[k][j] = A_orig[k][j] - soma
            
        for i in range(k + 1, n):
            soma = sum(L[i][j] * U[j][k] for j in range(k))
            L[i][k] = (A_orig[i][k] - soma) / U[k][k]
            
    # 2. Substituição Progressiva (L * y = b)
    y = [0.0] * n
    for i in range(n):
        soma = sum(L[i][j] * y[j] for j in range(i))
        y[i] = b_orig[i] - soma
        
    # 3. Substituição Regressiva (U * x = y)
    x = [0.0] * n
    x[n - 1] = y[n - 1] / U[n - 1][n - 1]
    for i in range(n - 2, -1, -1):
        soma = sum(U[i][j] * x[j] for j in range(i + 1, n))
        x[i] = (y[i] - soma) / U[i][i]
        
    return x, L, U

def thomas_tdma(e_in, f_in, g_in, r_in):
    """ Algoritmo de Thomas para Matrizes Tridiagonais """
    n = len(r_in)
    e, f, g, r = list(e_in), list(f_in), list(g_in), list(r_in)
    
    for k in range(1, n):
        e[k] = e[k] / f[k-1]
        f[k] = f[k] - e[k] * g[k-1]
        r[k] = r[k] - e[k] * r[k-1]
        
    x = [0.0] * n
    x[n-1] = r[n-1] / f[n-1]
    for k in range(n-2, -1, -1):
        x[k] = (r[k] - g[k] * x[k+1]) / f[k]
        
    return x

# -----------------------------------------------------------------------------
# MÓDULO 3: EQUAÇÕES DIFERENCIAIS ORDINÁRIAS (EDO)
# -----------------------------------------------------------------------------
def rk4(f_edo, t0, tf, y0, h):
    """ Método de Runge-Kutta 4ª Ordem """
    historico = []
    t = float(t0)
    y = float(y0)
    passos = int((tf - t0) / h)
    
    historico.append([t, y])
    
    for _ in range(passos):
        k1 = f_edo(t, y)
        k2 = f_edo(t + 0.5*h, y + 0.5*k1*h)
        k3 = f_edo(t + 0.5*h, y + 0.5*k2*h)
        k4 = f_edo(t + h, y + k3*h)
        
        y = y + (h / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4)
        t = t + h
        historico.append([round(t, 6), y])
        
    df = pd.DataFrame(historico, columns=['Tempo(t)', 'y(t)'])
    return y, df

# -----------------------------------------------------------------------------
# MÓDULO 4: OTIMIZAÇÃO UNIDIMENSIONAL
# -----------------------------------------------------------------------------
def razao_aurea(f_obj, a, b, tol=1e-6, max_iter=100):
    historico = []
    tau = (math.sqrt(5.0) - 1.0) / 2.0
    
    x1 = b - tau * (b - a)
    x2 = a + tau * (b - a)
    f1, f2 = f_obj(x1), f_obj(x2)
    
    for i in range(1, max_iter + 1):
        historico.append([i, a, b, (a+b)/2.0, abs(b-a)])
        if abs(b - a) < tol:
            break
            
        if f1 > f2:
            a = x1; x1 = x2; f1 = f2
            x2 = a + tau * (b - a)
            f2 = f_obj(x2)
        else:
            b = x2; x2 = x1; f2 = f1
            x1 = b - tau * (b - a)
            f1 = f_obj(x1)
            
    x_opt = (a + b) / 2.0
    df = pd.DataFrame(historico, columns=['Iter', 'a', 'b', 'x_otimo', 'Erro'])
    return x_opt, f_obj(x_opt), df

# =============================================================================
# ÁREA DE EDIÇÃO PARA A PROVA: ALTERE AS FUNÇÕES AQUI
# =============================================================================
def f_raiz(x):
    return x**2 - 4.0        # Ex: Achar raiz de x^2 - 4 = 0

def df_raiz(x):
    return 2.0 * x           # Derivada analítica para Newton

def f_edo(t, y):
    return -2.0 * y + t      # Ex: dy/dt = -2y + t

def f_custo(x):
    return (x - 2.5)**2 + 3.0 # Ex: Minimizar parabóla

# =============================================================================
# TESTES E SALVAMENTO DE DADOS (USANDO PANDAS E CSV)
# =============================================================================
if __name__ == "__main__":
    print("--- INICIANDO SOLUÇÕES NUMÉRICAS ---")
    
    # 1. Teste de Raiz e exportação para CSV
    raiz, df_bissecao = bissecao(f_raiz, 1.0, 3.0)
    print(f"\n[RAIZ] Bisseção: x = {raiz:.6f}")
    df_bissecao.to_csv("resultado_bissecao.csv", index=False)
    print("-> Tabela 'resultado_bissecao.csv' gerada com sucesso via Pandas.")
    
    # 2. Teste de Sistema Linear (Exportação usando módulo 'csv' nativo)
    A = [[3.0, -0.1, -0.2], [0.1, 7.0, -0.3], [0.3, -0.2, 10.0]]
    b = [7.85, -19.30, 71.40]
    
    x_lu, _, _ = decomposicao_lu(A, b)
    print(f"\n[SISTEMA LINEAR] Solução LU: {x_lu}")
    
    with open("resultado_sistema_lu.csv", mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Variavel", "Valor"])
        for i, val in enumerate(x_lu):
            writer.writerow([f"x{i+1}", val])
    print("-> Arquivo 'resultado_sistema_lu.csv' gerado via módulo CSV nativo.")
    
    # 3. Teste de EDO
    y_final, df_rk4 = rk4(f_edo, 0.0, 2.0, 1.0, 0.1)
    print(f"\n[EDO - RK4] Solução final em t=2.0: y = {y_final:.6f}")
    df_rk4.to_csv("resultado_rk4.csv", index=False)
    
    # 4. Teste de Otimização
    x_otimo, custo_minimo, df_otim = razao_aurea(f_custo, 0.0, 4.0)
    print(f"\n[OTIMIZAÇÃO] Razão Áurea: x_otimo = {x_otimo:.6f}, Custo = {custo_minimo:.6f}")
    df_otim.to_csv("resultado_otimizacao.csv", index=False)
