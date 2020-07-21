import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp

#Início
a = 0.0

#Fim (Comprimento)
C = 10.0

#Variáveis
L = 10.0
alpha = 1.0 
deltat = 0.5
nt = 5

#Tempo final que será calculado
tfinal = (nt - 1) * deltat

#Números de pontos (linha x, em função de i)
n = 5

if C == L:
    m = n
else:
    m = (n * L) / C

#Número de pontos (linha y, em função de j) PRECISA SER INTEIRO POR CAUSA DA FUNÇÃO LINSPACE
m = 5

#deltax = deltay = delta
deltax = C / n
deltay = L / m
delta = deltax

omega = (alpha * deltat) / (2 ** delta ** 2)

#Cria uma malha com os parâmetros (inicio, fim, n/m de "divisões" na malha)
x = np.linspace(a, C, n)
y = np.linspace(a, L, m)
t = np.linspace(a, tfinal, nt)

#Cria uma matriz de zeros NXN
A = np.zeros((n - 1, n - 1))
B = np.zeros((n - 1, n - 1))
C = np.zeros((n - 1, n - 1))

#Preenche a matriz fixa A
A[n - 2, n - 2] = 1.0 + 4.0 * omega

for i in np.arange(0, n - 2):
    A[i, i] = 1.0 + 4.0 * omega
    A[i, i + 1] = - omega
    A[i + 1, i] = - omega

#Preenche a matriz fixa B
for i in np.arange(0, n - 1):
    B[i, i] = omega

#Preenche a matriz fixa C
C[n - 3, n - 2] = - omega
C[n - 2, n - 3] = - omega
C[n - 2, n - 2] = 1.0 - 4.0 * omega

for i in np.arange(0, n - 2):
    C[i, i] = 1.0 - 4.0 * omega
    C[i, i + 1] = omega
    C[i + 1, i] = omega

#Cria as matrizes Gammas 1 (lado esquerdo) e 2 (lado direito)
G1 = np.kron(np.eye(m - 1), A) + np.kron(np.eye(m - 1, k = 1), -B) + np.kron(np.eye(m - 1, k = -1), -B)

G2 = np.kron(np.eye(m - 1), C) + np.kron(np.eye(m - 1, k = 1), B) + np.kron(np.eye(m - 1, k = -1), B)

#Preenche a matriz Temp (Lambda) para condição inicial, k = t = 0
#T(x,y,0) = q(y,t)
Temp = np.zeros(((n - 1) * (m - 1), 1))

for j in np.arange(0, m - 2):
    for i in np.arange(0, n - 2):
        k = i + j * (n - 1)
        Temp[k] = np.exp( - (x[i + 1] ** 2 + y[j + 1] ** 2) / L)
        
#Condições de contorno
#Contorno direito T(C,y,t) = g(y,t)
CD = np.zeros((m,nt)) 

#Contorno esquerdo T(0,y,t) = f(y,t)
CE = np.zeros((m,nt)) 

c0 = L ** 2

for k in range(nt):
  c1 = (c0 + 4.0 * alpha * t[k])
  c2 = c0 / c1
  for j in range(m):
    CD[j,k] = np.exp(-(L * 2 + y[j] * 2) / c2)
    CE[j,k] = np.exp(-(y[j] ** 2) / c2)

#Contorno inferior T(x,0,t) = h(x,t)
CI = np.zeros((n, nt)) 

#Contorno superior T(x,L,t) = p(x,t)
CS = np.zeros((n, nt)) 

for k in range(nt):
  c1 = (c0 + 4.0 * alpha * t[k])
  c2 = c0 / c1
  for j in range(n):
    CS[j,k] = np.exp(-(L * 2 + x[j] * 2) / c2)
    CI[j,k] = np.exp(-(x[j] ** 2) / c2)

#Preenche a matriz Theta para condição inicial, k = t = 0
T = np.zeros(((m - 1) * (n - 1), 1))

# for i in np.arange(0, )


























# #Gera o gráfico de linhas
# fig = plt.figure(1)
# ax = fig.add_subplot(111)
# for i in range(m):
# 	ax.plot(x[0:n], F[i])
# ax.set_xlabel('Passos de tempo')   
# ax.set_ylabel('Temperatura') 
# plt.show()