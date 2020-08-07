<<<<<<< HEAD
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

#Tamanho do passo de tempo
deltat = 0.5

#Número total de passos no tempo
nt = 5 

#Tempo final obtido ao fim de todas as iterações
tfinal = (nt - 1) * deltat

#Números de pontos (linha x, em função de i)
n = 5

if C == L:
    m = n
else:
    m = (n * L) / C

#deltax = deltay = delta
deltax = C / n
deltay = L / m
delta = deltax

omega = (alpha * deltat) / (2 * delta ** 2)

#Cria uma malha com os parâmetros (inicio, fim, n/m de "divisões" na malha)
x = np.linspace(a, C, n + 1)
y = np.linspace(a, L, m + 1)
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

for j in np.arange(0, m - 1):
    for i in np.arange(0, n - 1):
        k = i + j * (n - 1)
        Temp[k] = np.exp( - (x[i + 1] ** 2 + y[j + 1] ** 2) / L)
        
#Condições de contorno
#Contorno direito T(C,y,t) = g(y,t)
CD = np.zeros((m + 1,nt)) 

#Contorno esquerdo T(0,y,t) = f(y,t)
CE = np.zeros((m + 1,nt)) 

c0 = L ** 2

for k in range(nt):
    c1 = (c0 + 4.0 * alpha * t[k])
    c2 = c0 / c1
    for j in range(m + 1):
        CD[j,k] = c2 * np.exp(-(L ** 2 + y[j] ** 2) / c1)
        CE[j,k] = c2 * np.exp(-(y[j] ** 2) / c1)

#Contorno inferior T(x,0,t) = h(x,t)
CI = np.zeros((n + 1, nt)) 

#Contorno superior T(x,L,t) = p(x,t)
CS = np.zeros((n + 1, nt)) 

for k in range(nt):
    c1 = (c0 + 4.0 * alpha * t[k])
    c2 = c0 / c1
    for j in range(n + 1):
        CS[j,k] = c2 * np.exp(-(L ** 2 + x[j] ** 2) / c1)
        CI[j,k] = c2 * np.exp(-(x[j] ** 2) / c1)

#Preenche as matrizes dentro de Theta para condição inicial, k = t = 0
T1 = np.zeros(((m - 1) * (n - 1), 1))
T2 = np.zeros(((m - 1) * (n - 1), 1))

#Gera a primeira matriz (CE e CD)
for j in np.arange(0, m - 1):
    i, k = j * (n - 1), (j + 1) * (n - 1) - 1
    T1[i] = CE[j + 1,0] + CE[j + 1,1]
    T1[k] = CD[j + 1,0] + CD[j + 1,1]

#Gera a segunda matriz (CI e CS)
for j in np.arange(0, n - 1):
    k = (m - 2) * (n - 1) + j
    T2[j] = CI[j + 1,0] + CI[j + 1,1]
    T2[k] = CS[j + 1,0] + CS[j + 1,1]

#Cria a matriz Theta 
T = omega * T1 + omega * T2

#Resolve o sistema para o primeiro passo de tempo (k = t = 1)
D = G2.dot(Temp) + T
R = np.linalg.solve(G1, D)

#Gera matrizes para os próximos passos de tempo k = 2 até k = nt
h = 0
g = 0
for p in np.arange(2, nt):
    Temp = R
    h += 1
    for j in np.arange(0, m - 1):
        i, k = j * (n - 1), (j + 1) * (n - 1) - 1
        T1[i] = CE[j + 1,h] + CE[j + 1,h + 1]
        T1[k] = CD[j + 1,h] + CD[j + 1,h + 1]
    
    g += 1
    for j in np.arange(0, n - 1):
        k = (m - 2) * (n - 1) + j
        T2[j] = CI[j + 1,g] + CI[j + 1,g + 1]
        T2[k] = CS[j + 1,g] + CS[j + 1,g + 1]

    T = omega * T1 + omega * T2
    D = G2.dot(Temp) + T
    R = np.linalg.solve(G1, D)    

#Gera matrizes para montar o gráfico de isolinhas
X0, Y0 = np.meshgrid(x, y) 
X = X0.transpose()
Y = Y0.transpose()

#Gera a matriz S (resultados das coordenadas de X e Y)
S = np.zeros((n + 1, m + 1))

#Preenche a matriz S
h = 1
v = 0
g = 1
for j in range(m + 1):
    for i in range(n + 1):
        if j == 0:
            S[i, j] = CI[i, 1]   
        elif j > 0 and j < m:
            if i == 0:
                S[i, j] = CE[h, 1]   
                h += 1
            elif i == n:     
                S[n, j] = CD[g, 1]
                g += 1
            else:
                S[i, j] = R[v]
                v += 1      
        elif j == m:
            S[i, m] = CS[i, 1] 
    v = (j * (n - 1))
        
#Faz o gráfico de isolinhas    
map=plt.contourf(X,Y,S,cmap='jet')
cb=plt.colorbar(map)
cb.set_label("Temperatura (R)")
plt.title('Para o último passo de tempo (t = \u0394t = %.2f)' %(deltat))
plt.show()
=======
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
# m = 5

#deltax = deltay = delta
deltax = C / n
deltay = L / m
delta = deltax

omega = (alpha * deltat) / (2 * delta ** 2)

#Cria uma malha com os parâmetros (inicio, fim, n/m de "divisões" na malha)
x = np.linspace(a, C, n + 1)
y = np.linspace(a, L, m + 1)
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

for j in np.arange(0, m - 1):
    for i in np.arange(0, n - 1):
        k = i + j * (n - 1)
        Temp[k] = np.exp( - (x[i + 1] ** 2 + y[j + 1] ** 2) / L)
        
#Condições de contorno
#Contorno direito T(C,y,t) = g(y,t)
CD = np.zeros((m + 1,nt)) 

#Contorno esquerdo T(0,y,t) = f(y,t)
CE = np.zeros((m + 1,nt)) 

c0 = L ** 2

for k in range(nt):
    c1 = (c0 + 4.0 * alpha * t[k])
    c2 = c0 / c1
    for j in range(m + 1):
        CD[j,k] = c2 * np.exp(-(L ** 2 + y[j] ** 2) / c1)
        CE[j,k] = c2 * np.exp(-(y[j] ** 2) / c1)

#Contorno inferior T(x,0,t) = h(x,t)
CI = np.zeros((n + 1, nt)) 

#Contorno superior T(x,L,t) = p(x,t)
CS = np.zeros((n + 1, nt)) 

for k in range(nt):
    c1 = (c0 + 4.0 * alpha * t[k])
    c2 = c0 / c1
    for j in range(n + 1):
        CS[j,k] = c2 * np.exp(-(L ** 2 + x[j] ** 2) / c1)
        CI[j,k] = c2 * np.exp(-(x[j] ** 2) / c1)

#Preenche as matrizes dentro de Theta para condição inicial, k = t = 0
T1 = np.zeros(((m - 1) * (n - 1), 1))
T2 = np.zeros(((m - 1) * (n - 1), 1))

#Gera a primeira matriz (CE e CD)
for j in np.arange(0, m - 1):
    i, k = j * (n - 1), (j + 1) * (n - 1) - 1
    T1[i] = CE[j + 1,0] + CE[j + 1,1]
    T1[k] = CD[j + 1,0] + CD[j + 1,1]

#Gera a segunda matriz (CI e CS)
for j in np.arange(0, n - 1):
    k = (m - 2) * (n - 1) + j
    T2[j] = CI[j + 1,0] + CI[j + 1,1]
    T2[k] = CS[j + 1,0] + CS[j + 1,1]

#Cria a matriz Theta 
T = omega * T1 + omega * T2

#Resolve o sistema para o primeiro passo de tempo (k = t = 0)
D = G2 * Temp + T
R = np.linalg.solve(G1,D)

print(R)




















# #Gera o gráfico de linhas
# fig = plt.figure(1)
# ax = fig.add_subplot(111)
# for i in range(m):
# 	ax.plot(x[0:n], F[i])
# ax.set_xlabel('Passos de tempo')   
# ax.set_ylabel('Temperatura') 
# plt.show()
>>>>>>> 0cb098fdd50c95e76a06999fc242b7b1f9e1a9cb
