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

#Números de pontos (linha y, em função de j), condição para a existência de m
if C == L:
    m = n
else:
    m = (n * L) / C

#deltax = deltay = delta
deltax = C / n
deltay = L / m
delta = deltax

#Cálculo dado no documento
omega = (alpha * deltat) / (2 * delta ** 2)

#Cria uma malha com os parâmetros (inicio, fim, n/m de "divisões" na malha)
#Obs.: Os valores de n e m precisam ser inteiros, de acordo com a função linspace()
x = np.linspace(a, C, n + 1)
y = np.linspace(a, L, m + 1)
t = np.linspace(a, tfinal, nt)

#Cria as matrizes A, B e C, já preenchidas
A = np.kron(np.eye(n - 1), 1.0 + 4.0 * omega) + np.kron(np.eye(n - 1, k = 1), - omega) + np.kron(np.eye(n - 1, k = -1), - omega)
B = np.kron(np.eye(n - 1), omega)
C = np.kron(np.eye(n - 1), 1.0 - 4.0 * omega) + np.kron(np.eye(n - 1, k = 1), omega) + np.kron(np.eye(n - 1, k = -1), omega)

#Cria as matrizes Gammas 1 (lado esquerdo) e 2 (lado direito) já preenchidas
G1 = np.kron(np.eye(m - 1), A) + np.kron(np.eye(m - 1, k = 1), -B) + np.kron(np.eye(m - 1, k = -1), -B)
G2 = np.kron(np.eye(m - 1), C) + np.kron(np.eye(m - 1, k = 1), B) + np.kron(np.eye(m - 1, k = -1), B)

#Cria e preenche a matriz Temp (Lambda) para condição inicial, k = t = 0
###T(x,y,0) = q(y,t)###
Temp = np.zeros(((n - 1) * (m - 1), 1))

for j in np.arange(0, m - 1):
    for i in np.arange(0, n - 1):
        k = i + j * (n - 1)
        Temp[k] = np.exp( - (x[i + 1] ** 2 + y[j + 1] ** 2) / L)
        
###Condições de contorno###
#Contorno direito T(C,y,t) = g(y,t)
#Cria a matriz CD (Contorno Direito)
CD = np.zeros((m + 1,nt)) 

#Contorno esquerdo T(0,y,t) = f(y,t)
#Cria a matriz CE (Contorno Esquerdo)
CE = np.zeros((m + 1,nt)) 

#Preenche as matrizes CD e CE
c0 = L ** 2

for k in range(nt):
    c1 = (c0 + 4.0 * alpha * t[k])
    c2 = c0 / c1
    for j in range(m + 1):
        CD[j,k] = c2 * np.exp(-(L ** 2 + y[j] ** 2) / c1)
        CE[j,k] = c2 * np.exp(-(y[j] ** 2) / c1)

#Contorno inferior T(x,0,t) = h(x,t)
#Cria a matriz CI (Contorno Inferior)
CI = np.zeros((n + 1, nt)) 

#Contorno superior T(x,L,t) = p(x,t)
#Cria a matriz CS (Contorno Superior)
CS = np.zeros((n + 1, nt)) 

#Preenche as matrizes CS e CI
for k in range(nt):
    c1 = (c0 + 4.0 * alpha * t[k])
    c2 = c0 / c1
    for j in range(n + 1):
        CS[j,k] = c2 * np.exp(-(L ** 2 + x[j] ** 2) / c1)
        CI[j,k] = c2 * np.exp(-(x[j] ** 2) / c1)

#Cria as matrizes dentro de Theta para condição inicial, k = t = 0
T1 = np.zeros(((m - 1) * (n - 1), 1))
T2 = np.zeros(((m - 1) * (n - 1), 1))

#Preenche as matrizes T1 e T2
def matrizesTheta(h, g):
    #Preenche a primeira matriz (CE e CD) T1
    for j in np.arange(0, m - 1):
        i, k = j * (n - 1), (j + 1) * (n - 1) - 1
        T1[i] = CE[j + 1,h] + CE[j + 1,h + 1]
        T1[k] = CD[j + 1,h] + CD[j + 1,h + 1]

    #Preenche a segunda matriz (CI e CS) T2
    for j in np.arange(0, n - 1):
        k = (m - 2) * (n - 1) + j
        T2[j] = CI[j + 1,g] + CI[j + 1,g + 1]
        T2[k] = CS[j + 1,g] + CS[j + 1,g + 1]

matrizesTheta(0, 0)

def resolveSis(Temp):
    #Cria a matriz Theta 
    T = omega * T1 + omega * T2

    #Resolve o sistema para o primeiro passo de tempo (k = t = 1)
    D = G2.dot(Temp) + T
    R = np.linalg.solve(G1, D)
    return R

resultado = resolveSis(Temp)

#Gera matrizes para os próximos passos de tempo k = 2 até k = nt
h = 0
g = 0
for p in np.arange(2, nt):
    Temp = resultado
    h += 1
    matrizesTheta(h, g)
    
    g += 1
    matrizesTheta(h, g)

    resultado = resolveSis(Temp) 

###Gráfico de Isolinhas###
#Gera matrizes para montar o gráfico de isolinhas
# X0, Y0 = np.meshgrid(x, y) 
# X = X0.transpose()
# Y = Y0.transpose()
X, Y = np.meshgrid(x, y) 

#Gera a matriz S (resultados das coordenadas de X e Y)
S = np.zeros((n + 1, m + 1))

#Preenche a matriz S
h = 1
v = 0
g = 1
for j in range(m + 1):
    for i in range(n + 1):
        if j == 0:
            S[i, j] = CI[i, nt - 1]   
        elif j > 0 and j < m:
            if i == 0:
                S[i, j] = CE[h, nt - 1]   
                h += 1
            elif i == n:     
                S[n, j] = CD[g, nt - 1]
                g += 1
            else:
                S[i, j] = resultado[v]
                v += 1      
        elif j == m:
            S[i, m] = CS[i, nt - 1] 
    v = (j * (n - 1))
        
print(resultado)
        
#Gera o gráfico de isolinhas    
fig, ax = plt.subplots()
Z = np.flip(S.transpose(), 0)
map = ax.imshow(Z, interpolation='bilinear', cmap='jet',  extent = [x[0], x[n], y[0], y[m]])
cb=plt.colorbar(map)
cb.set_label("Temperatura (R)")
plt.title('t = %.2f' %(tfinal))
plt.show()