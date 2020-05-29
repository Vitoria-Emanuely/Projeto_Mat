import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#Início
a = 0.0

#Fim
b = 10.0

#Números de pontos (linha x, em função de i)
n = 50

#Cria uma malha com os parâmetros (inicio, fim, n de "divisões" na malha)
x = np.linspace(a, b, n)

deltax = (b - a) / (n - 1)

#Parâmetros da equação lambda
alfa = 1.0
deltat = 1.5

lambdaa = (alfa * deltat) / deltax ** 2

#Número de pontos (linha y, em função de j)
m = 5

#Cria uma matriz de zeros NXN
A = np.zeros((n - 2, n - 2))
B = np.zeros(n - 2)

#Preenche a matriz fixa A
A[0, 0] = -(2.0 * lambdaa + 1.0)
A[0, 1] = lambdaa
A[1, 0] = lambdaa
A[n - 3, n - 3] = -(2.0 * lambdaa + 1.0)
A[n - 3, n - 4] = lambdaa
A[n - 4, n - 3] = lambdaa

for i in np.arange(1, n - 3):
    A[i, i] = -(2.0 * lambdaa + 1.0)
    A[i, i + 1] = lambdaa
    A[i + 1, i] = lambdaa

#Matriz B: Primeiro passo de tempo j=0
B[0] = - 0.0 - lambdaa * 100.0
B[n - 3] = - 0.0 - lambdaa * 50.0

for i in np.arange(1, n - 3):
    B[i] = 0.0

#Resolve a matriz linear
C = np.linalg.solve(A,B)
F = C

#Matriz B: Próximos passos de tempo j=1 até j=m
for j in np.arange(1, m):
    B[0] = - C[0] - lambdaa * 100.0
    B[n - 3] = - C[n - 3] - lambdaa * 50.0

    for i in np.arange(1, n - 3):
        B[i] = - C[i]
    C = np.linalg.solve(A,B)   
    #Insere na matriz F as matrizes salvas temporariamente em C
    F = np.vstack([F, C])

print(F)

#Gera o DataFrame 
z = []
for i in np.arange(1, n - 1):
    z.append(i)

y = []
for i in np.arange(0, m):
    y.append(i)

df = pd.DataFrame(F, index=y, columns=z)
print(df)

#Preenche a matriz F com 100°C e 50°C
K = np.zeros((m, 1))
for i in np.arange(0, m):
    K[i] = 50.0

F = np.hstack([F, K]) 

for i in np.arange(0, m):
    K[i] = 100.0

F = np.hstack([K, F])     

#Gera o gráfico de linhas
fig = plt.figure(1)
ax = fig.add_subplot(111)
for i in range(m):
	ax.plot(x[0:n], F[i])
ax.set_xlabel('Passos de tempo')   
ax.set_ylabel('Temperatura') 
plt.show()