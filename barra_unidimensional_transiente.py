import numpy as np
import matplotlib.pyplot as plt

#Início
a = 0.0

#Fim
b = 10.0

#Números de pontos (linha x, em função de i)
n = 10

#Cria uma malha com os parâmetros (inicio, fim, n de "divisões" na malha)
x = np.linspace(a, b, n)

deltax = (b - a) / (n - 1)

#Parâmetros da equação lambda
alfa = 1.0
deltat = 0.5

lambdaa = (alfa * deltat) / deltax ** 2

#Número de pontos (linha y, em função de j)
m = 8

#Cria uma matriz de zeros NXN
A = np.zeros((n - 1, n - 1))
B = np.zeros(n - 1)

#Preenche a matriz fixa A
A[0, 0] = -(2.0 * lambdaa + 1.0)
A[0, 1] = lambdaa
A[1, 0] = lambdaa
A[n - 2, n - 2] = -(2.0 * lambdaa + 1.0)

for i in np.arange(1, n - 2):
    A[i, i] = -(2.0 * lambdaa + 1.0)
    A[i, i + 1] = lambdaa
    A[i + 1, i] = lambdaa

#Matriz B: Primeiro passo de tempo j=0
B[0] = - 0.0 - lambdaa * 100.0
B[n - 2] = - 0.0 - lambdaa * 50.0

for i in np.arange(1, n - 2):
    B[i] = 0.0

#Resolve a matriz linear
C = np.linalg.solve(A,B)
F = C

#Matriz B: Próximos passos de tempo j=1 até j=m
for j in np.arange(1, m):
    B[0] = - C[0] - lambdaa * 100.0
    B[n - 2] = - C[n - 2] - lambdaa * 50.0

    for i in np.arange(1, n - 2):
        B[i] = - C[i]
    C = np.linalg.solve(A,B)   
    F = np.vstack((F, C))

print(F)

#Gera o gráfico (arredonda as casas decimais)
# plt.plot(x, m)
# plt.show()