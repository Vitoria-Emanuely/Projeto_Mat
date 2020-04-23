import numpy as np
import matplotlib.pyplot as plt

#Início
a = 0

#Fim
b = 1

#Números de pontos
n = 5

#Cria uma malha com os parâmetros (inicio, fim, n de "divisões" na malha)
x = np.linspace(a, b, n)

deltax = (b - a) / (n - 1)

#Cria uma matriz de zeros NXN
A = np.zeros((n,n))
B = np.zeros(n)

#Valores fixos
A[0, 0] = 2.0
A[0, 1] = -1.0
A[n - 1, n - 1] = 2.0
A[n - 1, n - 2] = -1.0

#Preenche as matrizes A e B
for i in np.arange(1, n - 1):
    A[i, i] = 2.0
    A[i, i + 1] = -1.0
    A[i, i - 1] = -1.0

    B[i] = 100 * deltax**2 * (x[i] - 1)**2

#Resolve a matriz linear
U = np.linalg.solve(A,B)
print(U)

#Gera o gráfico (arredonda as casas decimais)
plt.plot(x, U)
plt.show()