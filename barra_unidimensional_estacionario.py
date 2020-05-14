import numpy as np
import matplotlib.pyplot as plt

#Início
a = 0.0

#Fim
b = 2.0

#Números de pontos
n = 5

#Cria uma malha com os parâmetros (inicio, fim, n de "divisões" na malha)
x = np.linspace(a, b, n)

deltax = (b - a) / (n - 1)

#Parâmetros da equação AK * (d²T(x)/dx²) - Ph * (T(x) - Tinf) = 0
A = 0.25
K = 3.8
T = 240.0
P = 2.5
h = 0.015
Tinf = 25.0

#Cria uma matriz de zeros NXN
X = np.zeros((n + 2, n + 2))
Y = np.zeros(n + 2)

#Valores fixos
#Contorno a esquerda
X[0, 0] = K / 2.0 * deltax
X[0, 2] = -K / (2.0 * deltax)
Y[0] = 127.324

#Contorno a direita
X[n + 1, n + 1] = - K / (2.0 * deltax)
X[n + 1, n - 1] = K / (2.0 * deltax)
Y[n + 1] = - h * Tinf

#Preenche as matrizes X e Y
for i in np.arange(1, n):
    X[i + 1, i] = (A * K) / deltax ** 2
    X[i, i] =  - ((2 * A * K / (deltax ** 2)) + P * h)
    X[i, i + 1] = (A * K) / deltax ** 2

    Y[i] = - P * h * Tinf

#Resolve a matriz linear
Z = np.linalg.solve(X,Y)

#Salva na matriz W a solução Z sem as entradas que não pertencem a malha (remove alguma casas decimais)
lista = []
for i, elem in enumerate(Z):
    if i != 0 and i != n + 1:
        lista.append(elem)

W = np.asarray(lista)
print(W)

#Gera o gráfico (arredonda as casas decimais)
plt.plot(x, W)
plt.show()
