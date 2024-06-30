#Importa as libs
import numpy as np 
import matplotlib.pyplot as plt 

def vortex_panel_function(Xb, Yb, c, angle_of_attack):
    Mp = len(Xb)  # Número total de pontos de dados no aerofólio.
    M = Mp - 1  # Número de painéis formados pelos pontos.
    angle_of_attack = np.deg2rad(angle_of_attack)  # Angulo de ataque em rad

    # Cria arrays de para cada paineis, comprimentos dos painéis, ângulos e etc
    x = np.zeros(M)
    y = np.zeros(M)
    s = np.zeros(M)
    theta = np.zeros(M)
    sine = np.zeros(M)
    cosine = np.zeros(M)
    RHS = np.zeros(M)
    
    # Calculo das variáveis criadas acima para cada painel
    for i in range(M):
        ip = i + 1  # Índice do ponto final do painel.
        x[i] = 0.5 * (Xb[i] + Xb[ip])  # Coordenada x do ponto de controle (meio do painel).
        y[i] = 0.5 * (Yb[i] + Yb[ip])  # Coordenada y do ponto de controle (meio do painel).
        s[i] = np.sqrt((Xb[ip] - Xb[i])**2 + (Yb[ip] - Yb[i])**2)  # Comprimento do painel.
        theta[i] = np.arctan2((Yb[ip] - Yb[i]), (Xb[ip] - Xb[i]))  # Ângulo do painel.
        sine[i] = np.sin(theta[i])  # Seno do ângulo do painel.
        cosine[i] = np.cos(theta[i])  # Cosseno do ângulo do painel.
        RHS[i] = np.sin(theta[i] - angle_of_attack)  # Termo do lado direito da equação

    # Inicializa matrizes de influência para os coeficientes normais e tangenciais.
    CN1 = np.zeros((M, M))
    CN2 = np.zeros((M, M))
    CT1 = np.zeros((M, M))
    CT2 = np.zeros((M, M))
    
    # Calcula os coeficientes de influência geométrica.
    for i in range(M):
        for j in range(M):
            if i == j:
                # Tratamento especial para o painel consigo mesmo (i == j)
                CN1[i, j] = -1.0
                CN2[i, j] = 1.0
                CT1[i, j] = 0.5 * np.pi
                CT2[i, j] = 0.5 * np.pi
            else:
                # Cálculo dos termos geométricos
                A = -(x[i] - Xb[j]) * cosine[j] - (y[i] - Yb[j]) * sine[j]
                B = (x[i] - Xb[j])**2 + (y[i] - Yb[j])**2
                C = np.sin(theta[i] - theta[j])
                D = np.cos(theta[i] - theta[j])
                E = (x[i] - Xb[j]) * sine[j] - (y[i] - Yb[j]) * cosine[j]
                F = np.log(1.0 + s[j] * (s[j] + 2 * A) / B)
                G = np.arctan2(E * s[j], B + A * s[j])
                P = (x[i] - Xb[j]) * np.sin(theta[i] - 2 * theta[j]) + (y[i] - Yb[j]) * np.cos(theta[i] - 2 * theta[j])
                Q = (x[i] - Xb[j]) * np.cos(theta[i] - 2 * theta[j]) - (y[i] - Yb[j]) * np.sin(theta[i] - 2 * theta[j])

                # Coeficientes de influência normal
                CN2[i, j] = D + 0.5 * Q * F / s[j] - (A * C + D * E) * G / s[j]
                CN1[i, j] = 0.5 * D * F + C * G - CN2[i, j]

                # Coeficientes de influência tangencial
                CT2[i, j] = C + 0.5 * P * F / s[j] + (A * D - C * E) * G / s[j]
                CT1[i, j] = 0.5 * C * F - D * G - CT2[i, j]

    # Inicializa a matriz de coeficientes de influência.
    AN = np.zeros((Mp, Mp))
    AT = np.zeros((M, Mp))
    for i in range(M):
        AN[i, 0] = CN1[i, 0]
        AN[i, Mp - 1] = CN2[i, M - 1]
        AT[i, 0] = CT1[i, 0]
        AT[i, Mp - 1] = CT2[i, M - 1]
        for j in range(1, M):
            AN[i, j] = CN1[i, j] + CN2[i, j - 1]
            AT[i, j] = CT1[i, j] + CT2[i, j - 1]

    # Ajusta a última linha da matriz de coeficientes de influência.
    AN[Mp - 1, 0] = 1.0
    AN[Mp - 1, Mp - 1] = 1.0
    for j in range(1, M):
        AN[Mp - 1, j] = 0.0
    RHS = np.append(RHS, 0.0)  # Adiciona um termo de condição de contorno à direita.

    # Resolve o sistema linear para obter a circulação nos painéis.
    Gama = np.linalg.solve(AN, RHS)
    V = np.zeros(M)  # Inicializa o array de velocidades tangenciais.
    CP = np.zeros(M)  # Inicializa o array de coeficientes de pressão.
    for i in range(M):
        V[i] = np.cos(theta[i] - angle_of_attack)
        for j in range(Mp):
            V[i] += AT[i, j] * Gama[j]
        CP[i] = 1.0 - V[i]**2  # Calcula o coeficiente de pressão.

    # Calcula a circulação total e os coeficientes de força.
    Circ = np.sum(V * s)
    cy = -np.sum(CP * s * cosine)
    cx = -np.sum(CP * s * -sine)
    cy /= c
    cx /= c
    c_d = cx * np.cos(angle_of_attack) + cy * np.sin(angle_of_attack)  # Coeficiente de arrasto.
    c_l = -cx * np.sin(angle_of_attack) + cy * np.cos(angle_of_attack)  # Coeficiente de sustentação.

    Gama = np.delete(Gama, (Mp - 1))  # Remove um valor específico de Gama para evitar singularidades.

    # Plota o coeficiente de pressão (C_p) ao longo da corda do aerofólio.
    plt.figure()
    plt.plot(x / c, CP)
    plt.title(f'-C_p vs x/c for α = {np.rad2deg(angle_of_attack):.2f}')
    plt.xlabel('x/c')
    plt.ylabel('C_p')
    plt.grid(True)
    plt.show()

    return c_l, c_d, V, s  # Retorna os coeficientes de sustentação e arrasto, a velocidade e o comprimento dos painéis.
