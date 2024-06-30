# Import de bibliotecas
import numpy as np
import matplotlib.pyplot as plt 

# Import do código do método dos painéis
from vortex_panel_function import vortex_panel_function 

def main():

    Vinf = 100 
    c = 1.0  # Corda do aerofólio

    # Armazena os dados do arquivo .dat (aerofólio)
    data = np.loadtxt('airfoil.dat')
    Xb = data[:, 0]  # Extrai as coordenadas X dos pontos do aerofólio.
    Yb = np.flip(data[:, 1])  # Extrai e inverte as coordenadas Y dos pontos do aerofólio.

    # Plota o aerofólio
    plt.figure()
    plt.plot(Xb, Yb)
    plt.xlabel('Xb')
    plt.ylabel('Yb')
    plt.title('Airfoil')
    plt.axis('equal')
    plt.show()

    # Ângulos de ataque a serem analisados
    angle_of_attack = [0, 4, 8, 12, 16]
    CL = [] 
    CD = []
    
    # Itera sobre cada ângulo de ataque para calcular os coeficientes de sustentação e arrasto.
    for alpha in angle_of_attack:
        print(f'At {alpha} deg:')
        c_l, c_d, _, _ = vortex_panel_function(Xb, Yb, c, alpha)  # Calcula os coeficientes usando a função vortex_panel_function (disponivel entre outro codigo)
        CL.append(c_l)  # Armazena o coeficiente de sustentação calculado.
        CD.append(c_d)  # Armazena o coeficiente de arrasto calculado.
        print(f'C_L = {c_l:.6f}')
        print(f'C_D = {c_d:.6f}\n') 

    # Dados do XFoil
    alpha_given = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0]
    CL_given = [0.0000, 0.1068, 0.2130, 0.3193, 0.4390, 0.5786, 0.7190, 0.8239, 0.9070, 0.9905, 1.0749, 1.1547, 1.2278, 1.2739, 1.3052, 1.3001, 1.2851]
    CD_given = [0.00559, 0.00570, 0.00609, 0.00677, 0.00778, 0.00909, 0.01044, 0.01170, 0.01303, 0.01460, 0.01650, 0.01898, 0.02169, 0.02538, 0.03065, 0.04155, 0.05695]

    # Plota o Cl em função do ângulo de ataque (α).
    plt.figure()
    plt.plot(angle_of_attack, CL, 'b-', linewidth=2, label='C_L (Método dos Painéis)')
    plt.plot(alpha_given, CL_given, 'r--', linewidth=2, label='C_L (XFoil)')
    plt.xlabel('α (graus)')
    plt.ylabel('C_L')
    plt.title('C_L vs α')
    plt.legend()
    plt.grid(True)
    plt.show()

    # Plota o Cd em função do ângulo de ataque (α).
    plt.figure()
    plt.plot(angle_of_attack, CD, 'b-', linewidth=2, label='C_D (Método dos Painéis)')
    plt.plot(alpha_given, CD_given, 'r--', linewidth=2, label='C_D (XFoil)')
    plt.xlabel('α (graus)')
    plt.ylabel('C_D')
    plt.title('C_D vs α')
    plt.legend()
    plt.grid(True)
    plt.show()

# Executa a função principal se o script for executado diretamente.
if __name__ == "__main__":
    main()
