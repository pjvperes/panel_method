import numpy as np

def test():
    data = np.loadtxt('airfoil.dat')
    Xb = data[:, 0]  # Extracts the X coordinates of the airfoil points.
    Mp = len(Xb)

    return Mp

if __name__ == "__main__":
    result = test()
    print(f"Number of points in the airfoil data: {result}")
