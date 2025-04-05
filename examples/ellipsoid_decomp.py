import numpy as np
import matplotlib.pyplot as plt


# I want to test what an ellipsoid's axis are given some E

E = [[1, 0.8], [0.8, 0.9]]

# Perform Cholesky decomposition
C = np.linalg.cholesky(E)
C_inv = np.linalg.inv(C)

print("C_inv: ", C_inv)

for i in range(len(C_inv)):
    print("axis: ", C_inv[:,i])

