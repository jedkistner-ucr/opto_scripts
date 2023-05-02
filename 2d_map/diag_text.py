import numpy as np

a = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]])
print(a)
print(np.diagonal(a))
print(np.diagonal(a, offset = 1))
print(np.diagonal(a, offset = -1))
a = np.rot90(a)
print(a)