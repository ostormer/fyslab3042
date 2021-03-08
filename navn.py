import numpy as np

t, x, y = np.genfromtxt("data/bane.csv", delimiter=",", skip_header=2, unpack=True)

print(t)
