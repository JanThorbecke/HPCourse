import time
import numpy as np

start_time = time.time()

data = np.ones(shape=(1000, 1000), dtype=float)

for i in range(1000):
    for j in range(1000):
        data[i][j] *= 1.0000001
        data[i][j] *= 1.0000001
        data[i][j] *= 1.0000001
        data[i][j] *= 1.0000001
        data[i][j] *= 1.0000001

end_time = time.time()

print("Run time = {}".format(end_time - start_time))

data = np.ones(shape=(1000, 1000), dtype=float)

start_time = time.time()

for i in range(5):
    data *= 1.0000001

end_time = time.time()

print("Run time = {}".format(end_time - start_time))
