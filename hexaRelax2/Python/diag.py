import numpy as np
import matplotlib.pyplot as plt

t = np.array([])
e  = np.array([])
file = open('run1/diagnostic', 'r')
while True:
    line = file.readline()
    if not line:
        break
    line = np.float32(line.split())
    t  = np.append(t,  line[0])
    e  = np.append(e,  line[1])


t = np.array(t)
e = np.array(e)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(t / 3600, e - e[0])

plt.show()
