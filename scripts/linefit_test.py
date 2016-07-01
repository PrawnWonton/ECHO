import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0,10,23)
y = np.random.rand(23)

m, b = np.polyfit(x,y,1)

plt.plot(x,y)
plt.plot(x, m*x+b, 'r')
plt.show()

print('line of fitty-ness: y = {:.3f}x + {:.4f}'.format(m,b))
