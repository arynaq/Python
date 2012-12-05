#isorthogonal?
import numpy as np
import matplotlib.pyplot as plt

a = np.random.randint(-5,5,2)
b = np.random.randint(-5,5,2)
a1= np.dot(a,b)/np.sqrt(b.dot(b)) * b
a2= a- a1
print a1
print a
print b
print a2

plt.axis([-5,5,-5,5])
plt.arrow(0,0,a[0], a[1],color="b",label="a")
plt.arrow(0,0,b[0], b[1],color="g",label="b")
plt.arrow(0,0,a1[0], a1[1],color="r",label="a1")
plt.arrow(0,0,a2[0], a2[1],color="c",label="a2")
plt.show()