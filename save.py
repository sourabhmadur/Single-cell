from neuron import h,gui
import numpy as np
from matplotlib import pyplot as plt
h.tstop = 10
a = h.Section()
b = h.Section()


vclamp = h.SEClamp(a(0.5))
vclamp.dur1 = 5
vclamp.amp1=-50.0
vclamp.rs=.00001

v1 = h.Vector()
v2 = h.Vector()
t= h.Vector()

t.record(h._ref_t)
v1.record(a(0.5)._ref_v)
v2.record(b(0.5)._ref_v)


h.run()

v1 = v1.to_python()
v2 = v2.to_python()
t = t.to_python()
v=[]
v.append(v1)
v.append(v2)
np.savetxt('test1.txt', v, fmt='%d')
np.savetxt('test2.txt', t, fmt='%d')

y = np.loadtxt('test1.txt', dtype=float)
x = np.loadtxt('test2.txt', dtype=float)

plt.plot(x,y[1])
plt.show()
# Out[4]: array([ True,  True,  True,  True], dtype=bool)
