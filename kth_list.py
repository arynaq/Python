import numpy as np

def kth_list(somelist,k):
	return somelist[::k], [x[i] for i in range(len(x)) if i%k]
	
x = range(10)
print kth_list(x,5)