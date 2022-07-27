import numpy as np
from matplotlib import pyplot as plt
import ctypes

def fft(x):
    if (type(x) != np.ndarray):
        x = np.array(x)
    F = np.zeros(x.shape, dtype=np.complex128)
    if(len(x) == 1):
        return x
    even = fft(x[::2])
    odd = fft(x[1::2])
    
    F[0] = x.sum()
    F[len(F)//2] = x[::2].sum()-x[1::2].sum()
    for k in range(1, len(F)//2):
        omega = np.exp(-2j*np.pi*k/len(x))
        F[k] = even[k] + omega*odd[k]
        F[len(F)-k] = F[k].conj()
    return F

    

if __name__ == "__main__":
    print("Hello World!")
    fft(np.array([3, 4, 5, 6]))
    t = np.linspace(1, 10, 1024*32)
    y = 3*np.sin(2*np.pi*5*t)
    Fy = fft(y)
