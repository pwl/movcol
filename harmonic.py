import matplotlib.pylab as plt
import os
import numpy as np
import math as m
import re

def plotn(f):
    data = np.loadtxt(f, dtype=np.dtype([('time',float),('dx',float)]))
    stime = data[-1]['time']
    n = int(re.search('n([0-9]+)',f).group(1))
    p, = plt.plot(
        -np.log(abs(stime-data['time'])),
        data['dx']*np.sqrt(stime-data['time']),
        label = n)
    return(p)

def main():
    path='data/d027/k003'
    at0 = sorted([os.path.join(path,f, 'at0.dat') for f in os.listdir(path)])

    plots = [plotn(f) for f in at0 ]

    plt.legend(loc=2)

    # plt.ylim([0,10])
    # plt.xlim([0,35])

    plt.show()


if __name__ == "__main__":
    main()
