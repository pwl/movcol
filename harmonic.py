import matplotlib.pylab as plt
import os
import numpy as np
import math as m

def plotn(f):
    data = np.loadtxt(f, dtype=np.dtype([('time',float),('dx',float)]))
    stime = data[-1]['time']
    p, = plt.plot(-np.log(abs(stime-data['time'])),data['dx']*np.sqrt(stime-data['time']))
    return(p)

def main():
    files = sorted([os.path.join('data/d007/k001',f,'at0.dat') for f in os.listdir('data/d007/k001')])
    # data = [np.loadtxt(f, dtype=np.dtype([('time',float),('dx',float)])) for f in files]

    # for f in files:
    #     plotn(f)

    plots = [plotn(f) for f in files ]

    plt.legend(plots,files)

    plt.ylim([0,20])
    plt.xlim([25,35])

    # plotn(data[-1],stime)

    plt.show()


if __name__ == "__main__":
    main()
