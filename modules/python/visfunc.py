import numpy as np
import matplotlib.pyplot as plt
import sys

class VisFunc:
    def __init__(self):
        self.func = self.defaultfunc
        #domain of form [max,min,step]
        self.domain = [0,1,0.1]
        self.display = plt.show()
        self.x = [0]
        self.y = [0]
        self.centre = 0.0
        self.k = 0.01

    def defaultfunc(self,xval):
        w =  0.5*float(self.k)*(float(xval) - float(self.centre))**2
        return np.exp(-w)

    def evalfunc(self):
        self.x = np.arange(self.domain[0],
                           self.domain[1],
                           self.domain[2])
        self.y = [self.func(xval) for xval in self.x]

    def plot(self):
        plt.plot(self.x,self.y)


def runplot(numwindows,windsep,k):
    visfunc = VisFunc()
    visfunc.k = k
    visfunc.domain = [0,float(windsep)*float(numwindows),1]
    for n in range(int(numwindows)):
        visfunc.centre = float(windsep)*n
        visfunc.evalfunc()
        visfunc.plot()

def run(numwindows,windsep,k):
    runplot(numwindows,windsep,k)
    plt.show()

    
def run_png(numwindows,windsep,k):
    runplot(numwindows,windsep,k)
    plt.savefig('biaspotentials.png', bbox_inches='tight')

if __name__ == '__main__':
    #USAGE: visfunc.py [number of windows] [window separation] [k]
    run(sys.argv[1],sys.argv[2],sys.argv[3])
    
