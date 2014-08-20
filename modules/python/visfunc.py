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
        plt.title('Bias factors vs Order Parameter Values')
        plt.xlabel('Order Parameter')
        plt.ylabel('Bias factor')
        plt.plot(self.x,self.y)


def runplot(numwindows,windsep,k,firstwindow):
    visfunc = VisFunc()
    visfunc.k = k
    visfunc.domain = [float(firstwindow),float(windsep)*float(numwindows) + float(firstwindow),float(windsep)/100.0]
    
    for n in range(int(numwindows)):
        visfunc.centre = float(firstwindow) + float(windsep)*n
        visfunc.evalfunc()
        visfunc.plot()

def run(numwindows,windsep,k,firstwindow):
    runplot(numwindows,windsep,k, firstwindow)
    plt.show()

    
def run_png(numwindows,windsep,k,firstwindow):
    runplot(numwindows,windsep,k,firstwindow)
    plt.savefig('biasfactors.png', bbox_inches='tight')

if __name__ == '__main__':
    #USAGE: visfunc.py [number of windows] [window separation] [k] [first window]
    run(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
    
