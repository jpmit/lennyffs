import initsim
import pickle
import os

class Windowmaker(object):
    """Produces window parameters to run with umbrella.py"""

    def __init__(self):

        self.params = initsim.getparams()
        self.windowcentres = [n*self.params['windowsep'] for n in range(self.params['numwindows'])]
        self.dir = os.path.abspath(__file__)
        #Add making master-pickles here

    def prep_windows(self):
        bashscript = open('wumbash.sh','w')
        for wcentre in self.windowcentres:
            self.params['N0'] = wcentre
            pickle.dump(self.params, open('params{0}.pkl'.format(wcentre),'w'))
            bashscript.write(str(os.path.split(self.dir)[0].join(['qsub -cwd -b y python ','/umbrella.py {0}\n'.format(wcentre)])))
        bashscript.close()
            
            
if __name__ == '__main__':
    windowmaker = Windowmaker()
    windowmaker.prep_windows()
    
