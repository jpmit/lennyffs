import initsim
import pickle
import os
import funcselector
import visfunc

class Windowmaker(object):
    """
    Produces window parameters to run with umbrella.py
    NOTE: the bash script produced is configured to run
          on the computer cluster
    """

    def __init__(self):

        self.params = initsim.getparams()
        self.windowcentres = [n*self.params['windowsep'] for n in range(self.params['numwindows'])]
        self.dir = os.path.abspath(__file__)
        funcman = funcselector.FuncSelector(self.params)
        self.writexyz = funcman.WriteXyzFunc()
        self.params['simulation'] = 'restart'
        #Add making master-pickles here

    def prep_windows(self):
        """Produces initial positions and parameter files"""
        bashscript = open('wumbash.sh','w')
        for wcentre in self.windowcentres:
            self.params['restartfile'] = 'initialpositions{0}.xyz'.format(wcentre)
            self.params['N0'] = wcentre
            self.params['nparseed'] = wcentre
            pickle.dump(self.params, open('params{0}.pkl'.format(wcentre),'w'))
            bashscript.write(str(os.path.split(self.dir)[0].join(['qsub -cwd -b y python ','/umbrella.py -w {0}\n'.format(wcentre)])))
            positions = initsim.initpositionsseed(self.params)
            self.writexyz('initialpositions{0}.xyz'.format(wcentre), positions, self.params)
        bashscript.close()

    def visualise_potentials(self):
        visfunc.run_png(self.params['numwindows'],
                        self.params['windowsep'],
                        self.params['k'])
    
            
            
if __name__ == '__main__':
    windowmaker = Windowmaker()
    windowmaker.prep_windows()
    windowmaker.visualise_potentials()
    
