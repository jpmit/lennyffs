import initsim
import pickle
import os
import funcselector
import visfunc
from sys import argv

class Windowmaker(object):
    """
    Produces window parameters to run with umbrella.py
    NOTE: the bash script produced is configured to run
          on the computer cluster
    """

    def __init__(self):

        self.params = initsim.getparams()
        self.windowcentres = [self.params['firstwindow'] + n*self.params['windowsep'] for n in range(self.params['numwindows'])]
        self.dir = os.path.abspath(__file__)
        funcman = funcselector.FuncSelector(self.params)
        self.writexyz = funcman.WriteXyzFunc()
        self.params['simulation'] = 'restart'
        # if queue_name_form supplied, the queue name will be this with the
        # window index appended
        self.queue_name_form = 'unnamed'
        # add making master-pickles here

    def prep_windows(self):
        """Produces initial positions and parameter files"""
        bashscript = open('wumbash.sh','w')
        for wcentre in self.windowcentres:
            self.params['restartfile'] = 'initialpositions{0}.xyz'.format(wcentre)
            self.params['umb_centre'] = wcentre
            self.params['nparseed'] = initsim.deduce_seed_size(self.params)
            pickle.dump(self.params, open('params{0}.pkl'.format(wcentre),'w'))

            # name option added to qsub command in bash script if queue_name_form supplied
            if self.queue_name_form == 'unnamed':
                namepart = ''
            else:
                namepart = ' -N \"{0}-win{1}\"'.format(self.queue_name_form,wcentre)
                
            bashscript.write(str(os.path.split(self.dir)[0].join(['qsub {0} -cwd -b y python '.format(namepart),'/umbrella.py -w {0}\n'.format(wcentre)])))
            positions = initsim.initpositionsseed(self.params)
            self.writexyz('initialpositions{0}.xyz'.format(wcentre), positions, self.params)
        bashscript.close()

    def visualise_potentials(self):
        visfunc.run_png(self.params['numwindows'],
                        self.params['windowsep'],
                        self.params['k'],
                        self.params['firstwindow'])
    
            
            
if __name__ == '__main__':
    windowmaker = Windowmaker()
    if argv[1:2] != []:
        windowmaker.queue_name_form = argv[1]
    windowmaker.prep_windows()
    #windowmaker.visualise_potentials()
    
