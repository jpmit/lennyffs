import initsim
import pickle
import os
import funcselector
#VISFUNC NOT CURRENTLY WORKING
#import visfunc
from sys import argv, exit

class Windowmaker(object):
    """
    Produces window parameters to run with umbrella.py
    NOTE: the bash script produced is configured to run
          on the computer cluster
    """

    def __init__(self):

        self.params = initsim.getparams()
        
        if len(self.params['numwindows']) == 1:
            self.windowcentres = [[self.params['firstwindow'][0] + n*self.params['windowsep'][0]] \
                                  for n in range(self.params['numwindows'][0])]
        elif len(self.params['numwindows']) == 2:
            self.windowcentres = [[self.params['firstwindow'][0] + n*self.params['windowsep'][0], \
                                   self.params['firstwindow'][1] + m*self.params['windowsep'][1]] \
                                  for m in range(self.params['numwindows'][1]) \
                                  for n in range(self.params['numwindows'][0])]
        else:
            print " > 2 order parameters not supported"
            exit(0)
                      
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
            
            strwc = ''
            for c in wcentre:
                strwc += str(c) + '-'
            strwc = strwc[:-1]
        
            self.params['restartfile'] = 'initialpositions{0}.xyz'.format(strwc)
            self.params['umb_centre'] = wcentre
            self.params['nparseed'] = initsim.deduce_seed_size(self.params)
            pickle.dump(self.params, open('params{0}.pkl'.format(strwc),'w'))

            # name option added to qsub command in bash script if queue_name_form supplied
            if self.queue_name_form == 'unnamed':
                namepart = ''
            else:
                namepart = ' -N \"{0}-win{1}\"'.format(self.queue_name_form,strwc)
                
            bashscript.write(str(os.path.split(self.dir)[0].join(['qsub {0} -cwd -b y python '.format(namepart),'/umbrella.py -w {0}\n'.format(strwc)])))
            positions = initsim.initpositionsseed(self.params)
            self.writexyz('initialpositions{0}.xyz'.format(strwc), positions, self.params)
        bashscript.close()

#THE VISFUNC MODULE DOESN'T CURRENTLY WORK
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
    
