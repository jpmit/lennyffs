import pylab as pl
import glob
import numpy as np

class HistMaker():
    def __init__(self, infile='surfhist.in', outfile='surfhist.out',
                 bin_width='NotSet', dimensions=1):
        # filenames stored as strings
        self.infile = infile
        self.outfile = outfile
        
        # error raised if dimension is not 1 or 2
        if dimensions in (1, 2):
            self.dim = dimensions
        else:
            raise ValueError('dimensions must be 1 or 2')

        # bin_width set to list
        if bin_width=='NotSet':
            bin_width = [1.0]*self.dim
        elif type(bin_width) != list:
            bin_width = [bin_width]
        if len(bin_width) != self.dim:
            raise ValueError('bin_width must be length dimension')
        else:
            self.bin_width = bin_width

        self.numbins = [1]*self.dim

        # no trimming implemented by infinite limits for syntax purposes
        self.tmin = [float('-inf')]*self.dim
        self.tmax = [float('inf')]*self.dim

        # loads 1 or 2 columns into op_array depending on self.dim
        self.op_array = pl.loadtxt(infile)[:,[n for n in range(1,self.dim+1)]]

    # set_lims will take 2 lists (or ints for 1D) defining the trim limits
    def set_lims(self, min_list, max_list):
        if type(min_list) != list:
            min_list = [min_list]

        if type(max_list) != list:
            max_list = [max_list]

        if len(min_list) != self.dim or len(max_list) != self.dim:
            raise ValueError('dimensions do not match')

        self.tmin = min_list
        self.tmax = max_list

    def makebins(self):
        # bin_maxes and bin_mins decide whether to set the max and min bins
        # to match the max/min data points or the trim limits - they are lists
        # of length self.dim
        bin_maxes = [0]*self.dim
        for i, m in enumerate(bin_maxes):
            bin_maxes[i] = min(self.tmax[0], pl.amax(self.op_array[:,i]))
        bin_mins = [0]*self.dim
        for i, m in enumerate(bin_mins):
            bin_mins[i] = max(self.tmin[0], pl.amin(self.op_array[:,i]))

        # number of bins in each dimension determined from bin_maxes/mins and
        # bin_width
        for i, n in enumerate(self.numbins):
            self.numbins[i] = int(pl.floor(bin_maxes[i]\
                                           /self.bin_width[i]) - \
                                  pl.floor(bin_mins[i]\
                                           /self.bin_width[i])) + 1

        self.bin_starts, self.bin_middles = [0]*self.dim, [0]*self.dim

        # bin_starts will contain a list for each dimension of the starting
        # value for each bin. Determined from max and min data, bin_width and
        # numbins.
        # MAYBE SHOULD USE bin_maxes, bin_mins INSTEAD OF MAX, MIN DATA
        for i in range(self.dim):
            self.bin_starts[i] = [(int((pl.amin(self.op_array[:,i])\
                                        /self.bin_width[i]))*\
                                 self.bin_width[i]) + n*self.bin_width[i] \
                                  for n in range(self.numbins[i])]

        # empty bins initialised - zeroed array of shape determined by numbins
        # in each dimension
        self.bin_counts = np.zeros(self.numbins)

        # bin_middles determined from bin_starts and bin_width
        # will contain a list for each dimension corresponding
        # EXACTLY to bin_starts
        for i,s in enumerate(self.bin_starts):
            self.bin_middles[i] = [bin + (self.bin_width[i]/2.0) \
                                   for bin in self.bin_starts[i]]

    def bin_data(self):
        for row in self.op_array:
            j = []
            # j will contain 1 or 2 indices depending on dimensionality of data
            for i, s in enumerate(self.bin_starts):
                j.append(int((row[i] - self.bin_starts[i][0])\
                              /self.bin_width[i]))

            # data checked against trim limits in all dimensions
            tcheck = True
            for i, t in enumerate(self.tmin):
                if row[i] < self.tmin[i] or row[i] > self.tmax[i]:
                    tcheck = False

            # data added to bin_counts if within limits
            # currently 1D and 2D cases handled separately
            if tcheck and self.dim == 1:
                self.bin_counts[j[0]] += 1

            if tcheck and self.dim == 2:
                self.bin_counts[j[0]][j[1]] += 1

    def debias(self, k, bcen):
        # k and bcen set to be 1D lists if integers        
        if type(k) != list:
            k = [k]
        if type(bcen) != list:
            bcen = [bcen]

        # dimensionality of k and bcen checked
        if len(k) != self.dim or len(bcen) != self.dim:
            raise ValueError('dimensions do not match')

        # quadratic-exponential bias factor calculated using data point, bcen
        # and k for each dimension - data then debiased in place
        for i, line in enumerate(self.surfhist):
            if line == 'LINEBREAK':
                pass
            else:
                quads = 0
                for j, kval in enumerate(k):

                    # bin_middle offset accounted for here - need to find a
                    # better solution
                    
                    quads += k[j]*((bcen[j] - (line[0] - (self.bin_width[j]\
                                                          /2.0)))**2)

                w = 0.5*quads
                W = pl.exp(-w)
                self.surfhist[i][self.dim] /= W

    def free_energy(self):
        # converts all counts into free energies(arbitrary constant offset)
        i = 0
        while i != len(self.surfhist):
            if self.surfhist[i] == 'LINEBREAK':
                # skips if LINEBREAK
                i += 1
            elif self.surfhist[i][self.dim] == 0.0:
                # removes point if count = 0
                self.surfhist.pop(i)
            else:
                # converts count to free energy in place
                self.surfhist[i][self.dim] = -1*pl.log(\
                    self.surfhist[i][self.dim])
                i += 1

    def make_hist(self):
        # converts bin_counts to histogram list of the form
        # [op1, fe_or_count] for 1D or [op1, op2, fe_or_count] for 2D
        # currently 1 and 2D cases handled separately
        self.surfhist = []
        if self.dim == 2:
            for i, x in enumerate(self.bin_middles[0]):
                for j, y in enumerate(self.bin_middles[1]):
                    self.surfhist.append([x, y, self.bin_counts[i][j]])
                self.surfhist.append('LINEBREAK')

        elif self.dim == 1:
            for i, x in enumerate(self.bin_middles[0]):
                self.surfhist.append([x, self.bin_counts[i]])


    def write_to_file(self):
        # just writes histogram list to file, delimited by space
        out = open(self.outfile, 'w')
        for line in self.surfhist:
            if line == 'LINEBREAK':
                out.write('\n')
            else:
                to_write = '{0} '.format(line[0])
                for val in line[1:]:
                    to_write += str(val) + ' '
                to_write += '\n'
                out.write(to_write)
        out.close()

# THIS IS THE DRIVER                    
def run(dim=1, bin_width='NotSet', debias=False, k='NotSet', fe=False,\
        trim=False, rtmin='NotSet', rtmax='NotSet'):
    """
Reads every opval* file (excluding opvalequil*) and processes using the
HistMaker class in this module.

Args
-------
- dim:        Dimension of data. Can be 1 or 2. Defaults to 1.

- bin_width:  Width of the histogram bins. Can be int or list of length 1 or 2.
              Length must match dim. Defaults to 1.0 for every dimension.
              
- debias:     Boolean. Switches debiasing on or off. Defaults to False.

- k:          k value(s) for debiasing. Can be float or list of length 1 or 2.
              Length must match dim. Defaults to 1.0 for every dimension. Only
              used if debias is True.

- fe:         Boolean. Switches Free Energy calculation on or off. Defaults to
              False.

- trim:       Boolean. Switches trimming on or off. Defaults to False.

- rtmin:      Relative trimming minimum. Relative to window centre. Can be
              float or list of length 1 or 2. Length must match dim. Defaults
              to negative infinity. Only used if trim is True.
              
- rtmax:      Relative trimming maximum. Relative to window centre. Can be
              float or list of length 1 or 2. Length must match dim. Defaults
              to positve infinity. Only used if trim is True.
    """

    if dim in (1, 2):
        pass
    else:
        raise ValueError('dim must be 1 or 2')

    # type checking and dimension checking blocks
    if bin_width=='NotSet':
        bin_width = [1.0]*dim
    elif type(bin_width) != list:
        bin_width = [bin_width]
    if len(bin_width) != dim:
        raise ValueError('bin_width must be list of length dim')

    if rtmin=='NotSet':
        rtmin = [float('-inf')]*dim
    elif type(rtmin) != list:
        rtmin = [rtmin]
    if len(rtmin) != dim:
        raise ValueError('rtmin must be list of length dim')

    if rtmax=='NotSet':
        rtmax = [float('inf')]*dim
    elif type(rtmax) != list:
        rtmax = [rtmax]
    if len(rtmax) != dim:
        raise ValueError('rtmax must be list of length dim')

    if k=='NotSet':
        k = [1.0]*dim
    elif type(k) != list:
        k  = [k]
    if len(bin_width) != dim:
        raise ValueError('k must be list of length dim')

    # cycles through every opval* but not opvalequil* file in directory
    for files in glob.iglob('opval[!equil]*'):
        hm = HistMaker(files, 'fhist' + files[5:], bin_width, dimensions=dim)
        hm.makebins()
        if debias or trim:
            # bias centre bcen set from opval file name
            try:
                bcen = [float(val) for val in files[5:-4].split('_')]
            except ValueError:
                if files[5:-4].split('_') == ['']:
                    bcen = [0.0]*dim
                else:
                    raise                   
        if trim:
            # trim limits set from bcen and rtmin/rtmax values
            trimmin = [bcen[i] + rtmin[i] for i in range(dim)]
            trimmax = [bcen[i] + rtmax[i] for i in range(dim)]
            hm.set_lims(trimmin, trimmax)
        hm.bin_data()
        hm.make_hist()
        if debias:
            hm.debias(k, bcen)
        if fe:
            hm.free_energy()
        hm.write_to_file()
        # returns the class - useful for debugging
    return hm

#if __name__ == '__main__':
    #cl = run(dim=2, bin_width=[1.0, 1.0], debias=True, k=[0.003, 0.003], fe=True,
             #rtmin=[-30, -30], rtmax=[30, 30], trim=True)
