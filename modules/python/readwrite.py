# readwrite.py
# 10th September 2010
# James Mithen

"""
Functions for reading / writing files
FUNCTIONS:
rxyz - read file in XYZ format
wxyz - write file in XYZ format
r2col - read 2 column file
r3col - read 3 column file
r4col - read 4 column file
r5col -read 5 column file
rncol -read n column file
w2col - write 2 column file
"""

import numpy as np

def rxyz(fname,retsymbols=False,splines=1):
    """Read a .xyz coordinate file."""
    fin = open(fname,'r')
    lines = fin.readlines()
    fin.close()
    npar = int(lines[0])
    positions = np.empty([npar,3])
    symbols = [None]*npar
    coordline = 1 + splines
    i = 0
    for line in lines[coordline:]:
        li = line.split()
        symbols[i] = li[0]
        positions[i,0] = float(li[1])
        positions[i,1] = float(li[2])
        positions[i,2] = float(li[3])
        i = i + 1
    if retsymbols:
        return positions,symbols
    else:
        return positions

def wxyz(fname,positions,symbols):
    """Write a .xyz coordinate file"""
    fout = open(fname,'w')
    npar = len(positions)
    fstr = '%d\n\n' %npar
    for i in range(npar):
        fstr = '%s%s %.8f %.8f %.8f\n' %(fstr,symbols[i],positions[i][0],
                                         positions[i][1],positions[i][2])
    fout.write(fstr)
    fout.close()
    return

def r1col(fname):
    """Read 1 column file of numeric data and return as numpy array"""
    fin = open(fname,'r')
    lines = fin.readlines()
    nlines = len(lines)
    fin.close()
    col1 = np.empty(nlines)
    i = 0 # count line num
    # invariant: we have stored i lines in col1
    for line in lines:
        if '#' in line:
            continue
        col1[i] = float(line)
        i += 1

    return col1[:i]
    

def r2col(fname,sep=None,**kwargs):
    """Read 2 column file of numeric data and return as numpy arrays"""
    sskip = kwargs.get('eskip',0) # skip lines at start
    eskip = kwargs.get('eskip',0) # skip lines at end
    fin = open(fname,'r')
    lines = fin.readlines()
    nlines = len(lines)
    fin.close()
    col1 = np.empty(nlines - sskip - eskip)
    col2 = np.empty(nlines - sskip - eskip)
    i = 0 # count line num
    # invariant: we have stored i lines in col1 and col2
    for line in lines[sskip:nlines-eskip]:
        if '#' in line:
            continue
        splin = line.split(sep)
        col1[i] = float(splin[0])
        try:
            col2[i] = float(splin[1])
            i += 1
        except:
            pass

    return col1[:i], col2[:i]

def r3col(fname,sep=None):
    """Read 3 column file of numeric data and return as numpy arrays"""
    fin = open(fname,'r')
    lines = fin.readlines()
    fin.close()
    col1, col2, col3 = [],[],[]
    for line in lines:
        if '#' in line:
            continue
        splin = line.split(sep)
        if len(splin) == 3:
            col1.append(float(splin[0]))
            col2.append(float(splin[1]))
            col3.append(float(splin[2]))
    return np.array(col1), np.array(col2), np.array(col3)

def r4col(fname,sep=None):
    """Read 4 column file of numeric data and return as numpy arrays"""
    fin = open(fname,'r')
    lines = fin.readlines()
    fin.close()
    col1, col2, col3, col4 = [],[],[],[]
    for line in lines:
        if '#' in line:
            continue
        splin = line.split(sep)
        col1.append(float(splin[0]))
        col2.append(float(splin[1]))
        col3.append(float(splin[2]))
        col4.append(float(splin[3]))
    return np.array(col1), np.array(col2), np.array(col3), np.array(col4)

def r5col(fname,sep=None):
    """Read 4 column file of numeric data and return as numpy arrays"""
    fin = open(fname,'r')
    lines = fin.readlines()
    fin.close()
    col1, col2, col3, col4, col5 = [],[],[],[],[]
    for line in lines:
        if '#' in line:
            continue
        splin = line.split(sep)
        col1.append(float(splin[0]))
        col2.append(float(splin[1]))
        col3.append(float(splin[2]))
        col4.append(float(splin[3]))
        col5.append(float(splin[4]))
    return np.array(col1), np.array(col2), np.array(col3), np.array(col4), np.array(col5)

def rncol(fname,n,sep=None):
    """Read n column file of numeric data and return as list of numpy arrays"""
    fin = open(fname,'r')
    lines = fin.readlines()
    fin.close()
    res = [[]]
    for i in range(n-1):
        res.append([])
    for line in lines:
        if '#' in line:
            continue  
        splin = line.split(sep)
        for i in range(n):
            res[i].append(float(splin[i]))
    return [np.array(col) for col in res]
 
def w2col(fname,col1,col2,headers=[]):
    """Write 2 column data of floats"""
    if headers:
        fstr = '# %s %s\n' %(headers[0],headers[1])
    else:
        fstr = ''
    for (i,j) in zip(col1,col2):
        fstr = '%s%.6f %.6f\n' %(fstr,i,j)
    outf = open(fname,'w')
    outf.write(fstr)
    outf.close()
    return
