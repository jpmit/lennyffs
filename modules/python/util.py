# util.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Functions to aid with the utilities in the util/ directory,
e.g. len_sk, len_neighbour, etc.
"""

from ffsfunctions import getboxdims, getpickparams

def boxdims_xyz(xyzfname):
    """Return lboxx, lboxy, lboxz from XYZ file fname.

    If dims are not in the XYZ file, try getting dims from params.pkl.
    If neither of these methods produce the box dims, raise UtilError.
    """

    lboxx, lboxy, lboxz = None, None, None
    try:
        params = getpickparams()
    except IOError:
        params = {}
    else:
        lboxx = params['lboxx']
        lboxy = params['lboxy']
        lboxz = params['lboxz']

    boxdims = getboxdims(xyzfname)
    if boxdims:
        # note we override any dims we got from the params.pkl file
        lboxx = boxdims[0]
        lboxy = boxdims[1]
        lboxz = boxdims[2]

    if lboxx is None:
        raise UtilError, 'box dims not found in XYZ file or params.pkl'

    else:
        return lboxx, lboxy, lboxz
    
def getndict(maxn2):
    """Used for len_sk.

     Return dictionary that holds value of nx**2 + ny**2 + nz**2 as
     key and list of lists as value.
     """
    combod =  {k : [] for k in range(maxn2 + 1)}
    
    for nx in range(maxn2 + 1):
        nx2 = nx**2
        for ny in range(maxn2 + 1):
            ny2 = ny**2
            nint2 = nx2 + ny2
            if (nint2 > maxn2):
                # next nx value
                break
            for nz in range(maxn2 + 1):
                n2 = nint2 + nz**2
                if (n2 > maxn2):
                    # next ny value
                    break
                # add all combinations to dictionary
                combod[n2].append([nx, ny, nz])
                if (nx > 0):
                    if (ny > 0):
                        if (nz > 0): # nx, ny, nz all > 0
                            combod[n2].append([nx, -ny, nz])
                            combod[n2].append([nx, ny, -nz])
                            combod[n2].append([nx, -ny, -nz])
                        else: # nx, ny only > 0
                            combod[n2].append([nx, -ny, nz])
                    elif (nz > 0): # nx and nz only > 0
                        combod[n2].append([nx, ny, -nz])
                    else: # nx only > 0
                        pass
                else:
                    if (ny > 0):
                        if (nz > 0): # ny, nz only > 0
                            combod[n2].append([nx, -ny, nz])

    # remove 0 value
    del combod[0]
    
    # remove any values with no k vectors
    return {k : v for k, v in combod.items() if combod[k] != []}
