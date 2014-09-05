# define some hist.py run use case argument dictionaries below
# arguments can be ommited but must be supplied at hist.py command line
USE_CASES = {
    'STANDARD_1D' : {
        'dim' : 1,
        'bin_width' : 1.0,
        'debias' : True,
        'k' : 0.003,
        'fe' : True,
        'rtmin' : -30,
        'rtmax' : 30,
        'trim' : True,
        },

    'STANDARD_2D' : {
        'dim' : 2,
        'bin_width' : [1.0, 1.0],
        'debias' : True,
        'k' : [0.003, 0.003],
        'fe' : True,
        'rtmin' : [-30, -30],
        'rtmax' : [30, 30],
        'trim' : True,
        },
    }

DEFAULT = USE_CASES['STANDARD_1D']

def arg_merge(override_dict, default_dict):
    merged_dict = {}
    for key in default_dict.keys() + override_dict.keys():
        if key in override_dict:
            merged_dict[key] = override_dict[key]
        else:
            merged_dict[key] = default_dict[key]

    return merged_dict

def arg_typer(in_arg_dict):
    arg_dict = in_arg_dict
    
    for key in arg_dict:
        arg = arg_dict[key]

        if key in ('dim',):
            arg_dict[key] = int(arg)
    
        elif key in ('bin_width', 'k', 'rtmin', 'rtmax'):
            arg_dict[key] = arg.translate(None,'{[(])}').split(',')
            for i, el in enumerate(arg_dict[key]):
                arg_dict[key][i] = float(arg_dict[key][i])

        elif key in ('debias', 'fe', 'trim'):
            arg_dict[key] = bool(arg)

        else:
            raise ValueError('Unknown arguments')

    return arg_dict

def arg_list_to_dict(arg_list):
    if len(arg_list) % 2 != 0:
        raise ValueError('Keyword/Arguments must be paired')
    i = 0
    arg_dict = {}
    while i < len(arg_list):
        arg_dict[arg_list[i]] = arg_list[i+1]
        i += 2

    return arg_typer(arg_dict)
        
        
    
    
    
