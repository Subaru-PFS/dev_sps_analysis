##
## misceallaneous or utils function
##


# str2bool : see https://cmsdk.com/python/parsing-boolean-values-with-argparse.html
# useful for argparse 
def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
    

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

# def find_nearest(X, value):
#    return X[np.unravel_index(np.argmin(np.abs(X - value)), X.shape)]
