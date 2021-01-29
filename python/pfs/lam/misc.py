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