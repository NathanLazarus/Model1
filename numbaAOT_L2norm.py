from numba.pycc import CC
import numpy

cc = CC('l2norm_compiled')
# Uncomment the following line to print out the compilation steps
#cc.verbose = True

@cc.export('l2norm', 'f8(f8, f8)')
def l2norm(x, y):
    return ((x - y) ** 2) #.sum(-1)

if __name__ == "__main__":
    cc.compile()