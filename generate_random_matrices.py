import random

def build_rand_mapping(filename,binSeed,nBins,N):
    f = open(filename,'w')
    random.seed(binSeed)
    binExpansion = 3
    # bins = [[] for i in range(nBins)]
    f.write('{\n')
    for fileIdx in range(N):
        bins = random.sample(xrange(nBins),binExpansion)
        f.write('{')
        for bin in bins[:-1]:
            f.write(hex(bin) + ", ")
        if fileIdx == N-1:
            f.write(hex(bins[-1]) + " } }")
        elif fileIdx % 4 == 3:
            f.write(hex(bins[-1]) + " },\n")
        else:
            f.write(hex(bins[-1]) + " }, ")
    f.close()
    
    
if __name__=='__main__':
    filename = 'test'
    binSeed = 1
    nBins = 11
    N = 10000
    
    build_rand_mapping(filename,binSeed,nBins,N)