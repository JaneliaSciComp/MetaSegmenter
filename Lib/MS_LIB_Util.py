import numpy
import math

# ------------------------------------------------------------------------

def get_VI_counts(im1, im2):
    N1  = {}
    N2  = {}
    N12 = {}

    # Extract counts
    print "Creating list of labels1 ..."
    labels1 = numpy.sort(numpy.unique(im1))
    print " ... len(labels1)=", len(labels1), " labels1=", labels1
    print "Creating list of labels2 ..."
    labels2 = numpy.sort(numpy.unique(im2))
    print " ... len(labels2)=", len(labels2), " labels2=", labels2
    print "Computing counts N1 ..."
    for k1 in labels1:
        count = (im1 == k1).sum()
        if not k1 in N1.keys():
            N1[str(k1)]  = count
        else:
            N1[str(k1)] += count
    print "Computing counts N2 ..."
    for k2 in labels2:
        count = (im2 == k2).sum()
        if not k2 in N2.keys():
            N2[str(k2)]  = count
        else:
            N2[str(k2)] += count
    print "Computing counts N12 ..."
    for k1 in labels1:
        for k2 in labels2:
            k12 = str(k1) + "_" + str(k2)
            count = ((im1 == k1) & (im2 == k2)).sum()
            if not k12 in N12.keys():
                N12[k12]  = count
            else:
                N12[k12] += count
    print "len(labels1)=", len(labels1), " len(keys1)=", len(N1.keys()), \
         " len(labels2)=", len(labels2), " len(keys2)=", len(N2.keys())

    print "len(N1.values())=", len(N1.values()), " len( N2.values())=",  len( N2.values()), " len(N12.values())=", len(N12.values())

    return (N1, N2, N12)

# -------------------------------------------------------------------------

def counts2VI2D(N1, N2, N12):
    print "N1=", N1
    print "N2=", N2
    print "N12=", N12
    P1   = {}
    P2   = {}
    P12  = {}
    N1s  = 0
    N2s  = 0
    N12s = 0

    print "Computing P1  ..."
    for k1 in N1.keys():
        N1s     += N1[k1]
    for k1 in N1.keys():
        P1[k1]   = float(N1[k1])/float(N1s)
    print "P1.values()=", P1.values()

    print "Computing P2  ..."
    for k2 in N2.keys():
        N2s     += N2[k2]
    for k2 in N2.keys():
        P2[k2]   = float(N2[k2])/float(N2s)
    print "P2.values()=", P2.values()

    print "Computing P12 ..."
    for k12 in N12.keys():
        N12s    += N12[k12]
    for k12 in N12.keys():
        P12[k12] = float(N12[k12])/float(N12s)
    print "P12.values()=", P12.values()
#   print "N1s=", N1s, " N2s=", N2s, " N12s=", N12s
    print "len(N1.keys())=", len(N1.keys())
    print "len(N2.keys())=", len(N2.keys())
    print "len(N12.keys())=", len(N12.keys())

    # Compute VI
    H1 = 0
    H2 = 0
    I  = 0
    for k in P1.keys():
        if P1[k] > 0:
            H1 += -P1[k]*math.log(P1[k])
        else:
            print "Warning: P1[", k, "]=", P1[k]

    print "H1=", H1

    for k in P2.keys():
        if P2[k] > 0:
            H2 += -P2[k]*math.log(P2[k])
        else:
            print "Warning: P2[", k, "]=", P2[k]
    print "H2=", H2

    nsteps = 0
    max_steps = len(N12.keys())
    VI_U = 0   # under-segmentation
    VI_O = 0   # over-segmentation
    for k1 in N1.keys():
        for k2 in N2.keys():
            k12 = str(k1) + "_" + str(k2)
            try:
                if P1[k1] > 0 and P2[k2] > 0 and P12[k12] > 0:
                    I += P12[k12]*math.log(P12[k12]/P1[k1]/P2[k2])
                    nsteps += 1
                    VI_U += -P12[k12]*math.log(P12[k12]/P1[k1])
                    VI_O += -P12[k12]*math.log(P12[k12]/P2[k2])
            except:
                print "NOTE: failed to compute I term for P1=", P1[k1], " P2=", P2[k2], " P12=", P12[k12]
                continue
           
    print "I=", I

    VI = (H1 + H2 - 2*I)

    return (VI, VI_U, VI_O)

# -------------------------------------------------------------------------

def counts2VI3D(N1, N2, N12):
    # Compute probabilities
    P1   = {}
    P2   = {}
    P12  = {}
    N1s  = 0
    N2s  = 0
    N12s = 0

    print "Computing P1  ..."
    for k1 in N1.keys():
        N1s     += N1[k1]
    for k1 in N1.keys():
        P1[k1]   = float(N1[k1])/float(N1s)
#   print "k1=", k1, " N1=", N1[k1], " P1[k1]=", P1[k1]

    print "Computing P2  ..."
    for k2 in N2.keys():
        N2s     += N2[k2]
    for k2 in N2.keys():
        P2[k2]   = float(N2[k2])/float(N2s)
#   print "k2=", k2, " N2=", N2[k2], " P2[k2]=", P1[k2]

    print "Computing P12  ..."
    for k12 in N12.keys():
        N12s    += N12[k12]
    for k12 in N12.keys():
        P12[k12] = float(N12[k12])/float(N12s)
#   print "k12=", k12, " N12=", N12[k12], " P12[k12]=", P12[k12]
    print "N1s=", N1s, " N2s=", N2s, " N12s=", N12s
    print "len(N1.keys())=", len(N1.keys())
    print "len(N2.keys())=", len(N2.keys())
    print "len(N12.keys())=", len(N12.keys())

    # Compute VI3D
    H1 = 0
    H2 = 0
    I  = 0


    print "\nComputing H1  ..."
    for k in N1.keys():
        if P1[k] > 0:
            H1 += -P1[k]*math.log(P1[k])
        else:
            print "Warning: P1[", k, "]=", P1[k]
    print "H1=", H1

    print "\nComputing H2  ..."
    for k in N2.keys():
        if P2[k] > 0:
            H2 += -P2[k]*math.log(P2[k])
        else:
            print "Warning: P2[", k, "]=", P2[k]
    print "H2=", H2
    VI_U = 0   # under-segmentation
    VI_O = 0   # over-segmentation

    print "\nComputing I  ..."
    for k1 in N1.keys():
        for k2 in N2.keys():
            k12 = str(k1) + "_" + str(k2)
            try:
                if P1[k1] > 0 and P2[k2] > 0 and P12[k12] > 0:
                    I += P12[k12]*math.log(P12[k12]/P1[k1]/P2[k2])
                    VI_U += -P12[k12]*math.log(P12[k12]/P1[k1])
                    VI_O += -P12[k12]*math.log(P12[k12]/P2[k2])
            except:
                continue
    print "I=", I

    VI3D = H1 + H2 - 2*I

    return (VI3D, VI_U, VI_O)


