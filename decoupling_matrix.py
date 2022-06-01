## Copyright 2022 JB Lagrange, STFC RAL

import numpy

def swap_eigen_order(eigenvalue, eigenvector, i,j):
    conj_change = [0,0,0,0]
    conj_mu = eigenvalue[i]
    eigenvalue[i] = eigenvalue[j];
    eigenvalue[j] = conj_mu
    for n in range(4): conj_change[n] = eigenvector[n,i]
    #print "temp"
    #print conj_change
    for n in range(4): eigenvector[n,i] = eigenvector[n,j]
    for n in range(4): eigenvector[n,j] = conj_change[n]

def decoupling_matrix(matrix):
    #print matrix
    eigenvalues, eigenvectors = numpy.linalg.eig(matrix)
    #print "original eigenvalues:"
    #for i in range(4): print eigenvalues[i]
    #print "original eigenvectors:"
    #print eigenvectors
    
    #pair the conjugated eigenvectors and eigenvalues
    k=0
    for i in range(1,4):
        if abs(eigenvalues[i]-numpy.conj(eigenvalues[0]))<1.e-9:
            k=i
            break
    if k==0:
        raise ValueError("problem in eigenvalues, not conjugates?")
    swap_eigen_order(eigenvalues, eigenvectors, 1, k)
    
    #print "after pairing eigenvalues:"
    #for i in range(4): print eigenvalues[i]
    #print "eigenvectors:"
    #print eigenvectors
    
    for j in range(0,4,2):
        # Unity test of absolute value
        if abs(abs(eigenvalues[j])-1)>1.e-3: raise ValueError("eigenvalue is not unity value:")
        #put >0 imaginary part first
        if abs(numpy.imag(eigenvalues[j]))>1.e-12 and numpy.imag(eigenvalues[j])<0:
            swap_eigen_order(eigenvalues, eigenvectors, j, j+1)
    #print "after unity test and imag>0 eigenvalues:"
    #for i in range(4): print eigenvalues[i]
    #print "eigenvectors:"
    #print eigenvectors
    for j in range(0,4,2):
        norm = 0
        for n in range(0,4,2): 
            #print numpy.real(eigenvectors[n,j]),"*",numpy.imag(eigenvectors[n+1,j]),"-",numpy.imag(eigenvectors[n,j]),"*",numpy.real(eigenvectors[n+1,j])
            norm += numpy.real(eigenvectors[n,j])*numpy.imag(eigenvectors[n+1,j]) - numpy.imag(eigenvectors[n,j])*numpy.real(eigenvectors[n+1,j])
        #print "ori norm:", norm
        #if norm<0, swap the order of the concerned eigenvector pair
        if norm<0:
            print("norm<0")
            norm = -norm
            swap_eigen_order(eigenvalues, eigenvectors, j, j+1)
        norm = numpy.sqrt(norm)
        #print "final norm", norm
        if abs(norm)<1.e-9:
            raise ValueError("normalisation value<1.e-9")
        
        for n in range(4):
            eigenvectors[n,j] /= norm
            eigenvectors[n,j+1] /= norm
        
    #print "after 1st normalisation eigenvalues:"
    #for i in range(4): print eigenvalues[i]
    #print "eigenvectors:"
    #print eigenvectors
    par_t_evector_inv = [[0+0j for i in range(4)] for j in range(4)]
    par_t_evector_inv = numpy.array(par_t_evector_inv)
    #print par_t_evector_inv
    for j in range(0,4,2):
        #print eigenvectors[j+1,j], eigenvectors[j,j]
        ratio = eigenvectors[j+1,j]/eigenvectors[j,j]
        beta = 1./numpy.imag(ratio)
        print("j=", j)
        #print "ratio:", ratio
        #print "beta:", beta
        if beta<0:
            print("j=",j, " beta<0, swap of eigenvectors and eigenvalues")
            #print "before swapping:"
            #print "eigenvectors"
            #print eigenvectors
            swap_eigen_order(eigenvalues, eigenvectors, j, j+1)
            #print "eigenvectors"
            #print eigenvectors
            ratio = eigenvectors[j+1,j]/eigenvectors[j,j]
            beta = 1./numpy.imag(ratio)
            #print "ratio: ", ratio, "beta", beta
            if beta<0: raise ValueError("beta still <0, decoupling not done")
        print("beta:", beta)
        phase = numpy.angle(eigenvectors[j,j])
        exp_phi = numpy.exp(1j*phase)
        exp_mphi = numpy.exp(-1j*phase)
        alpha = -beta*numpy.real(ratio)
        tune = numpy.arctan(numpy.imag(eigenvalues[j])/numpy.real(eigenvalues[j]))/(2*numpy.pi)
        #print "phase: ", phase
        #print "exp_phi: ", exp_phi
        #print "exp_mphi: ", exp_mphi
        #print "alpha: ", alpha
        print("tune: ", tune, "1-tune", 1-tune)
        par_t_evector_inv[j  ,j  ] = (-alpha-1j)*exp_mphi/(numpy.sqrt(beta));
        par_t_evector_inv[j+1  ,j] = -(-alpha+1j)*exp_phi/(numpy.sqrt(beta));
        par_t_evector_inv[j,j+1  ] = -numpy.sqrt(beta)*exp_mphi;
        par_t_evector_inv[j+1,j+1] = numpy.sqrt(beta)*exp_phi;
    print(par_t_evector_inv)
    r = numpy.dot(eigenvectors,par_t_evector_inv)
    print("before final norm")
    print(r)
    #final normalisation
    r /= -2j
    #for i in range(4):
    #    for j in range(4):
    #        r[i,j] /= -2j
    print("R:")
    print(r)
    print("det:", numpy.linalg.det(r))
    r_inv = numpy.linalg.inv(r)
    t = numpy.dot(r_inv, numpy.dot(matrix, r))
    print("T:")
    print(t)
    print("det:", numpy.linalg.det(t))

    return r, r_inv, t

def swap_eigenvalues_order(eigenvalue, i,j):
    conj_change = [0,0,0,0]
    conj_mu = eigenvalue[i]
    eigenvalue[i] = eigenvalue[j];
    eigenvalue[j] = conj_mu

def get_tunes(matrix):
    qu = 0
    qv = 0
    #print matrix
    eigenvalues = numpy.linalg.eigvals(matrix)
    #print "original eigenvalues:"
    #for i in range(4): print eigenvalues[i]
    
    #pair the conjugated eigenvectors and eigenvalues
    k=0
    for i in range(1,4):
        if abs(eigenvalues[i]-numpy.conj(eigenvalues[0]))<1.e-9:
            k=i
            break
    if k==0:
        raise ValueError("problem in eigenvalues, not conjugates?")
    swap_eigenvalues_order(eigenvalues, 1, k)
    
    #print "after pairing eigenvalues:"
    #for i in range(4): print eigenvalues[i]
    
    for j in range(0,4,2):
        # Unity test of absolute value
        if abs(abs(eigenvalues[j])-1)>1.e-3: raise ValueError("eigenvalue is not unity value:")
        #put >0 imaginary part first, 0<tunes<0.5
        if abs(numpy.imag(eigenvalues[j]))>1.e-12 and numpy.imag(eigenvalues[j])<0:
            swap_eigenvalues_order(eigenvalues, j, j+1)
    #print "after unity test and imag>0 eigenvalues:"
    #for i in range(4): print eigenvalues[i]
    #print "eigenvectors:"
    qu = numpy.arctan(numpy.imag(eigenvalues[0])/numpy.real(eigenvalues[0]))/(2*numpy.pi)
    qv = numpy.arctan(numpy.imag(eigenvalues[2])/numpy.real(eigenvalues[2]))/(2*numpy.pi)
    return qu, qv




matrix = [
    [ 0.768455,  27.219055, 0.926560,  -0.001729],
    [ -0.046495, 0.768996,  0.000013,  -0.924621],
    [ 0.923850,  -0.001668, 0.768089,  27.145483],
    [ -0.000013, -0.926456, -0.046638, 0.767885 ],
]
#r = None
#r_inv = None
#t = None
#decoupling_matrix(matrix, r, r_inv, t)

qu, qv = get_tunes(matrix)
print(qu, qv)
