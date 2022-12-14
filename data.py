import pyGPs
import sobol_seq
import numpy as np
import os
import numpy as np


benchmark = "benchmark"
epsilon  = 1e-6

def kriging_instance(seed=None, sobol=True, N=100, d=2, m=10):
    """
    generate the design matrix A of a optimal design instance
    obtained by truncating 
    * m eigenvalues 
    * of an RBF kernel in d dimensions
    * evaluated at N >> m candidate points
    * the candidate points are either from a regular grid or a sobol sequence
    
    returns:
    * the coordinates `xxx` of the candidate support points (an N x d)-matrix 
    * A matrix `A` such that the elementary design matrix of the ith candidate point xxx[i] is equal to 
      the tensor product of the ith column of `A` with itself:
          M[i] = np.tensordot(A.T[i],A.T[i],0)
      Note that `M[i]` is of size (m x m) and `A` has shape m x N
    """


    np.random.seed(seed)

    if sobol:
        xxx = sobol_seq.i4_sobol_generate(d,N)
    else:
        N1ovd = N**(1./d)
        N1ovd_rounded = int(np.round(N1ovd))
        N = N1ovd_rounded**d
        if N1ovd_rounded != N1ovd:
            print(f'N needs to be a dth power when sobol=False, so it will be rounded to {N}') 
        xx = np.linspace(0,1,N1ovd_rounded)
        msh = np.meshgrid(*([xx] * d))
        xxx = msh[0].ravel()
        for k in range(1, d):
            xxx = np.c_[xxx, msh[k].ravel()]

    assert N>m, 'not enough support points'

    #log-uniform between ~ 0.3 and 0.5
    log_ells = 0.5*np.random.rand(d)-1.2
    sigma = 1.
    log_sigma = np.log(sigma)

    Ker = pyGPs.cov.RBFard(log_ell_list=list(log_ells), log_sigma=log_sigma)
    K = Ker.getCovMatrix(xxx,xxx, 'train')

    #take a unit mu
    n = K.shape[0]
    mu = np.ones(n)

    #tmp,to test
    #mu = np.ones(n) + np.random.rand(n)

    D = np.diag(mu)
    Lb, U = np.linalg.eigh(D.dot(K).dot(D))
    V = np.diag(mu**(-0.5)).dot(U)

    Phi = V[:, -m:]
    Lm = Lb[-m:]
    print('{0:.3f}% of spectrum retained'.format(100*sum(Lm)/sum(Lb)))

    sg2 = mu-(np.linalg.norm((Phi*Lm**0.5), axis=1)**2)
    A = (Phi * Lm ** 0.5).T / sg2 ** 0.5
    
    return xxx, A


def compute_block_design_data(t=6,sz=2,levi=False):
    """
    computes the design matrix for a block design problem (with blocks of size 2) on `t` treatments
    returns:
        * list_of_pairs: list of pairs of treatments (these are the 'candidate points' of the design problem)
        * A matrix `A` such that the elementary design matrix of the ith pair of treatments list_of_pairs[i] is equal to 
          the tensor product of the ith column of `A` with itself:
              M[i] = np.tensordot(A.T[i],A.T[i],0)
          Note that `M[i]` is of size (m x m) and `A` has shape m x N, where N=t*(t-1)/2 is the number of pairs of treatments
          and m=t-1 is the number of unknown parameters (one less than number of treatments, since only differences can be estimated,
          so we have to fix the value of some arbitrary treatment to some constant)
    """
    list_of_pairs=[]
    
    AA2=[]
    ind=0
    I,J,O=[],[],[]
    for i in range(t):
        for j in range(t):
            if j>i:
                list_of_pairs.append((i,j))
                #pair2id[i,j]=ind
                gij= np.zeros(t)
                gij[i]=1
                gij[j]=-1
                AA2.append(gij[:-1])
    A = np.array(AA2).T
    return list_of_pairs, A
    

def norm_instance(N, m):
    A =  np.random.normal(0, 1 / N**(0.5), (N, m))
    return A

# generate normal distribution instances
def generate_norm():
    # clean 
    for filename in os.listdir(benchmark):
        if filename[0] == "n":
            os.remove(benchmark + "/" + filename)
    Ns = [50, 60, 70]
    ms = [0.4, 0.6]

    m = {}
    candidates = {}
    As = {}
    addcaps = [0, 1, 2, 3, 4]
    for N in Ns:
        for m in ms:
            m = int( m * N)
            for addcap in addcaps:
                #print(m)
                A = norm_instance( N,  m)
                cap = m + addcap
                file_name = "normal" + "_" + str(N) + "_" + str(m) + "_" + str(cap) +  ".design"
                print(file_name)
                descriptor = str(N) + " "  + str(m)  +  " " + str(cap) + " " + str(epsilon) 
                f = open(benchmark + "/" + file_name, "x")
                f.write(descriptor +  '\n')
                for i in range(N):
                    for j in range(m):
                        f.write(str(A[i][j]) + " ")
                    f.write("\n")
                f.close()  

# generate kriging instances

def generate_kriging():
    # clean 
    for filename in os.listdir(benchmark):
        if filename[0] == "k":
            os.remove(benchmark + "/" + filename)
    Ns = [50, 70]
    ds = [2, 4]
    ms = [0.4, 0.6]

    m = {}
    candidates = {}
    As = {}
    addcaps = [0, 1, 2, 3, 4]
    for N in Ns:
        for d in ds:
            for m in ms:
                m = int( m * N)
                for addcap in addcaps:
                    #print(m)
                    xxx, A = kriging_instance(None, True, N, d, m)
                    candidates[(N,d,m)] = xxx
                    As[(N,d,m)] = A
                    cap = m + addcap
                    file_name = "krigging" + "_" + str(N) + "_" +  str(d)  +  "_"  + str(m) + "_" + str(cap)  +  ".design"
                    descriptor = str(N) + " "  + str(m) + " " +  str(cap) + " " + str(epsilon) 
                    f = open(benchmark + "/" + file_name, "x")
                    f.write(descriptor +  '\n')

                    for i in range(N):
                        for j in range(m):
                            f.write(str(A[j][i]) + " ")
                        f.write("\n")
                    f.close()   

            
            
    

# generate block of size 2 instances

def generate_block():
    for filename in os.listdir(benchmark):
        if filename[0] == "b":
            os.remove(benchmark + "/" + filename)
    candidates = {}
    As = {}
    ts = [10, 11, 12]
    addcaps = [0, 1, 2, 3, 4]
    for t in ts:
        for addcap in addcaps:
            list_of_pairs, A = compute_block_design_data(t)
            candidates[t] = list_of_pairs
            As[t] = A
            N= int(t*(t-1)/2)
            m = t-1
            cap = m + addcap
            file_name = "block2" + "_" + str(N) + "_" +  str(t)  +  "_"  + str(m) + "_" + str(cap) + ".design"
            descriptor = str(N) + " "  + str(m) + " " +  str(cap) + " " + str(epsilon) 
            f = open(benchmark + "/" + file_name, "x")
            f.write(descriptor +  '\n')
            for i in range(N):
                for j in range(m):
                    f.write(str(A[j][i]) + " ")
                f.write("\n")
            f.close()   


#generate_kriging()

generate_block()

generate_norm()
