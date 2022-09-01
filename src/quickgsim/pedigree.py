"""


"""
from . import np

from collections import OrderedDict
from calcnrm import *

try:
    from numba import njit
    from numba.typed import List
except ImportError:
    njit = lambda x:x
    List = list()


class Pedigree:
    
    def __init__(self):
        self.ped = OrderedDict()
        
    def read_ped_from_file(self,inFile,sep=" ",header=False):
        sort_set_toList = lambda x: sorted(list(x))
        sires,dams = set(),set()
        tmp = {}
        with open(inFile) as fin:
            if header:
                next(fin)
            for row in fin:
                a,s,d = map(int,row.strip().split(sep))
                tmp[a] = [s,d]
                sires.add(s)
                dams.add(d)
        basesires = sort_set_toList(sires.difference(tmp.keys()))
        basedams = sort_set_toList(dams.difference(tmp.keys()))
        for s in basesires:
            self.ped[s] = [0,0]
        for d in basedams:
            self.ped[d] = [0,0]
        self.ped.update(tmp)
        print(f"Pedigree inlcudes {len(self.ped)} animals. In total, there are {len(sires)} sires and {len(dams)} dams and includes {len(basesires)} and {len(basedams)} base-sires/dams")
        return sires,dams,basesires,basedams
        
    def core_ped(self,parents):
        """Returns a reduced pedigree with only the parents"""
        out,sires,dams = OrderedDict(),set(),set()
        print(f"Will reduce pedigree size of {len(self.ped)} by keeping only the {len(parents)} parents provided")
        for a,(s,d) in self.ped.items():
            if a in parents:
                out[a] = [s,d]
                sires.add(s)
                dams.add(d)
        print(f"redcuded pedigree has now {len(out)} animals coming from {len(sires)} sires and {len(dams)} dams")
        return out,sires,dams
    
    def get_Amat(self,ped,method="numba"):
        unknown_parent,sires,dams = 0,List(),List()
        recoded_tags = {a:i for i,a in enumerate(ped.keys())}
        n = len(ped)
        for a,(s,d) in ped.items():
            sires.append(recoded_tags.get(s,n))
            dams.append(recoded_tags.get(d,n))
        if method == "numba":
            A = calcAmat(sires,dams,n)
        elif method == "fortran":
            A = calcnrm(sires,dams,n)
        else:
            A = calcAmat_python(sires,dams,n)
        return A[:-1,:-1]
        
    def get_inbreeding(self,Amat):
        return Amat.diagonal()
        

        
@njit
def calcAmat(sires,dams,n):
    N = n+1
    A = np.zeros((N,N),dtype=np.float32)
    for i in range(n):
        A[i,i] = 1 + 0.5*A[sires[i],dams[i]]
        for j in range(i+1,n):
            if j > n:
                break
            A[i,j] = 0.5*(A[i,sires[j]] + A[i,dams[j]])
            A[j,i] = A[i,j]
    return A

def calcAmat_python(sires,dams,n):
    N = n+1
    A = np.zeros((N,N),dtype=np.float32)
    for i in range(n):
        A[i,i] = 1 + 0.5*A[sires[i],dams[i]]
        for j in range(i+1,n):
            if j > n:
                break
            A[i,j] = 0.5*(A[i,sires[j]] + A[i,dams[j]])
            A[j,i] = A[i,j]
    return A
