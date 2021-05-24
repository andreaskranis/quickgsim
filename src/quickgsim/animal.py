from . import RAN_GEN, np

class Animal:
    
    def __init__(self,tag,sex,genotype,pcr=[],mcr=[],seed=None):
        self.tag = tag
        self.sex = sex
        self.genotype = genotype
        self.paternal_crossovers = pcr
        self.maternal_crossovers = mcr

    
    def get_haplo_genotype(self,strand):
        g1,N = self.genotype.as_haplo(strand)
        return self.genotype.vectorise_snp_gen(g1,N)
        
    def get_array_genotype(self):
        g1,N = self.genotype.as_snps()
        return  self.genotype.vectorise_snp_gen(g1,N)
        

