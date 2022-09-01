"""


"""

from . import RAN_GEN,np

class Trait:

    def __init__(self,name="trait",mu=0,rs=RAN_GEN):
        self.name = name
        self.mu = mu
        self.rs=rs
        self.stdA,self.stdE,self.stdD,self.stdI = map(np.sqrt,[self.Va,self.Ve,self.Vd,self.Vi])

    
    def _calc_h2(self):
        return self.Va/(self.Va+self.Vd+self.Ve)

    def _draw_sample(self,mu,stdev):
        return self.rs.normal(mu,stdev)
    
    def draw_bv(self):
        return self._draw_sample(0,self.stdA)

    def draw_res(self):
        return self._draw_sample(0,self.stdE)

    def get_phen(self,bv,dbv=0,fixed=0):
        return self.mu + fixed + bv + dbv + self.draw_res()


class pedTrait(Trait):
    
    def __init__(self,name,Va,Ve,Vd=0,Vi=0,mu=0,rs=RAN_GEN):
        super().__init__(name,mu,RAN_GEN)
        self.Va = Va
        self.Ve = Ve
        self.Vd = Vd
        self.Vi = Vi
        #self.stdA,self.stdE,self.stdD,self.stdI = np.sqrt(Va),np.sqrt(Ve),np.sqrt(Vd),np.sqrt(Vi)
        self.rs = rs
        self.h2 = self._calc_h2()
        self.bv_type = "pedigreeBased"
        super().__init__(name,mu,rs)
    
    def get_progeny_bv(self,bv_sire,bv_dam,fsire=0,fdam=0):
        """Progeny BV is parental average +/ the mendelian sampling (MS)
           I assume that var(MS)=0.5*Va*(1- (fsire+fdam)/2)
        """
        return 0.5*(bv_sire+bv_dam) + self._draw_sample(mu=0,stdev=0.5*self.stdA*(1-(1-(fsire+fdam)/2)))
    
    
class genomicTrait(Trait):
    
    def __init__(self,name,Va,Ve,genome,mu=0,rs=RAN_GEN):
        self.bv_type = "genomeBased"
        self.qtls ={}
        self.nqtls = 0
        self.Va = Va
        self.Ve = Ve
        self.Vd = 0
        self.Vi = 0
        self.h2 = self._calc_h2()
        super().__init__(name,mu,rs)
        
    def generateQTL_effects(self,genome):
        for chrom in self.qtls:
            nq = len(self.qtls[chrom])
            additive = self.rs.gamma(shape=0.2, scale=5, size=nq) * self.rs.permutation(np.repeat([1,-1],nq))[:nq]
            for i,snp_idx in enumerate(self.qtls[chrom]):
                genome.chroms[chrom].variants[snp_idx].add_effects.update({self.name:additive[i]})            
    
    def assign_snps_QTLs(self,genome,nqtls):
        """Decide which Variants are QTLs for the trait """
        self.nqtls = nqtls
        n_qtls_chrom = {c:int(round(nqtls*genome.chroms[c].length/genome.total_len('bp'),0)) for c in genome.chroms}
        for c in n_qtls_chrom:
            if n_qtls_chrom[c] < genome.chroms[c].nvars:
                self.qtls[c] = sorted(self.rs.choice(range(genome.chroms[c].nvars),n_qtls_chrom[c],replace=False))
            else:
                self.qtls[c] = range(genome.chroms[c].nvars)    
        
    def get_progeny_bv(self,genome,genotype):
        """ [NOT implemented]
            It will use info on genome to get QTLs and calculate BV
        """
        bv = 0
        return bv  ##sum of QTL effects
    
    def get_Va_fromBVs(self,bvs):
        return bvs.mean(),bvs.var()
    
    
    def adjust_effects_va(self,founders,genome):
        """Adjusts the effects to have bv.var=Va. 
           Founders is a list of Animal() objects
        """
        bv = np.zeros(len(founders))
        gens = []
        for i,an in enumerate(founders):
            gens.append(an.genotype.as_snps()[0])
            for chrom in self.qtls:
                bv[i] += sum(self.bv_per_chrom(gens[-1],genome,chrom))

        ## now adjust the effects to ensure BVs of founders have mean=0 and var=self.Va
        m,v = self.get_Va_fromBVs(bv)
        bias = -1*m/self.nqtls   ## a value to adjust each effect to ensure that bv.mean=0           
        coef = np.sqrt(self.Va/v) ## the scaling coef to adjust effects in order bv.var=self.Va
        for chrom in self.qtls:
            for s in self.qtls[chrom]:
                genome.chroms[chrom].variants[s].add_effects[self.name] += bias
                genome.chroms[chrom].variants[s].add_effects[self.name] *= coef            
        
        ## ensure that the Va from the BVs is now equal to self.Va 
        bv1 = np.zeros(len(founders))
        for i in range(len(founders)):
            for chrom in self.qtls:
                bv1[i] += self.bv_per_chrom(gens[i],genome,chrom,total=True)

        m1,v1 = self.get_Va_fromBVs(bv1)     
        print(f"TRAIT {self.name} -> var_before:{v:.3f},var_after:{v1:.3f},target:{self.Va} [means:{m:.2f} -> {m1:.2f}]")

    
    def bv_per_chrom(self,gen_dict,genome,c,total=False):
        """Return the BV for a chromosome"""
        if total:
            return sum(np.array([v.get_effect(self.name) for v in genome.chroms[c].get_variants_attrs(self.qtls[c])]) * gen_dict[c][self.qtls[c]])
        else:
            return np.array([v.get_effect(self.name) for v in genome.chroms[c].get_variants_attrs(self.qtls[c])]) * gen_dict[c][self.qtls[c]]