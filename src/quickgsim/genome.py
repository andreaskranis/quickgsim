"""


"""


from dataclasses import dataclass,field
from typing import List, Dict
from . import RAN_GEN, Genotype, np


@dataclass
class Variant:
    chrom : str
    pos : int
    cm_pos : float = 0
    add_effects : Dict[str,float] = field(default_factory=dict)
    
    def get_name(self):
        return f"{self.chrom}_{self.pos}"



@dataclass
class Chrom:
    name : str
    length : int
    morgans: float
    nvars : int = 0
    variants: List[float] = field(default_factory=list)
    pos_idx : Dict[int,int] = field(default_factory=dict)
    xover_p: List[float] = field(default_factory=list) or None
    
    
    def add_variant(self,pos=None,cm_pos=None):
        self.variants.append(Variant(self.name,pos,cm_pos))
        self.pos_idx[pos] = self.nvars
        self.nvars += 1
        if cm_pos>=0:
            self.xover_p.append(cm_pos) 
        else:
            self.xover_p.append(pos) 
        
    def finalise_chrom_configuration(self):
        if self.xover_p:
            self.xover_p = np.array(self.xover_p)
            self.nvars = len(self.variants)

            if not self.morgans:
                self.morgans = max(self.xover_p)/100
        else:
            self.xover_p = [v.pos for v in self.variants]  ## if no cm is available, that physical position is a 1:1 proxy for cm
        ## normalise the map distances to probabilities [NOTE probably exclude the first and last SNP]
        self.cm_to_p()
           
    def cm_to_p(self,skip_first_last=True):
        ##check https://numpy.org/doc/stable/reference/generated/numpy.diff.html
        self.xover_p = np.diff(self.xover_p,prepend=0)/(self.morgans*100)
        if skip_first_last:
            self.xover_p = self.xover_p[1:-1]/self.xover_p[1:-1].sum()
        else:
            self.xover_p /= self.xover_p.sum()
            
##NOTE: 
#    def p_to_cm(self):
#        self.cm_pos = np.cumsum(self.cm_pos*self.morgans*100)
    
    def get_variant_by_pos(self,pos,mv=None) -> Variant:
        """Returns a variant in a specific <pos>ition"""
        return self.pos_idx.get(pos,None)
    
    def get_variant_by_ordinal(self,nth_pos,mv=None) -> Variant:
        """Returns the n-th variant"""
        try:
            return self.variants[nth_pos]
        except IndexError:
            return mv
        
        


class Genome:
    SWITCH = {0:1,1:0}
    
    def __init__(self,rs=None):
        self.chroms = {}
        self.rs = rs if rs else RAN_GEN
        
    def add_chrom(self,chrom_name,length,morgans):
        if chrom_name not in self.chroms:
            self.chroms[chrom_name] = Chrom(chrom_name,length,morgans)               
            
    def add_variant(self,chrom_name,pos,cm_pos=None):
        try:
            self.chroms[chrom_name].add_variant(pos,cm_pos)
        except KeyError:
            print(f"**WARN: chromosome {chrom_name} has not been registered to current genome object. Use add_chrom() before adding variant {chrom_name}_{pos}")
    
    def __repr__(self):
        info  = []
        for c in self.chroms:
                info.append(f"- Chromosome {c} contains {self.chroms[c].nvars} variants and is {self.chroms[c].morgans} MAP-UNITS long")
        return "\n".join(info)



    ##<Recombination> 
    def recomb_events(self,chrom_name,obligatory=1,max_events_asPropChromLen=None):
        """Determines the number of cross-overs in a chromosome drawing from a poisson distribution
           Allows to specify the maximun number of events as a proportion of the chromosome length (could use 3 or 4 to limit outliers)
        """
        if chrom_name in self.chroms:
            r = self.rs.poisson(self.chroms[chrom_name].morgans)
            if max_events_asPropChromLen:
                max_recomb = self.chroms[chrom_name].morgans * max_events_asPropChromLen
                if r > max_recomb:
                    r = max_recomb
            return r if r >= 1 else obligatory
        return None
    
    def place_recomb(self,chrom,n_recomb):
        """Determines the location of the crossovers drawing from a uniform distribution and using as weights the chrom.xover_p"""
        n_markers = chrom.nvars
        #Ensure that placing xoversa makes sense: I cannot do with less than 3 markers or if xovers are more than non-edge SNPs
        if (n_markers > 2) and (n_recomb < n_markers - 2):        
            crossovers = self.rs.choice(np.arange(1,n_markers-1),n_recomb,p=chrom.xover_p,replace=False)  ##avoid first and last markers to have clear recomb
            crossovers.sort()
            return crossovers.astype(int)
        return []
        
    ##<Mutation>
    def mutate_gamete(self,gamete):
        """would switch alleles ONLY, not introducing new mutations"""
        return gamete
    ##</Mutation>
    
    def get_gamete(self,genotype):
        gamete = {c:[] for c in genotype.iterate_chroms()}
        crossovers = {c:[[],[]] for c in genotype.iterate_chroms()}
        
        for c in genotype.iterate_chroms():
            n_markers = self.chroms[c].nvars
            gamete[c] = np.empty(n_markers,dtype=int)
            
            ## determine recombination paramters for current chromosome c
            n_recomb = self.recomb_events(c,obligatory=1)
            breaks = self.place_recomb(self.chroms[c],n_recomb)
            st_strand = self.rs.integers(0,2)
            
            ## Collate the gamete and record break points
            gamete[c][:breaks[0]] = genotype[c][st_strand][:breaks[0]]
            last_break = breaks[0]
            for i,b in enumerate(breaks[1:]):
                new_strand = self.SWITCH[(st_strand+i)%2]
                gamete[c][last_break:b] = genotype[c][new_strand][last_break:b]
                last_break = b
            
            ## Glue the last haplotype block after the last cross-over
            if len(breaks)%2 == 0:
                last_strand = st_strand
            else:
                last_strand = self.SWITCH[st_strand]
            gamete[c][last_break:] = genotype[c][last_strand][last_break:]
            crossovers[c] = [[st_strand,last_strand],breaks]
        
        return gamete,crossovers

    def get_zygote(self,pat_gam,mat_gam):
        """BEWARE: No checks for dimension compatability """
        genotype = Genotype(self.chroms)
        
        for c in self.chroms.keys():
            genotype[c][0] = pat_gam[c]
            genotype[c][1] = mat_gam[c]
        return genotype