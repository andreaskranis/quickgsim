# -*- coding: utf-8 -*-
"""Implementing the Animal class

"""
from dataclasses import dataclass,field
from typing import List, Dict

@dataclass
class Animal:
    """Dataclass implementation to provide alternative access to Genotype() 

    Attributes:
        tag (str): An identifier for the animal. Useful for retrieving in 
                   backend
        sex (int): The sex of the animal [1:male, 2:female]
        genotype(dict): A dictionary compatible to a Genotype (ie a UserDict)
        pcr(list): paternal xovers [optional]
        mcr(list): maternal xovers [optional]
        phens(dict): phenotypes for traits [optional]
        bvs(dict): breeding values for traits [optional] 
        f(float): inbreeding coefficient [optional]
    """
    tag: str
    sex : int
    genotype : Dict
    pcr: List[float] = field(default_factory=list)           #optional to store paternal xovers
    mcr: List[float] = field(default_factory=list)           #optional to store maternal xovers
    phens: Dict[str,float] = field(default_factory=dict)     #optional to store phenotypes for traits
    bvs: Dict[str,float] = field(default_factory=dict)       #optional to store breeding values (BVs) for traits
    f: float=0.0                                             #optional inbreeding                             
    
    def __repr__(self):
        return f"Animal(tag:{self.tag})"
    
    def get_haplo_genotype(self,strand):
        return self.genotype.as_haplo(strand)

    def get_array_genotype(self):
        return self.genotype.as_snps()

    def get_phenotype(self,trait,genome,fixed_effect=0):
        if not self.phens.get(trait.name,None):
            #self.phens[trait.name] = trait.mu + self.get_bv(trait,genome) + fixed_effect + trait.draw_res()        
            self.phens[trait.name] = trait.get_phen(self.get_bv(trait,genome),dbv=0,fixed=fixed_effect)
            
        return self.phens[trait.name]
    
    def get_bv(self,trait,genome):
        if not self.bvs.get(trait.name,None):
            bv = 0
            snp_gen,n = self.genotype.as_snps()
            for chrom in snp_gen:
                if trait.qtls.get(chrom,None):
                    bv += trait.bv_per_chrom(snp_gen,genome,chrom,total=True)
            self.bvs[trait.name] = bv
        return self.bvs[trait.name]
    
    def add_pheno(self,name,val):
        self.phens[name] = val
    
    def add_bv(self,name,val):
        self.bvs[name] = val
    
    def get_Ncrossovers(self,parent='paternal'):
        if parent == 'paternal':
            return sum(len(self.pcr[c][1]) for c in self.pcr)
        elif parent == 'maternal':
            return sum(len(self.mcr[c][1]) for c in self.mcr)
        else:
            return 0

        
        
