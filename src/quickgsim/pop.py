
from . import mate

class Pop:
    """ 
    
    """
    def __init__(self,traits,backend=None):
        if not isinstance(traits,list):
            traits = [traits]
        self.traits = {tr.name:tr for tr in traits}
        self.backend = backend
        self.animals = {}
    
    def register_animals(self,animals):
        for a in animals:
            for trname,tr in self.traits.items():
                a.bvs[trname] = tr.draw_bv()
            self.animals[a.tag] = a
    
    def produce_newAnimal(self,sire,dam,kid_tag,genome,kid_sex=None,mod_len_genomes={'paternal':0.3,'maternal':0},fixed_effects={},progeny_inbreeding=0):        
        ## or sire adn dam, get a "phenotype" for mod_len_genomes
        
        
        new_animal = mate(sire,dam,kid_tag,genome,kid_sex,mod_len_genomes=mod_len_genomes)
        for trname,tr in self.traits.items():
            ## depending on the type of the trait, calcualte the breeding value (bv) of the trait
            if tr.bv_type == "pedigreeBased":
                bv = tr.get_progeny_bv(sire.bvs[trname],dam.bvs[trname],fsire=sire.f,fdam=dam.f)
            else:
                bv = tr.draw_bv()
            ## update the BVs
            new_animal.bvs[trname]  = bv
            new_animal.phens[trname] = tr.get_phen(bv,fixed_effects.get(trname,0))
        new_animal.f = progeny_inbreeding
        
        if self.backend:
            self.backend.store_animal(tag,sire,dam,sex,new_animal)
        
        self.animals[kid_tag] = new_animal
        return new_animal
            
            
            
    
