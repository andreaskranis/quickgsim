"""
"""
from . import Animal,np,tqdm


def create_founders(genders,gens,genome):
    """Can use either 'real' genotypes or artificial"""
    founders = {}
    a1 = set(genders.keys())
    a2 = set(gens.keys())
    print(f"Genotypes and gender info for {len(a2)} and {len(a1)} founders was received")
    
    common = a1.intersection(a2)
    print(f"** Note: {len(common)} instances of Animal() will be created")
    for a in common:
        founders[a] = Animal(a,genders[a],gens[a])
    return founders


def mate(sire,dam,kid_tag,genome,kid_sex=None,mod_len_genomes={'paternal':0,'maternal':0}):
    pg,pcr = genome.get_gamete(sire.genotype,mod_len_genome=mod_len_genomes.get('paternal',0))
    mg,mcr = genome.get_gamete(dam.genotype,mod_len_genome=mod_len_genomes.get('maternal',0))
    
    if not kid_sex:
        kid_sex = genome.rs.integers(1,3,dtype=int)
    zygote = genome.get_zygote(pg,mg)
        
    return Animal(kid_tag,kid_sex,zygote,pcr,mcr)  ##new kid


def drop_pedigree(pop,genome,ped,store_backend=None):
    for kid,sire,dam,sex in tqdm.tqdm(ped):
        progeny = mate(pop[sire],pop[dam],kid,genome,sex)
        pop[kid] = progeny
        if store_backend:
            store_backend.store_animal(kid,sire,dam,sex,progeny)
            
    if store_backend:
        store_backend.finalise()


