
"""Provides functions to generate SNP panels with fixed characteristics.

The modules also provides functions to validate the imported data to ensure
that are suitable for quicgsim.

Todo:
    * For module TODOs

"""

from .genome import Genome
from .genotype import Genotype


class GenomeSetupError(Exception):
    def __init__(self,message):
        super().__init__(message)
        



def _fill_equidistant_variants(instructions,genome):
    """Evenly-spaced variants across chromosomes 
    
    Args:
        instructions (dict): A dictionary with key the chrom name and value a list with two elements: (i) number of variants (ii) the length in morgans
        genome (Genome): The genome object
    """
    for chrom in instructions:
        nsnps,length,morgans = instructions[chrom]
        if nsnps < 3:
            raise GenomeSetupError(f"Cannot initiate chromosome {chrom} with less than 3 SNPs!")
        variant_len_space = length/(nsnps-1)
        variant_mor_space = morgans/(nsnps-1)
        print(f"will now add {nsnps} variants in chromosome {chrom} [len in bp:{length}/morgans:{morgans}] equally spaced every {variant_len_space} bp, ie {variant_mor_space} map units")
        genome.add_chrom(chrom,length,morgans)
        
        for i in range(nsnps):
            genome.add_variant(chrom,pos=i*variant_len_space,cm_pos=i*variant_mor_space*100)    #mult by 100 to convert to centi-morgans
        print(f" chrom {chrom} has {len(genome.chroms[chrom].variants)} variants")
    for chrom in genome.chroms:
            genome.chroms[chrom].finalise_chrom_configuration()

def equidistant_genome(instructions,genome=None,rs=None):
    """A genome with variants placed in equal distances 
    
    Args:
        instructions (dict): A dictionary with key the chrom name and value a list with two elements: (i) number of variants (ii) the length in morgans
        genome (Genome): An optional instance of a genome object. If none is given, it will use the genome.Genome() baseclass
        rs (np.random.default_rng, optional): A random number generator, conveniently from randomise()
    """
    if genome:
        g = genome
    else:
        g = Genome(rs=rs)
    _fill_equidistant_variants(instructions,g)
    return g
            
    
    
def generate_ancestors():
    pass



