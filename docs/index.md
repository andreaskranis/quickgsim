Hello Sim!

Welcome to the **QuickgSim** Documentation

Quikgsim was designed to be easily extensible to enable users to add functionality to fit specific requirements. The class diagram below depicts the overall structure between the classes that make up QuickgSim.

``` mermaid
classDiagram

Genotype *-- Animal: has a genotype
Phenotype *-- Animal: optional phenotype

Variant *-- Chrom: has many in self.variants

Chrom *-- Genome: has many in self.chroms
Chrom: +name
Chrom: +length
Chrom: +morgans
Chrom: +nvars
Chrom: +variants
Chrom: +pos_idx
Chrom: +xover_pos



Genome --> Genotype : blueprint for genotype
Genome: +chroms
Genome: +rs

Genome <|-- CustomGenome: modify recomb and/or mut
CustomGenome: +recomb_events()
CustomGenome: +place_recomb()
CustomGenome: +mutate_gamete()



SQLITE_bend --> Animal : persistant storage 



```

### Extend functionality
As seen in the class digram, it is possible to extend and modify the functionality of the *Genome* class. In the `examples/` folder the notebook *customGenome.ipynb* shows two examples on how to capitalise on inheritance to tweak the model used to simulate recombination. The first example shows how to approximate the differences in recombination rate in telomeres and centromeres. The second example uses the Bionomial distribution to determine the number of cross-overs instead of the Poisson in the base Class *Genome*.



### safe to ignore for nopw
UserDict <|-- Genome
UserDict <|-- Genotype

DataClass <|-- Animal
DataClass <|-- Chrom

