# GeneticDrift
C program that simulates genetic drift and determines fixation time of an allele. 
Illustrates how impact of genetic drift depends on population size.
The complete documentation is in the actual program file, but the key aspects are repeated here:

Purpose:
Generate data on how the impact of genetic drift varies 
by population size. Specifically, the program determines the 
average time to fixation of an allele for a user-specified 
set of population sizes. Each user-specified population size is
simulated a user-specified number of times, and the average time 
to fixation reported. The results were used in an Evolutionary 
Genetics lecture at Santa Clara University.

Biological Specifications:
Two alleles (I have a separate version that does three alleles)
User-specified genotype frequencies
    Thus user-specified allele frequencies
    Default is 30% homozygous for "allele 1", 40% heterozygous, 30% homozygous for "allele 2"
         Thus default allele frequencies are 50% allele 1, 50% allele 2
    Default is to do 10 simulations each of population sizes 10, 100, and 1000
    Random mating
    Fixed population size
    I.e., everyone survives, reproduces, and has two offspring
    
Limitations and Possibilities:
    Currently, must change parameters in the code itself (though they are at least all globals at start of program).
    Population sizes need to be even numbers, but this is not checked.
    Does not save results, just prints them to screen.
    No GUI.
    No error handling to speak of.
    No unit tests.
    
