/************************************
GENETIC DRIFT SIMULATOR (version 3.1)
*************************************
By David DeGusta (david.degusta@gmail.com)
Written July 3, 2014 to July 4, 2014
    with minor upgrades from June 29, 2015 to July 1, 2015
In C, using Codelite IDE and GCC compiler
MIT License (see below)

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
    
    
>>PSUEDOCODE STRUCTURE OF MAIN FUNCTION<<

For loop to iterate across array of population sizes user wants to simulate (Population_Size)
    For loop where each loop represents one simulation of a given population size
        Pass the population size and an empty population array to the create_initial_population function
            That function fills population array with elements, each element represents an individual
        Set fixation time (fixtime) to zero generations, will count up until an allele reaches fixation
        While check_fixation function returns "not fixed" flag:
            Pass population to make_new_generation_function
                That function randomly "mates" the individuals to produce a new generation
            Increase the fixation time counter (fixtime) by one
        Pass fixation time counter, population size, and simulation number to write_to_result function
            That function stores that information in the raw_results array
Pass raw_results array to calculate_final_results 
    That function calculates the overall results (averages etc) from the raw fixation times
Print results to screen
end 


>>FUNCTIONS<<

create_initial_population 
    called from main
    takes size (integer) and pointer to population (integer) array
    populates population array with members
    n.b. order of individuals in the array is not random
    but doesn't matter because make_new_generation randomizes it
	
make_new_generation
    called from main
    takes size (integer) and pointer to population (integer) array
    returns modified population array, which represents the next generation of individuals
    does so by randomizing order of elements in population array
    then pairs adjacent elements and sends those pair to make_offspring function (twice)
    adds resulting offspring to new population (integer) array
    once done, replaces population array with contents of new population array

make_offspring
    called from make_new_generation
    takes two "parental" elements (integers) from population array
    returns one "offspring" element (integer) 
    so from two parental genotypes, function makes one random offspring genotype

fixation_check
    called from main
    takes size (integer) and pointer to population (integer) array
    checks to see if an allele has become fixed (frequency of 100%)
    in other words, checks to see if all elements of the population array are identical
    returns 1 if all integers identical (fixed), else returns 0
	
write_to_result
    called from main
    takes fixtime (int), popsizecounter (int), replicatecounter (int), pointer to raw_results (integer) array
    writes the fixtime for the completed simulation to the correct position in the raw_results array
    returns modified raw_results (integer) array

calculate_final_results
    called from main
    takes pointer to raw_results (integer) array and uses Replicates (integer constant) and Population_Sizes (integer constants) array
    for each population size, takes fixation times and calculates the average, median, standard deviation, minimum, and maximum
    populates the final_results global struct with those data and returns it
        
sorting_comparison
    called from calculate_final_results
    takes two integer elements and returns the difference between them
    required part of built-in qsort function, which is used to sort array to determine median fixation time
    
random_number
    called from make_new_generation and make_offspring
    takes two integers, the minimum and maximum desired values
    returns one random integer in range between min to max (including endpoints)

	
>>KEY VARIABLES & DATA STRUCTURES<<

final_results
    structure / global
    seven fields
        popsize (int) - the population size
        simulations (int) - the number of replicates (simulations) of that population size
        meanfixtime (float) - mean time to fixation across all simulations of that population size
        standarddeviation (float) - sample standard deviation of meanfixtime
        medianfixtime (float) - median time to fixation across all simulations of that population size
        fastestfix (int) - minimum time to fixation out of all simulations of that population size
        slowestfix (int) - maximum time to fixation out of all simulations of that population size
    each population size is a record called "pop[n]" where n is 0 to number of population sizes investigated
    
fixtime
    integer / main
    records how many generations have passed in a given simulations
    
Homo_Allele_One
    float (constant) / global
    Proportion of initial population that is homozygous
    for allele one (genotype frequency, in decimal format)

Homo_Allele_Two
    float (constant) / global
    Proportion of initial population that is homozygous
    for allele two (genotype frequency, in decimal format)
    
new_generation
    integer array / make_new_generation
    same size as population array
    temporarily holds new generation
    once completed, replaces existing population array
    
numpopsizes 
    integer / main, calculate_final_results
    the number of different population sizes to be simulated
    calculated by dividing the size of Population_Sizes array by the size of one integer 
	
popsizecounter
    integer / write_to_result
    counter showing which element in the Population_Sizes array is being investigated
    
population
    integer array / main
    represents the population of individuals that are being simulated
    each individual is represented by an integer
    11 = homozygous for allele one
    12 = heterozygous
    22 = homozygous for allele two
    
Population_Sizes
    integer array (constants) / global
    population sizes to investigate 
    sizes must be even numbers for results to be accurate
    
raw_results
    integer array / main
    stores the results (generations to fixation) for all simulations
    use one array for all results
    know that first [Replicates] entries are the 
    results for first Population_Size simulated and so forth
    
replicatecounter
    integer / write_to_result
    counter showing which replication is being done
    
Replicates 
    integer (constant) / global
    how many simulations to do of each population size
	
size
    integer / main
    particular population size being investigated
    pulled from Population_Sizes array

    
>>OTHER VARIABLES<<

all_ones
    integer (boolean) / fixation_check
    flag for whether allele "1" has become fixed in population
    1 = yes (fixed), 0 = no (not fixed)
    
all_twos
    integer (boolean) / fixation_check
    flag for whether allele "2" has become fixed in population
    1 = yes (fixed), 0 = no (not fixed)
    
avgfixtime
    float / calculate_final_results
    the average (mean) fixation time across all replicates of a given population size
    calculated as fixtotal/Replicates

dad_gamete  
    integer / make_new_generation
    the allele that the "father" passes to his offspring, either 1 or 2
    see also mom_gamete
    
deviation   
    float / calculate_final_results
    the standard deviation of avgfixtime
    calculated as sample standard deviation (n-1), not population standard deviation (n)
    
done    
    integer / main
    records time that program finishes, used to calculate runtime of program
    
elem_a
    integer / sorting_comparison
    one element of results array to be compared as part of built-in qsort function
    
elem_b    
    integer / sorting_comparison
    other element of results array to be compared as part of built-in qsort function
    
fixed
	integer (boolean)
	flag whether an allele has become fixed (1) or not (0)
    
fixtime 
    integer / calculate_final_results   
    an individual fixation time pulled from raw_results
    used to find minfix, maxfix, and calculate avgfixtime and deviation
    
fixtotal
    integer / calculate_final_results
    sum of all fixation times for all replicates of a given population size
    used to calculate average fixation time for that population size
    
hetero
    integer / create_initial_population
    the initial (generation zero) number of heterozygous individuals
    calculated as size - (homo1 + homo2)
    so if there is an "extra" individual, it becomes a heterozygote
    i.e., if size = 101 and Homo_Allele_One = Homo_Allele_Two = 0.4
    then you would end up with homo1 = 40, homo2 = 40, hetero = 21
    
homo1
    integer / create_initial_population
    the initial (generation zero) number of individuals homozygous for allele 1
    calculated as Homo_Allele_One * size

homo2
    integer / create_initial_population
    the initial (generation zero) number of individuals homozygous for allele 2
    calculated as Homo_Allele_Two * size
    
launch
    integer / main
    records time that program starts, used to calculate runtime of program
    
max 
    integer / random_number 
    the largest desired possible value of the random number to be generated

maxfix 
    integer / calculate_final_results
    the longest fixation time across all replicates of a given population size
    
medarray
    integer array / calculate_final_results
    array of all fixation times for a given population size (so size of array = Replicates)
    used to calculate the median fixation time
    
medfixtime  
    float / calculate_final_results 
    median fixation time for a given population size
    
min 
    integer / random_number
    the smallest desired possible value of the random number to be generated

minfix  
    integer / calculate_final_results   
    the shortest fixation time across all replicates of a given population size

mom_gamete  
    integer / make_new_generation
    the allele that the "mother" passes to her offspring, either 1 or 2
    see also dad_gamete
    
offspring   
    integer / make_new_generation
    the genotype of the offspring, so the combination of mom_gamete and dad_gamete
    either 11, 12, or 22
    
relem
    integer / make_new_generation
    randomly selected element of population array, used in shuffling its elements

relem_value
    integer / make_new_generation
    temporary placeholder for value of relem, used in shuffling elements of population array
    
results_place   
    integer / write_to_result
    the correct position in the raw_results array for the current fixation time result
    calculated as (popsizecounter * Replicates) + replicatecounter
    
rnum
    integer / random_number, make_new_generation
    random number
    in random_number function, it is generated using built-in rand function
    in make_new_generation function, it is obtained from random_number function 
        and used to decide which allele a heterozygote passes to its offspring
    
runtimeseconds
    double / main
    runtime of program, in seconds
    calculated as runtimeticks/CLOCKS_PER_SEC (constant from time.h library)

runtimeticks    
    double / main
    runtime of program, in clock ticks
    calculated as done - launch
    
seed    
    time / main    
    used to seed random number generator

totaldeviation  
    float / calculate_final_results 
    sum of the individual deviations from the average fixation time
    used to calculate the standard deviation
    
variance
    float / calculate_final_results
    the variance of the fixation times  
    used to calculate the standard deviation


>>COUNTER VARIABLES<<

ctrPop
    integer / main
    iterates across Population_Size array

ctrRep  
    integer / main
    iterates number of Replicates for each population size being simulated
    
ctrResults  
    integer / main
    iterates across final_results struct to report results
	
Outside the main function, loop counters are generally just called "ctr"


>>LICENSE<<

The MIT License (MIT)

Copyright (c) 2014-2015 David DeGusta

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

The Software is provided "as is", without warranty of any kind, express or
implied, including by not limited to the warranties of merchantability,
fitness for a partiular purpose and noninfringement. In no even shall the 
authors or copyright holders be liable for any claim, damages or other
liability, whether in an action of contract, tort or otherwise, arising from,
out of or in connection with the software or the use or other dealings in the
Software.
*/

/***************************
**** HEADERS ***************
***************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


/***************************
**** GLOBALS ***************
***************************/
#define Replicates 10//number of times each population size is simulated

int Population_Sizes [ ] = {10, 100, 1000}; //the population sizes to be simulated, each population size must be an EVEN NUMBER

#define Homo_Allele_One 0.3 //initial genotype frequency of homozygous for allele one ("11")
#define Homo_Allele_Two 0.3 //initial genotype frequency of homozygous for allele two ("22")

struct final_results //where the final results are stored
{
    int popsize;
    int simulations;
    float meanfixtime;
    float standarddeviation; // note: this is the sample SD, not the population SD
    float medianfixtime;
    int fastestfix;
    int slowestfix;
};
struct final_results pop[sizeof(Population_Sizes)/sizeof(int)];


/***************************
**** PROTOTYPES ************
***************************/
void create_initial_population (int size, int * population);
void make_new_generation (int size, int * population);
int make_offspring (int mom, int dad);
int fixation_check (int size, int * population);
void write_to_result (int fixtime, int popsizecounter, int replicatecounter, int * raw_results);
void calculate_final_results (int * raw_results);
int sorting_comparison (const void * elem_a, const void * elem_b);
int random_number (int min, int max);


/***************************
**** MAIN ******************
***************************/
int main(int argc, char **argv)
{
    clock_t launch = clock(); //start timer to measure how long program takes to run
    

    printf ("\n \n Status");
    printf ("\n ---------\n");
	
    time_t seed; //seed the random number generator
    srand(time(&seed)); //should only be done once per program run, thus doing it in main
	
    int numpopsizes; //number of entries in Population_Sizes (global array of int consts)
    numpopsizes = (sizeof(Population_Sizes))/(sizeof(int));
    int size; //placeholder for specific population size being used
	
    int fixtime; //counter to track number of generations
	
    int raw_results [(numpopsizes * Replicates)]; //array to hold fixtimes
	
    for (int ctrPop = 0; ctrPop < numpopsizes; ctrPop++)
        // iterate through each specified population size
    {
        for (int ctrRep = 0; ctrRep < Replicates; ctrRep++)
            // conduct the specified number of simulations (replications) of drift for the given population size
        {
            size = Population_Sizes[ctrPop]; // pull the population size being used
            int population [size]; // create an empty array for the population
            create_initial_population (size, population); // populate it
			
            //can comment out the next line to eliminate status report as program runs
            printf ("\n \tSimulating Run %d of Population Size %d ...", (ctrRep+1), size);
			
            fixtime = 0; // start at generation 0
			
            while (fixation_check (size, population) != 1)
                //make new generations one at a time until an allele is fixed
                {
                make_new_generation (size, population);
                fixtime +=1;
                }
			
                write_to_result (fixtime, ctrPop, ctrRep, raw_results); // store the fixation time for this run
        }
    }
    calculate_final_results (raw_results);
    clock_t done = clock(); //stop timer that is measuring how long program takes to run
    double runtimeticks = (done - launch); //calculate how long program took in clock ticks
    double runtimeseconds = runtimeticks/CLOCKS_PER_SEC; //convert program duration from clock ticks to seconds using system constant
    printf ("\n \n Finished! \n \n");
    printf ("\n \n Results");
    printf ("\n ----------\n");
    printf ("\n %d replicates for each population size", Replicates);
    printf ("\n With initial genotype frequencies of:");
    printf ("\n \t %.1f homozygous for 'allele 1'", Homo_Allele_One);
    printf ("\n \t %.1f homozygous for 'allele 2'", Homo_Allele_Two); 
    printf ("\n \t %.1f heterozygous \n \n", 1-(Homo_Allele_One + Homo_Allele_Two));
    
    for (int ctrResults = 0; ctrResults < numpopsizes; ctrResults ++)
    {
        printf ("\n For a Population Size of %d:", pop[ctrResults].popsize);
        printf ("\n \t \t Average Time to Fixation = %.0f generations (sd = %.0f)", pop[ctrResults].meanfixtime, pop[ctrResults].standarddeviation);
        printf ("\n \t \t Median Time to Fixation = %.0f generations", pop[ctrResults].medianfixtime);
        printf ("\n \t \t Range = %d to %d generations", pop[ctrResults].fastestfix, pop[ctrResults].slowestfix);
        printf ("\n \n \n");
    }
    printf ("\n");
    printf ("Runtime of program = about %.1f seconds (%.0f clock ticks) \n \n", runtimeseconds, runtimeticks);

return 0;
}


/************************************************
**** create_initial_population ******************
*************************************************/
void create_initial_population (int size, int * population)
    // takes population size (int) and pointer to population (int) array
    // populates population array with members
    // Genotype frequencies taken from the two global constants
    // Homo_Allele_One (homozygous "1", 11) and Homo_Allele_Two (homozygous "2", 22)
    // Order of elements in array doesn't matter because will be randomized in next step (make_new_generation)
{
    for (int ctr = 0; ctr < size; ctr++) 
        //fills array with zeroes
    {
        population[ctr] = 0;
    }
	
    int homo1; //number of homozygous "1" individuals (11's)
    int homo2; //number of homozygous "2" individuals (22's)
    int hetero; //number of heterozygote individuals (12's)
	
    homo1 = Homo_Allele_One * size;
    homo2 = Homo_Allele_Two * size;
    hetero = size - (homo1 + homo2); // so any "extra" individuals are heteros

    for (int ctr = 0; ctr < homo1; ctr++)
        //adds specified number [homo1] of homozygous 1 individuals (11's) to population
    {
        population[ctr] = 11;
    }
    
    for (int ctr = 0; ctr < homo2; ctr++)
        //adds specified number [homo2] of homozygous 2 individuals (22's) to population
    {
        population[ctr + homo1] = 22;
    }
    
    for (int ctr = 0; ctr < hetero; ctr++)
        //adds specified number [hetero] of heterozygote individuals (12's) to population
    {
        population[ctr + homo1 + homo2] = 12;
    }

return;
}


/************************************************
**** make_new_generation ************************
*************************************************/
void make_new_generation (int size, int * population)
    // takes population size (int) and pointer to population (int) array
    // returns modified population (int) array
    // does so by randomizing order of starting population (the population array)
    // then pairs up adjacent individuals
    // sends those pair to make_offspring function (twice)
    // adds result to new population (int) array
    // when finished, replace contents of original population array with those of new population array
{
    int ctr;
    for (ctr = 0; ctr < size - 1; ctr++) 
        //randomize order of individuals in the population array by shuffling the elements
        //this loop modified from code by Ben Pfaff (benpfaff.org/writings/clc/shuffle.html)
    {
        int relem = ctr + random_number(1,(size-1-ctr)); //relem is randomly selected element
        int relem_value = population[relem]; //relem_value is just a temporary placeholder for value of relem
        population[relem] = population[ctr]; //replace random element with current element
        population[ctr] = relem_value; //replace current element with random element
    }
    	
    int new_generation [size]; // new generation to populate
			
    for (int ctr = 0; ctr < size; ctr = ctr + 2)
        //iterate through existing population array, pair up individuals,
        //make offspring, and put them into new_generation array
    {
        new_generation [ctr] = make_offspring (population[ctr], population[ctr+1]);
        new_generation [(ctr+1)] = make_offspring (population[ctr], population[ctr+1]);
    }
	
    for (int ctr = 0; ctr < size; ctr++)
        //put new_generation individuals into population array
    {
        population [ctr] = new_generation [ctr];
    }
	
return;
}


/************************************************
**** make_offspring *****************************
*************************************************/
int make_offspring (int mom, int dad)
    // takes two "parental" elements (integers) from population array
    // returns one "offspring" element (integer) 
    // so from two parental genotypes, function makes one random offspring genotype
{
    int mom_gamete; //allele that mom passes on
    int dad_gamete; //allele that dad passes on
    int offspring; //combination of mom_gamete and dad_gamete
	
    int rnum; //random number to decide which allele gets passed on in heteros
	
//Figure out which gamete mom passes on
    if (mom == 11)
        //homozygote case
    {   mom_gamete = 1;}
    else if (mom == 22)
        //homozygote case
    {   mom_gamete = 2;}
    else
        //heterozygote case
	{
        rnum = random_number (1, 10);
        (rnum <= 5) ? (mom_gamete = 1) : (mom_gamete = 2);
    }
	
//Figure out which gamete dad passes on	
    if (dad == 11)
        //homozygote case
    {   dad_gamete = 1;}
    else if (dad == 22)
        //homozygote case
    {   dad_gamete = 2;}
    else
        //heterozygote case
    {
        rnum = random_number (1, 10);
        (rnum <= 5) ? (dad_gamete = 1) : (dad_gamete = 2);
    }

//Combine mom and dad gametes into one individual

    if ((mom_gamete == 1) && (dad_gamete == 1))
    {
        offspring = 11;
    }

    else if ((mom_gamete == 2) && (dad_gamete == 2))
    {
        offspring = 22;
    }
	
    else
    {
        offspring = 12;
    }

return offspring;
}


/************************************************
**** fixation_check *****************************
*************************************************/
int fixation_check (int size, int * population)
    // takes size (integer) and pointer to population (integer) array 
    // returns 1 if all elements identical (fixed), else returns 0
{	
    int all_ones = 1; // flag for whether "1" allele has become fixed. 1 = yes, 0 = no.
    int all_twos = 1; // flag for whether "2" allele has become fixed. 1 = yes, 0 = no.
	
    for (int ctr = 0; ctr < size; ctr++)
    {
        if (population [ctr] == 12)
            //if any heterozygote is found, then neither allele fixed
        {
            all_ones = 0;
            all_twos = 0;
        }
		
        else if (population [ctr] == 11)
            //if homozygote 1 is found, then "2" allele is not fixed
        {
            all_twos = 0;
        }	
		
        else if (population[ctr] == 22)
            //if homozygote 2 is found, then "1" allele is not fixed
        {
            all_ones = 0;
        }
		
        if ((all_ones == 0) && (all_twos == 0))
            //once you know neither allele is fixed, bail out and return false
        {
            return (0);
        }
			
    }
	
    if ((all_ones == 1) || (all_twos == 1))
        // should only reach here if an allele is fixed, but just to be sure, tested again
        {
            return (1);
        }

printf ("\n Error occured in fixation_check function. \n" ); //all cases should be handled by above, so if not, flag error
return (0);
}


/************************************************
**** write_to_result ****************************
*************************************************/
void write_to_result (int fixtime, int popsizecounter, int replicatecounter, int * raw_results)
    // takes fixtime (int), popsizecounter (int), replicatecounter (int), pointer to raw_results (int) array
    // writes fixtime to the right place in raw_results (int array) 
{	
    int results_place; //the element of the array where the results belong
    results_place = (popsizecounter * Replicates) + replicatecounter; //Replicates is global constant
    raw_results [results_place] = fixtime;
    printf (" Fixation happened in %d generations", fixtime); //comment this out if you don't want live status updates to appear
return;
}
	

/************************************************
**** calculate_final_results ********************
*************************************************/	
void calculate_final_results (int * raw_results)
    // takes pointer to raw_results (integer) array and uses Replicates (int constant) and Population_Sizes (int array constants)
    // for each population size, takes fixation times and calculates the average, median, standard deviation, minimum, and maximum
    // populates the final_results global struct with those data and returns it
{	
    int numpopsizes; //number of entries in Population_Sizes (global array of int consts)
    numpopsizes = (sizeof(Population_Sizes))/(sizeof(int));
	
    int fixtotal; //sum of all times to fixation for all replicates of a given population size
    float avgfixtime; //mean fixation time for all replicates of a given population size
    float deviation; //standard deviation of avgfixtime; n.b., uses sample standard deviation (n-1), not population version (n)
    float variance; //variance of the fixation times for a given population size, used to calculate the standard deviation
    float totaldeviation; //sum of the individual deviations of the fixation times from the average, used to calculate the standard deviation
    float medfixtime; //median fixation time for all replicates of a given population size
    int medarray [Replicates]; //array of fixation times for a given population size, used to calculate median fixation time
    int minfix; //shortest fixation time across all replicates of a given population size
    int maxfix; //longest fixation time across all replicates of a given population size
    int fixtime; //an individual fixation time pulled from raw results

    for (int ctrP = 0; ctrP < numpopsizes; ctrP++)
        //iterate through each population size
    {	
        fixtotal = 0; //resets to zero with start of each new population size
        minfix = raw_results [(ctrP * Replicates)]; //set to value of first replicate as starting point for a new population size
        maxfix = raw_results [(ctrP * Replicates)]; //set to value of first replicate as starting point for a new population size
			
        for (int ctrR = 0; ctrR < Replicates; ctrR++)
            //iterate through all the replicates of a given pop size and summarize results
        {	
            fixtime = raw_results[((ctrP * Replicates) + ctrR)];
            fixtotal += fixtime;
            minfix = (fixtime < minfix) ? (fixtime) : (minfix); //if fixtime is smallest so far, it becomes new minimum
            maxfix = (fixtime > maxfix) ? (fixtime) : (maxfix); //if fixtime is smallest so far, it becomes new minimum
            medarray [ctrR] = fixtime;
        }
        
        avgfixtime = ((float)(fixtotal))/Replicates; // calculate average
        totaldeviation = 0;
          
        for (int ctrSD = 0; ctrSD < Replicates; ctrSD++)
            //iterate through the medarray to calculate the total deviation
        {
            totaldeviation += (medarray[ctrSD] - avgfixtime) * (medarray[ctrSD] - avgfixtime); 
        }
          
        variance = totaldeviation/(Replicates-1); //the population variance; to change it to the sample variance, make denominator just (Replicates)
        deviation = sqrt (variance); //the standard deviation is the square root of the variance
        qsort (medarray, Replicates, sizeof(int), sorting_comparison); //sorts medarray from small to large to prepare for determining median value
        
        if ((Replicates % 2) == 0) //if number of replicates is an even number
        {
            medfixtime = ((float)(medarray[(Replicates/2)] + medarray[(Replicates/2)-1])/2); //average the two middle values to get median of even # of elements
        }
        else //if number of replicates is an odd number
        {
            medfixtime = medarray [(Replicates - 1) /2]; //takes the middle value to get the median of an odd # of elements
        }
        
        //write results for the given population size into final_results structure
        pop[ctrP].popsize = Population_Sizes[ctrP];
        pop[ctrP].simulations = Replicates;
        pop[ctrP].meanfixtime = avgfixtime;
        pop[ctrP].standarddeviation = deviation;
        pop[ctrP].medianfixtime = medfixtime;
        pop[ctrP].fastestfix = minfix;
        pop[ctrP].slowestfix = maxfix;
    }
return;
}


/************************************************
**** sorting_comparison ****************************
*************************************************/
int sorting_comparison (const void * elem_a, const void * elem_b)
    // takes two elements of array passed to it from qsort
    // returns difference between them, which is used by qsort to order elements
    // which is done in order to calculate the median value of an array of fixation times in the calculate_final_results function
{
return ( *(int*)elem_a - *(int*)elem_b );
}


/************************************************
**** random_number ******************************
*************************************************/
int random_number (int min, int max)
    // takes two non-negative integers
    // returns a random number in that range (endpoints included)
{
    int rnum;
    rnum = (min - 1);
    while (rnum < min)
	{
        rnum = (rand() % (max + 1));
    }
return rnum;
}