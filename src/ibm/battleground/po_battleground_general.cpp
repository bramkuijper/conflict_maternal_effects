#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <cmath>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "auxiliary.h"


//#define NDEBUG

using namespace std;

// random number generator 
// see http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generation.html#Random-Number-Generation 
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *r; // gnu scientific rng 

const int Npatches = 2000; 
const int Npp = 1; // individuals per patch
const int numgen = 20000;
const int Clutch = 100;

// mutation rates
double mu = 0;
double sdmu = 0;

// maladaptation costs 
double c[2] = { 0, 0 };

// prob of environmental change
double sigma[2] = { 0, 0 };

// frequency of patch 1
double p = 0;

// dispersal probability
double d = 0;

// probability of local mating
double l = 0;

// cost of z2 offspring
double gamma1 = 1.0;
double gamma2 = 1.0;

// cost of z1 offspring 
double beta1 = 1.0;
double beta2 = 1.0;

// tally of dispersers
int Ndisp = 0;

// runtime for stats
time_t total_time; 

int generation = 0;

int seed = -1;

bool maternal_expression = true;

// skip the number of generations in output
// to prevent output files from becoming too large
int skip = 10;

// haploid individual
struct Individual
{
    // environment-dependent signalling probs
    double s[2];

    bool phen;
};

struct Patch
{
    Individual locals[Npp]; // all the local breeders

    // philopatric offspring
    Individual phils[2 * Npp * 10 * Clutch];     
    // (note that dispersing offspring
    // ends up in global pool, see below)

    // total number of kids 
    int Nkids; 

    // variable that allows for correct
    // rounding of the number of immigrants
    // per patch (see below)
    int immigrant_bonus;

    // local environmental state
    bool envt;
};

// generate the complete population
Patch MetaPop[Npatches];
Individual Dispersers[Npatches * Npp * 10 * Clutch];

// give the outputfile a unique name
// by using the create_filename function (see 
// "bramauxiliary.h")
string filename("sim_conflict");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

// initialize the command line arguments to vary 
// the parameters
void init_arguments(int argc, char *argv[])
{
    mu = atof(argv[1]); // mutation prob
    sdmu = atof(argv[2]); // mutational distribution
    c[0] = atof(argv[3]); 
    c[1] = atof(argv[4]); // survival slope
    p = atof(argv[5]); // initial freq envt 1
    d = atof(argv[6]); // dispersal 
    l = atof(argv[7]); // probability of mating locally
    sigma[0] = atof(argv[8]); //prob envt 1 changes
    sigma[1] = atof(argv[9]); // prob envt 2 changes
    gamma1 = atof(argv[10]); // production cost of z2 in envt e1
    gamma2 = atof(argv[11]); // production cost of z2 in envt e2
    beta1 = atof(argv[12]); // production cost of z1 in envt e1
    beta2 = atof(argv[13]); // production cost of z1 in envt e2
    maternal_expression = atoi(argv[14]); // dispersal at t=0
}

void init_pop()
{
    // start the time
    total_time = time(NULL);

    // obtain a seed from current nanosecond count
	seed = get_nanoseconds();

    // set the seed to the random number generator
    // stupidly enough, this can only be done by setting
    // a shell environment parameter

    // set up the random number generators
    // (from the gnu gsl library)
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);

    // go through all patches
    for (int i = 0; i < Npatches; ++i)
    {
        // initialize the population with an environment at the
        // desired frequency
        MetaPop[i].envt = gsl_rng_uniform(r) < p;

        for (int j = 0; j < Npp; ++j)
        {
            for (int env = 0; env < 2; ++env)
            {
                MetaPop[i].locals[j].s[env] = 0.5;
            }

            MetaPop[i].locals[j].phen = MetaPop[i].envt;
        }
    }
}

// mutate an allele
double Mut(double h)
{
    h+=gsl_rng_uniform(r) < mu ? gsl_ran_gaussian(r, sdmu) : 0;
    h = h < 0 ? 0 : h > 1.0 ? 1.0 : h;
    return(h);
}

// allocate a kid and give it genes
void Create_Kid(Individual &mother, Individual &father, Individual &Kid, bool current_envt)
{
    for (int i = 0; i < 2; ++i)
    {
        Kid.s[i] = gsl_rng_uniform(r) < 0.5 ? Mut(mother.s[i]) : Mut(father.s[i]);
    }

    // determine who is expressing the trait
    double prob = maternal_expression ? mother.s[current_envt] : Kid.s[current_envt];

    // prob is the probability that the kid's phenotype is 0, rather than 1
    Kid.phen= gsl_rng_uniform(r) < prob ? 0 : 1;
}

// mate and create kids across all patches...
void make_juveniles()
{
    Ndisp = 0;

    // save current envt before changing it
    bool current_envt;

    // amount of resources for each parent
    double resources;

    // cost of an individual offspring
    double cost;

    int father;
    int father_patch;

    double current_gamma = 0;
    double current_beta = 0;

    // generate offspring on each patch
    for (int i = 0; i < Npatches; ++i)
    {
        // reset local offspring counter
        MetaPop[i].Nkids = 0;

        // variable that takes into account
        // any beyond-the-minimum number of immigrants
        // (see later)
        MetaPop[i].immigrant_bonus = 0;

        // change the environment
        // save the previous one as we still need to produce offspring
        current_envt = MetaPop[i].envt;

        if (gsl_rng_uniform(r) < sigma[MetaPop[i].envt])
        {
            MetaPop[i].envt = !MetaPop[i].envt;
        }

        current_gamma = current_envt ? gamma2 : gamma1;
        current_beta = current_envt ? beta2 : beta1;

        for (int j = 0; j < Npp; ++j)
        {

            // resources should be sufficient to at least
            // make Nclutch offspring
            resources = Clutch;

            resources *= current_envt ? gamma2 + beta2 : gamma1 + beta1;
                
            // create kids and reduce parental resources
            while (resources > 0)
            {
                // sample random father
                father = gsl_rng_uniform_int(r, Npp);

                if (gsl_rng_uniform(r) < l)
                {
                    // find random father within patch
                    father_patch = i;
                }
                else
                {
                    // find random father within remote patch
                    father_patch = gsl_rng_uniform_int(r, Npatches);
                }

                Individual Kid; 

                // create kid
                // maternal signal etc dependent on the current environment (before any change)
                Create_Kid(MetaPop[i].locals[j],
                        MetaPop[father_patch].locals[father],
                        Kid, 
                        current_envt);
                
                // cost of producing this offspring
                cost = Kid.phen ? current_gamma : current_beta;

                // offspring too costly for the few resources
                // that remain
                if (resources < cost)
                {
                    // calculate the probability that offspring
                    // will be produced anyway
                    if (gsl_rng_uniform(r) > resources / cost)
                    {
                        // to few resources, no offspring produced
                        break;
                    }

                    resources = 0;
                }
                else
                {
                    resources -= cost;
                }

                // disperse or not
                if (gsl_rng_uniform(r) < d)
                {
                    Dispersers[Ndisp++] = Kid;
                }
                else
                {
                    // already perform survival of philopatric offspring
                    // survival dependent on the novel envt (after change)
                    if (Kid.phen == MetaPop[i].envt || gsl_rng_uniform(r) < 1.0 - c[MetaPop[i].envt])
                    {
                        MetaPop[i].phils[MetaPop[i].Nkids++] = Kid;
                    }
                }
            } //clutch
        }//Npp
    }//Npatches

    assert(Ndisp < Npatches * Npp * Clutch * 10);
}

// replacement of adults with juveniles
void replace_adults()
{
    int rand_disperser;
    int surviving_immigrants;
    int arriving_immigrants;

    // okay, we have Ndisp/Npatches dispersers per patch
    int dispersers_per_patch = floor((double) Ndisp / Npatches);
    //cout << "Ndisp: " << Ndisp << ", Npatches: " << Npatches << " disp per patch: " << dispersers_per_patch << endl;

    // however, we need to correctly round this rational number
    // to the actual number of immigrants for 
    // a given patch. To this end, 
    // we just randomly distribute the dispersers that remain after the
    // previous rounding over the different patches
    for (int i = 0; i < Ndisp - Npatches * dispersers_per_patch; ++i)
    {
        // randomly picked patch receives additional immigrant
        MetaPop[gsl_rng_uniform_int(r,Npatches)].immigrant_bonus++;
    }

    
    // now replace local breeders on each patch
    for (int i = 0; i < Npatches; ++i)
    {
        if (dispersers_per_patch + MetaPop[i].immigrant_bonus > 0)
        {
            assert(Ndisp > 0);
        }

        arriving_immigrants = dispersers_per_patch + MetaPop[i].immigrant_bonus;

        Individual incoming[arriving_immigrants];

        surviving_immigrants = 0;
        
        // first perform patch-dependent survival of immigrants
        for (int j = 0; j < arriving_immigrants; ++j)
        {
            // store disperser
            rand_disperser = gsl_rng_uniform_int(r, Ndisp);

            // does disperser survive in the local envt?
            if (Dispersers[rand_disperser].phen == MetaPop[i].envt || gsl_rng_uniform(r) < 1.0 - c[MetaPop[i].envt])
            {
                incoming[surviving_immigrants++] = Dispersers[rand_disperser];
            }
            // remove disperser from the global pool
            Dispersers[rand_disperser] = Dispersers[Ndisp-1];
            --Ndisp;
        }

        for (int j = 0; j < Npp; ++j)
        {
            assert((double) surviving_immigrants / (surviving_immigrants + MetaPop[i].Nkids) >= 0.0 && (double) surviving_immigrants / (surviving_immigrants + MetaPop[i].Nkids) <= 1.0);

            if (gsl_rng_uniform(r) < (double) surviving_immigrants / (surviving_immigrants + MetaPop[i].Nkids))
            {
                int rand_disp = gsl_rng_uniform_int(r,surviving_immigrants);
                MetaPop[i].locals[j] = incoming[rand_disp];
                incoming[rand_disp] = incoming[surviving_immigrants - 1];
                --surviving_immigrants;
            }
            else
            {
                int rand_phil = gsl_rng_uniform_int(r,MetaPop[i].Nkids);
                MetaPop[i].locals[j] = MetaPop[i].phils[rand_phil];
                MetaPop[i].phils[rand_phil] = MetaPop[i].phils[MetaPop[i].Nkids-1];
                --MetaPop[i].Nkids;
            }
        }
    }

    assert(Ndisp==0);
}

void write_data_headers()
{
    DataFile << "generation;s1;s2;vars1;vars2;p1;" << endl;
}

void write_data()
{
    double means0 = 0;
    double means1 = 0;
    double sss0 = 0;
    double sss1 = 0;

    double p2 = 0;

    for (int i = 0; i < Npatches; ++i)
    {
        for (int j = 0; j < Npp; ++j)
        {
            means0 += MetaPop[i].locals[j].s[0];
            means1 += MetaPop[i].locals[j].s[1];
            sss0 += MetaPop[i].locals[j].s[0] * MetaPop[i].locals[j].s[0];
            sss1 += MetaPop[i].locals[j].s[1] * MetaPop[i].locals[j].s[1];
        }

        p2 += MetaPop[i].envt;
    }

    p2 /= Npatches;
    
    means0 /= (Npatches * Npp);
    means1 /= (Npatches * Npp);
        
    DataFile << generation << ";" << means0 << ";" << means1 << ";" 
                                << sss0 / (Npatches * Npp)  - means0 * means0 << ";" 
                                << sss1 / (Npatches * Npp)  - means1 * means1 << ";"  
                                << 1.0 - p2 << ";" <<
                                endl;
}

void write_parameters()
{
    DataFile << endl << endl << "patch;" << Npatches << endl
                << "npp;" << Npp << endl
                << "numgen;" << numgen << endl
                << "mu;" << mu << endl
                << "sdmu;" << sdmu << endl
                << "gamma1;" << gamma1 << endl
                << "gamma2;" << gamma2 << endl
                << "beta1;" << beta1 << endl
                << "beta2;" << beta2 << endl
                << "c0;" << c[0] << endl
                << "c1;" << c[1] << endl
                << "s12;" << sigma[0] << endl
                << "s21;" << sigma[1] << endl
                << "p;" << p << endl
                << "l;" << l << endl
                << "d;" << d << endl
                << "expression;" << (maternal_expression ? "mother" : "offspring") << endl
                << "runtime;" << total_time << endl;
}


int main(int argc, char * argv[])
{
    init_arguments(argc,argv);
    init_pop();

    write_data_headers();

    for (generation = 0; generation < numgen; ++generation)
    {
        make_juveniles();

        if (generation % skip == 0)
        {
            write_data();
        }

        replace_adults();
    }

    write_data();
    write_parameters();
}
