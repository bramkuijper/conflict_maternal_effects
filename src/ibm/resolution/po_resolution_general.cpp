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
const int Clutch = 200;

// mutation rates
double mu_s = 0;
double mu_q = 0;
double sdmu = 0;

// initial values
double initqS = 0;
double initqNS = 0;
double inits1 = 0;
double inits2 = 0;

// cost of z2 offspring
double gamma1 = 1.0;
double gamma2 = 1.0;

// cost of z1 offspring 
double beta1 = 1.0;
double beta2 = 1.0;

// maladaptation costs 
double c[2] = { 0, 0 };

// prob of environmental change
double sigma[2] = { 0, 0 };

// frequency of patch 1
double p = 0;

// dispersal probability
double d = 0;

// probability of nonlocal mating 
double l = 0;

// tally of dispersers
int Ndisp = 0;

// runtime for stats
time_t total_time; 

int generation = 0;

int seed = -1;

// skip the number of generations in output
// to prevent output files from becoming too large
int skip = 10;

// haploid individual
struct Individual
{
    // environment-dependent signalling probs
    double s[2][2];
    double qS[2];
    double qNS[2];

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
Individual Dispersers[Npatches * Npp * 10 * Clutch * 2];

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
    mu_s = atof(argv[1]); // mutation prob
    mu_q = atof(argv[2]); // mutation prob
    sdmu = atof(argv[3]); // mutational distribution
    c[0] = atof(argv[4]); 
    c[1] = atof(argv[5]); // survival slope
    p = atof(argv[6]); // freq envt 1
    d = atof(argv[7]); // dispersal 
    l = atof(argv[8]); // probability of mating locally
    sigma[0] = atof(argv[9]); //prob envt 1 changes
    sigma[1] = atof(argv[10]); // prob envt 2 changes
    gamma1 = atof(argv[11]); // production cost of z2 in envt e1
    gamma2 = atof(argv[12]); // production cost of z2 in envt e2
    beta1 = atof(argv[13]); // production cost of z1 in envt e1
    beta2 = atof(argv[14]); // production cost of z1 in envt e2
    initqS=atof(argv[15]);
    initqNS=atof(argv[16]);
    inits1=atof(argv[17]);
    inits2=atof(argv[18]);
}

void init_pop()
{
    // start the time
    total_time = time(NULL);

    // obtain a seed from current nanosecond count
	seed = get_nanoseconds();

    // set the seed to the random number generator

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
        MetaPop[i].envt = gsl_rng_uniform(r) > p;

        for (int j = 0; j < Npp; ++j)
        {
            for (int allele_i = 0; allele_i < 2; ++allele_i)
            {
                MetaPop[i].locals[j].s[allele_i][0] = inits1/2;
                MetaPop[i].locals[j].s[allele_i][1] = inits2/2;
                    
                MetaPop[i].locals[j].qNS[allele_i] = initqNS/2;
                MetaPop[i].locals[j].qS[allele_i] = initqS/2;
            }

            MetaPop[i].locals[j].phen = MetaPop[i].envt;
        }
    }
}

// mutate an allele
double MutS(double h)
{
    h+=gsl_rng_uniform(r) < mu_s ? gsl_ran_gaussian(r, sdmu) : 0;

    if (h < 0)
    {
        h = 0;
    }
    else if (h > 0.5)
    {
        h = 0.5;
    }
    return(h);
}

double MutQ(double h)
{
    h+=gsl_rng_uniform(r) < mu_q ? gsl_ran_gaussian(r, sdmu) : 0;

    if (h < 0)
    {
        h = 0;
    }
    else if (h > 0.5)
    {
        h = 0.5;
    }
    return(h);
}

// allocate a kid and give it genes
void Create_Kid(Individual &mother, Individual &father, Individual &Kid, bool current_envt)
{
    for (int i = 0; i < 2; ++i)
    {
        Kid.s[0][i] = MutS(mother.s[gsl_rng_uniform_int(r,2)][i]);
        Kid.s[1][i] = MutS(father.s[gsl_rng_uniform_int(r,2)][i]);
    }


        
    Kid.qS[0] = MutQ(mother.qS[gsl_rng_uniform_int(r,2)]);
    Kid.qS[1] = MutQ(father.qS[gsl_rng_uniform_int(r,2)]);

    Kid.qNS[0] = MutQ(mother.qNS[gsl_rng_uniform_int(r,2)]);
    Kid.qNS[1] = MutQ(father.qNS[gsl_rng_uniform_int(r,2)]);

    double s = mother.s[0][current_envt] + mother.s[1][current_envt];

    assert(s >= 0);
    assert(s <= 1.0);
    double qS = Kid.qS[0] + Kid.qS[1];
    double qNS = Kid.qNS[0] + Kid.qNS[1];

    double prob = s * qS + (1-s) * qNS;

    // Kid.phenotype should be 0 in envt 1  and 1 in envt 2
    Kid.phen= gsl_rng_uniform(r) < prob ? 0 : 1;
}

// mate and create kids across all patches...
void make_juveniles()
{
    Ndisp = 0;

    // save current envt before changing it
    bool current_envt;

    // total amount of resources available to
    // mother to make offspring
    double resources;

    // cost of individual offspring
    double cost;
    
    double current_gamma = 0;
    double current_beta = 0;
    
    
    int father;
    int father_patch;

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
                // maternal signal given to offspring
                // depends on the current 
                // environment (before any env'tal change)
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
                        // too few resources, no offspring produced
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
                    //
                    // we'll deal with dispersing offspring later
                    if (Kid.phen == MetaPop[i].envt 
                            || gsl_rng_uniform(r) < 1.0 - c[MetaPop[i].envt])
                    {
                        MetaPop[i].phils[MetaPop[i].Nkids++] = Kid;
                    }
                }
            } 
        }//Npp
    }//Npatches

    assert(Ndisp < Npatches * Npp * 10 * Clutch * 2);
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
            // store randomly chosen disperser from the stack of dispersers
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

    // all dispersers should have dispersed...
    assert(Ndisp==0);
}

void write_data_headers()
{
    DataFile << "generation;s1;s2;qNS;qS;vars1;vars2;varqNS;varqS;prob1;prob2;p1;" << endl;
}

void write_data()
{
    double means0 = 0;
    double means1 = 0;
    double meanqNS = 0;
    double meanqS = 0;
    double sss0 = 0;
    double sss1 = 0;
    double ssqNS = 0;
    double ssqS = 0;
    double p1 = 0;

    double prob0 = 0;
    double prob1 = 0;

    for (int i = 0; i < Npatches; ++i)
    {
        for (int j = 0; j < Npp; ++j)
        {
            means0 += MetaPop[i].locals[j].s[0][0] + MetaPop[i].locals[j].s[1][0];
            means1 += MetaPop[i].locals[j].s[0][1] + MetaPop[i].locals[j].s[1][1];

            meanqS += MetaPop[i].locals[j].qS[0]+ MetaPop[i].locals[j].qS[1];
            meanqNS += MetaPop[i].locals[j].qNS[0]+ MetaPop[i].locals[j].qNS[1];
            
            sss0 += (MetaPop[i].locals[j].s[0][0] + MetaPop[i].locals[j].s[1][0]) * (MetaPop[i].locals[j].s[0][0] + MetaPop[i].locals[j].s[1][0]);
            sss1 += (MetaPop[i].locals[j].s[0][1] + MetaPop[i].locals[j].s[1][1]) *  (MetaPop[i].locals[j].s[0][1] + MetaPop[i].locals[j].s[1][1]);
            ssqS += (MetaPop[i].locals[j].qS[0]+ MetaPop[i].locals[j].qS[1]) * (MetaPop[i].locals[j].qS[0]+ MetaPop[i].locals[j].qS[1]);
            ssqNS += (MetaPop[i].locals[j].qNS[0]+ MetaPop[i].locals[j].qNS[1]) * (MetaPop[i].locals[j].qNS[0]+ MetaPop[i].locals[j].qNS[1]);

            prob0 += (MetaPop[i].locals[j].s[0][0] + MetaPop[i].locals[j].s[1][0]) * (MetaPop[i].locals[j].qS[0]+ MetaPop[i].locals[j].qS[1]) +
                            (1.0 - (MetaPop[i].locals[j].s[0][0] + MetaPop[i].locals[j].s[1][0])) * (MetaPop[i].locals[j].qNS[0]+ MetaPop[i].locals[j].qNS[1]);

            prob1 += (MetaPop[i].locals[j].s[0][1] + MetaPop[i].locals[j].s[1][1]) * (MetaPop[i].locals[j].qS[0]+ MetaPop[i].locals[j].qS[1]) +
                            (1.0 - (MetaPop[i].locals[j].s[0][1] + MetaPop[i].locals[j].s[1][1])) * (MetaPop[i].locals[j].qNS[0]+ MetaPop[i].locals[j].qNS[1]);

        }
            
        p1 += !MetaPop[i].envt;
    }
    
    means0 /= (Npatches * Npp);
    means1 /= (Npatches * Npp);
    meanqS /= (Npatches * Npp);
    meanqNS /= (Npatches * Npp);
    prob0 /= Npatches * Npp;
    prob1 /= Npatches * Npp;

    p1 /= Npatches;
        
    DataFile << generation << ";" << means0 << ";" << means1 << ";"  << meanqNS << ";" << meanqS << ";"
                                << sss0 / (Npatches * Npp)  - means0 * means0 << ";" 
                                << sss1 / (Npatches * Npp)  - means1 * means1 << ";"  
                                << ssqNS / (Npatches * Npp)  - meanqNS * meanqNS << ";" 
                                << ssqS / (Npatches * Npp)  - meanqS * meanqS << ";"  
                                << prob0 << ";"
                                << prob1 << ";"
                                << p1 << ";"
                                << endl;
}

void write_parameters()
{
    DataFile << endl << endl << "patch;" << Npatches << endl
                << "npp;" << Npp << endl
                << "numgen;" << numgen << endl
                << "mu_s;" << mu_s << endl
                << "mu_q;" << mu_q << endl
                << "sdmu;" << sdmu << endl
                << "gamma1;" << gamma1 << endl
                << "gamma2;" << gamma2 << endl
                << "beta1;" << beta1 << endl
                << "beta2;" << beta2 << endl
                << "c1;" << c[0] << endl
                << "c2;" << c[1] << endl
                << "s12;" << sigma[0] << endl
                << "s21;" << sigma[1] << endl
                << "p;" << p << endl
                << "l;" << l << endl
                << "seed;" << seed << endl
                << "d;" << d << endl
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
