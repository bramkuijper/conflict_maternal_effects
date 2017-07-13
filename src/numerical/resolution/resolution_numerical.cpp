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


// give the outputfile a unique name
// by using the create_filename function (see 
// "bramauxiliary.h")
string filename("iter_conflict");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

double bound(double const val)
{
    return(val < 0 ? 0 : (val > 1.0 ? 1.0 : val));
}
double ent(double p)
{
    if (p == 0)
    {
        return(0);
    }
    if (p==1.0)
    {
        return(0);
    }

    return(- p * (log2(p)+(1-p) * log2(1-p)));
}
double siginfo(double p1, double p2)
{
    return(1 - (((p1 + p2)/2) * ent(p1/ ( p1 + p2)) + (1 - ((p1 + p2)/2)) * ent((1-p1)/ ((1-p1)+(1-p2)))));
}

int main(int argc, char **argv)
{
    double s1, s1tplus1;
    double s2, s2tplus1;
    double qNS, qNStplus1;
    double qS, qStplus1;

    double p1 = atof(argv[1]);
    double n = atof(argv[2]);
    double d = atof(argv[3]);
    double sigma12 = atof(argv[4]);
    double sigma21 = atof(argv[5]);
    double gam1 = atof(argv[6]);
    double gam2 = atof(argv[7]);
    double beta1 = atof(argv[8]);
    double beta2 = atof(argv[9]);
    double c1 = atof(argv[10]);
    double c2 = atof(argv[11]);
    double l = atof(argv[12]);

    // eigenvector related thiings
    double ev, v1, v2, u1;
    double dWds1, dWds2, dWdqS, dWdqNS;
    double dWds1offspring, dWds2offspring, dWds1mother, dWds2mother;

    double s1mother, s2mother;
    double s1offspring, s2offspring;

    DataFile << "s1;s2;z1_e1_off;z1_e2_off;z1_e1_mom;z1_e2_mom;zqS;qNS;p1;l;d;n;gam1;gam2;beta1;beta2;c1;c2;sigma12;sigma21;ev;u1;v1;v2;siginfo;" << endl;

            // first calculate battleground
            s1 = 0.5;
            s2 = 0.5;

            for (size_t iter = 0; iter < 10000000; ++iter)
            {
                // first calculate eigenvectors et al 
                EIGENVECS_BATTLE

                // then evaluate selection gradients mother
                SELGRADS_MOM

                // update values
                s1tplus1 = s1 + 0.01 * dWds1mother;
                s2tplus1 = s2 + 0.01 * dWds2mother;

                if (s1tplus1 < 0)
                {
                    s1tplus1 = 0;
                } else if (s1tplus1 > 1.0)
                {
                    s1tplus1 = 1.0;
                }
                
                if (s2tplus1 < 0)
                {
                    s2tplus1 = 0;
                } else if (s2tplus1 > 1.0)
                {
                    s2tplus1 = 1.0;
                }

                if (
                        (fabs(s2tplus1 - s2) < 1e-08 &&
                        fabs(s1tplus1 - s1) < 1e-08)
                        ||
                        iter == 10000000 - 1
                )
                {
                    s2 = s2tplus1;
                    s1 = s1tplus1;
                    break;
                }
                
                //cout << "s1mom: " << s1 << " s2mom: " << s2 << endl;

                assert(isnan(s2) == 0);
                assert(isnan(s1) == 0);
                
                s1 = s1tplus1;
                s2 = s2tplus1;
            }

            s1mother = s1;
            s2mother = s2;

            // now the offspring

            s1 = 0.5;
            s2 = 0.5;

            
            for (size_t iter = 0; iter < 10000000; ++iter)
            {
                // first calculate eigenvectors et al 
                EIGENVECS_BATTLE

                // then evaluate selection gradients offspring
                SELGRADS_OFF

                // update values
                s1tplus1 = s1 + 0.01 * dWds1offspring;
                s2tplus1 = s2 + 0.01 * dWds2offspring;

                if (s1tplus1 < 0)
                {
                    s1tplus1 = 0;
                } else if (s1tplus1 > 1.0)
                {
                    s1tplus1 = 1.0;
                }
                
                if (s2tplus1 < 0)
                {
                    s2tplus1 = 0;
                } else if (s2tplus1 > 1.0)
                {
                    s2tplus1 = 1.0;
                }

                if (
                        (fabs(s2tplus1 - s2) < 1e-08 &&
                        fabs(s1tplus1 - s1) < 1e-08)
                        ||
                        iter == 10000000 - 1
                )
                {
                    s2 = s2tplus1;
                    s1 = s1tplus1;
                    break;
                }
                //cout << "s1off: " << s1 << " s2off: " << s2 << endl;

                assert(isnan(s2) == 0);
                assert(isnan(s1) == 0);
                
                s2 = s2tplus1;
                s1 = s1tplus1;
            }

            s1offspring = s1;
            s2offspring = s2;



            s1 = 0.9;
            s2 = 0.1;
            qS = 0.9;
            qNS = 0.1;

            //cout << "optima done" << endl;

            for (size_t iter = 0; iter < 1000000; ++iter)
            {
                if (iter % 1000 == 0)
                {
                //cout << iter << " s1: " << s1 << " s2: " << s2 << " qNS: " << qNS << " qS: " << qS << endl;
                }
                // first calculate eigenvectors et al 
                EIGENVECS_RESOLUTION

                // then evaluate selection gradients mother
                SELGRADS_RESOLUTION

                // update values
                s1tplus1 = s1 + 0.01 * dWds1;
                s2tplus1 = s2 + 0.01 * dWds2;
                qNStplus1 = qNS + 0.01 * dWdqNS;
                qStplus1 = qS + 0.01 * dWdqS;

                s1tplus1 = bound(s1tplus1);
                s2tplus1 = bound(s2tplus1);
                qStplus1 = bound(qStplus1);
                qNStplus1 = bound(qNStplus1);

                if (
                        (
                        fabs(s2tplus1 - s2) < 1e-08 &&
                        fabs(s1tplus1 - s1) < 1e-08 &&
                        fabs(qNStplus1 - qNS) < 1e-08 &&
                        fabs(qStplus1 - qS) < 1e-08
                       )
                        ||
                        iter == 10000000 - 1
                )
                {
                    s2 = s2tplus1;
                    s1 = s1tplus1;
                    qNS = qNStplus1;
                    qS = qStplus1;
                    DataFile <<  s1 << ";" << s2 << ";" << s1offspring << ";" << s2offspring << ";" << s1mother << ";" << s2mother << ";" << qS << ";" << qNS << ";" << p1 << ";" << l << ";" << d << ";" << n << ";" << gam1 << ";" << gam2 << ";" << beta1 << ";" << beta2 << ";" << c1 << ";" << c2 << ";" << sigma12 << ";" << sigma21 << ";" << ev << ";" << u1 << ";" << v1 << ";" << v2 << ";" << siginfo(s1,s2) << ";" << endl;
                    break;
                }

                //cout << "s1: " << s1 << " s2: " << " qNS: " << qNS << " qS: " << qS << endl;

                assert(isnan(s2) == 0);
                assert(isnan(s1) == 0);
                assert(isnan(qNS) == 0);
                assert(isnan(qS) == 0);
                
                s2 = s2tplus1;
                s1 = s1tplus1;
                qNS = qNStplus1;
                qS = qStplus1;
            }
}
