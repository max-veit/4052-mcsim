#include <iostream>
#include <cstdlib>
#include <time.h>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <fstream>

using namespace std;

// Define pi
const double pi = 3.14159;

const double deg_per_rad = 180.0/pi;


int main(int argc, char** argv)
{
    const char *usagemsg = "Usage: %s [-n numevents] [-l plength] [-w pwidth]"
        "[-h height] [-x dist | -i inputfile] [-o outputfile]\n"
        "\tWhere plength and pwidth are the paddles' lengths and widths,"
        "\n\theight is the height of the top paddle off the ground. The paddle"
        "\n\tdistances must either be specified on the command line, or in"
        "\n\tthe specified inputfile as a newline-separated list of paddle"
        "\n\tx-positions, beginning with the number of paddles."
        "\n\tIf the output file is not specified, it defaults to stdout.";

    //;-------------------------------------------
    //; Physical Properties of System
    // To be read from the command-line arguments, with defaults in case
    // the are left unspecified.

    //;---------------------------------------
    //; Panel Dimensions
    double L = 70;  //Panel length (cm) along north-south (y-axis)
    double W = 19;  //Panel width (cm) along east-west (x-axis)

    //;-----------------------------------------
    //  Panel Locations
    // x is in the direction of east and west. y is in the direction
    // of north and south. z is the height above the ground.

    // Upper Panel Center Location
    double x0, y0, z0;
    x0 = 0;        // in (cm)
    y0 = 0;        // in (cm)
    z0 = 124.5;      // in (cm) Height above ground

    //  Lower Panel Center Location
    double x1, y1, z1;
    x1 = 136;        // in (cm), distance from upper panel
    y1 = 0;         // in (cm) Keep panels in a line
                    // so that there is no deviation in
                    // the north-south direction.
    z1 = 0;         // in (cm) on ground

    // Provisions for specifying multiple lower panels
    int pidx;
    int n_lpanels = 0;
    double *lower_x = NULL;

    ifstream fin;
    ostream *fout = NULL;

    // Simulation parameters
    long trial_max = 10000000; //10 million trials by default

    // Parse the arguments
    char* prog_name = argv[0];
    argv++;
    argc--;

    while (argc > 0) {
        if ((strcmp(argv[0], "-n") == 0) && argc > 1) {
            trial_max = atoi(argv[1]);
        } else if ((strcmp(argv[0], "-l") == 0) && argc > 1) {
            L = atof(argv[1]);
        } else if ((strcmp(argv[0], "-w") == 0) && argc > 1) {
            W = atof(argv[1]);
        } else if ((strcmp(argv[0], "-h") == 0) && argc > 1) {
            z0 = atof(argv[1]);
        } else if ((strcmp(argv[0], "-x") == 0) && argc > 1) {
            x1 = atof(argv[1]);
        } else if ((strcmp(argv[0], "-i") == 0) && argc > 1) {
            fin.open(argv[1], ifstream::in);
        } else if ((strcmp(argv[0], "-o") == 0) && argc > 1) {
            fout = new ofstream(argv[1], ofstream::out | ofstream::trunc);
        }
        // Advance to the next set of arguments.
        argc -= 2;
        argv += 2;
    }

    // Parse the input file, if it exists.
    if (fin.is_open()) {
        fin >> n_lpanels;
        lower_x = (double*)malloc(n_lpanels*sizeof(double));
        for (pidx = 0; pidx < n_lpanels; pidx++) {
            fin >> lower_x[pidx];
        }
        fin.close();
    } else {
        // Just use the single panel x specified/defaulted to.
        n_lpanels = 1;
        lower_x = (double*)malloc(n_lpanels*sizeof(double));
        lower_x[0] = x1;
    }


//;----------------------------------------------------
// Random number generator

    //;------------------------------------
    // Choose theta_max (p_max)

        double p_max;

        p_max = pi / 2;

    //;-------------------------------------
    // Set seed

        srand(time(0));



    //;---------------------------
    // Start Simulation

    long trial_num = 0;
    // Different counts for each panel. Make sure they're initialized to 0.
    int *counts = (int*)calloc(n_lpanels, sizeof(int));


    double x_upper, y_upper;    // Random spot on upper panel
    double phi, p;              // Random orientation of muon
    double x, y;                // Where muon lands on ground

    // Luckily, setting the bytes to zero works for doubles as well.
    double *sum_p =
        (double*)calloc(n_lpanels, sizeof(double));   //Radians
    double *sum_p_squared =
        (double*)calloc(n_lpanels, sizeof(double));   //Radians


    // The nesting order of the for-loops doesn't really matter, except
    // maybe for cache efficieny (relevant only for a large number of panels)
    for (pidx = 0; pidx < n_lpanels; pidx++)
    {
    for (trial_num = 0; trial_num < trial_max; trial_num++)
    {
        //;----------------------------------------
        // Choose random spot on upper panel;

        // double x_upper, y_upper;

        x_upper = W*rand()/RAND_MAX - W/2;
        y_upper = L*rand()/RAND_MAX - L/2;

        //;----------------------------------------
        // Choose random direction to travel with phi in the
        // azimuthal direction (measured from east), and p in the zenithal direction

        // double phi, p;

        phi = 2*pi*rand()/RAND_MAX;
        p = acos( pow ( ( 1 - (1 - pow(cos(p_max), 3.0) ) * rand()/RAND_MAX ) , 1.0 / 3) );

        //;-----------------------------------------------
        // Find where muon from top panel lands.
            //Required information:
            //Azimuthal direction 'phi'
            //Zenithal direction 'p'
            //Height above ground 'z0'
            //x_location 'x_upper'
            //y_location 'y_upper'

        //double x, y;    //Location on ground where muon passes through


        x = z0 * tan(p) * cos(phi) + x_upper;
        y = z0 * tan(p) * sin(phi) + y_upper;

        //;---------------------------------------------
        // See if muon passes through detector

        if( abs(y1 - y) < L/2 && abs (x1 - x) < W / 2)
        {
            //Muon passes through
            counts[pidx]++;
            sum_p[pidx] += p;
            sum_p_squared[pidx] += (p*p);
        }

        else
        {
            //Muon doesn't pass through
        }
    }

    } // END double for loops over paddles & trials

    double th_avg;
    double th_var;

    // Now output the results. Go to a file if specified, otherwise use stdout.
    if (fout == NULL) {
        fout = &cout;
    }
    for (int pidx = 0; pidx < n_lpanels; pidx++) {
        double th_avg = sum_p[pidx] / counts[pidx];
        double th_var = sum_p_squared[pidx] / counts[pidx] - th_avg*th_avg;
        *fout << "Paddle No.,Distance,Angle,Average Angle," <<
            "Angle Stdev,Counts" << endl;
        *fout << pidx << "," <<
            lower_x[pidx] << "," <<
            atan(x1/z0)*deg_per_rad << "," <<
            th_avg*deg_per_rad << "," <<
            sqrt(th_var)*deg_per_rad << "," <<
            counts[pidx] << endl;
    }

    free(lower_x);
    free(counts);
    free(sum_p);
    free(sum_p_squared);

    /*
    fout << "Counts = " << counts << endl
         << "Distance = " << x1 << endl
         << "Angle = " << 180*atan(x1/z0)/pi << endl
    << "Average angle = " << 180*th_avg / pi << endl
    << "Angle variance = " << th_var*180*180 / (pi*pi) << endl
    << "Angle standard deviation = " << sqrt(th_var)*180/pi << endl;
    */

}



