#include <iostream>
#include <cstdlib>
#include <time.h>
#include <stdio.h>
#include <cmath>

using namespace std;

// Define pi
const double pi = 3.14159;


int main(int argc, char** argv)
{
    char *usagemsg = "Usage: %s [-n num] [-l length] [-w width] [-h height]
//;-------------------------------------------
//; Physical Properties of System

    //;---------------------------------------
    //; Panel Dimensions

        double L, W;         //Length, Width

        L = 70;      //Panel length (cm) along north-south (y-axis)
        W = 19;      //Panel width (cm) along east-west (x-axis)

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
    long trial_max = 10000000; //10 million trials
    int counts = 0;


    double x_upper, y_upper;    // Random spot on upper panel
    double phi, p;              // Random orientation of muon
    double x, y;                // Where muon lands on ground

    double sum_p = 0;           //Radians
    long double sum_p_squared = 0;   //Radians


    while (trial_num < trial_max)
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
        counts++;
        sum_p = p + sum_p;
        sum_p_squared = (p*p) + sum_p_squared;
    }

    else
    {
        //Muon doesn't pass through
    }

    trial_num++;

    }

    double th_avg = sum_p/counts;
    double th_var = sum_p_squared/counts - th_avg*th_avg;

    cout << "Counts = " << counts << endl
         << "Distance = " << x1 << endl
         << "Angle = " << 180*atan(x1/z0)/pi << endl
    << "Average angle = " << 180*th_avg / pi << endl
    << "Angle variance = " << th_var*180*180 / (pi*pi) << endl
    << "Angle standard deviation = " << sqrt(th_var)*180/pi << endl;

}


