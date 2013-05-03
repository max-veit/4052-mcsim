#!/usr/bin/env python

import numpy as np
from numpy import pi,cos,arccos
import random

def mk_cosk_dist(k, th_max):
    """Return a function that transforms two uniformly-distributed
    random variates into two angles, the zenithal one distributed
    as cos^k(theta). More specifically, returns a function that applies
    the inverses of the CDFs of the marginal distributions of zenithal
    and azimuthal angle, where the zenithal angle is distributed as
    cos^k(theta) and azimuthal is uniformly distributed.

    The generated function returns a tuple (theta, phi) that gives the
    zenithal and azimuthal angles, respectively.

    The function uses only operations that are also NumPy ufuncs, and is
    thus suitable to be called with NumPy vector/array arguments.
    """

    # Normalization constant for the theta distribution
    th_norm = 1 - (cos(th_max))**(k+1)

    # The distribution function (inverse CDF) itself
    def cosk_dist(u_th, u_ph):
        theta = arccos((1 - th_norm*u_th)**(1 / (k+1)))
        phi = u_ph / (2*pi)
        return (theta, phi)

# The default distribution to use
cos2_dist = mk_cosk_dist(2, pi/2)

def get_counts(panel_pos, n=1E7, panel_ht=123.5, panel_wd=19, panel_len=70,
        dist=cos2_dist, blksize=1024):
    """Generate random particle events with a given angular distribution,
    originating at the top detector panel, and count how many intercept the
    bottom panel(s).

    Arguments:
    panel_pos       The horizontal positions (in the x-direction) of the
                    bottom panels relative to the top panel
    n               Number of random trials to generate
    panel_ht        The height of the top panel above the bottom ones
    panel_wd        The width of each detector panel
    panel_len       The length of each detector panel
    dist_trans      Angular distribution to use to generate events.
                    Must be specified in the form of a function that takes
                    two uniformly-distributed random variates and returns
                    two random variates (theta and phi) jointly distributed
                    as the desired angular distribution. The first returned
                    value should correspond to angle off the zenith, the second
                    to azimuthal angle measured from the positive x-direction.
    blksize         The sizes of trial vectors to use to make this function
                    more efficient with vector processing.
    """
    # Do magic here.
