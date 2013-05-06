#!/usr/bin/env python

import numpy as np
from numpy import pi,cos,arccos,sin,tan,sum,size
from numpy.random import random_sample
from math import ceil

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
    # Not applicable for k=-1
    assert k != -1

    # Normalization constant for the theta distribution
    th_norm = 1 - (cos(th_max))**(k+1)

    # The distribution function (inverse CDF) itself
    def cosk_dist(u_th, u_ph):
        theta = arccos((1 - th_norm*u_th)**(1 / (k+1)))
        phi = u_ph / (2*pi)
        return (theta, phi)
    return cosk_dist

# The default distribution to use
cos2_dist = mk_cosk_dist(2, pi/2)

def cos2_dist(u_th, u_ph):
    theta = arccos((1 - u_th)**(1/3))
    phi = u_ph / (2*pi)
    return (theta, phi)

def get_counts(panel_pos, n=1E7, panel_ht=123.5, panel_wd=19, panel_len=70,
        dist=cos2_dist, blksize=1024):
    """Generate random particle events with a given angular distribution,
    originating at the top detector panel, and count how many intercept the
    bottom panel(s).

    A note on the coordinate system: The z-axis is the vertical direction,
    the x-axis is the axis along which the bottom detectors are aligned, and
    the y-axis is the axis along the panels' lengths. The z=0 plane is the
    ground plane (on which the bottom detectors rest), and x=0, y=0 corresponds
    to the lower left corner of the top detector (that is, the top detector
    is contained within the (+x, +y) quadrant).

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
    nblks = ceil(n / blksize)
    # Account for case when n is not a multiple of blksize
    if nblks * blksize > n:
        endblksize = n - (nblks-1) * blksize
    else :
        endblksize = blksize
    print("Using {} blocks of size {} (end: {})".format(nblks, blksize, endblksize))

    hits = np.zeros_like(panel_pos)
    th_acc = np.zeros_like(panel_pos)
    th2_acc = np.zeros_like(panel_pos)

    for blk in range(nblks):
        if (blk == (nblks - 1)):
            blksize = endblksize

        # Pick random angles and make sure they are appropriately distributed
        theta = random_sample(blksize)
        phi = random_sample(blksize)
        (theta, phi) = dist(theta, phi)
        # Also pick a random position on the top detector
        xt = random_sample(blksize) * panel_wd
        yt = random_sample(blksize) * panel_len
        # And determine where the intersections of the rays with the
        # groundplane are.
        xb = xt + cos(phi) * panel_ht*tan(theta)
        yb = yt + sin(phi) * panel_ht*tan(theta)

        # Count the number of intersections
        it = np.nditer([panel_pos, hits, th_acc, th2_acc],
                op_flags=[['readonly'],
                          ['readwrite'],
                          ['readwrite'],
                          ['readwrite']],
                )
        for pos, phits, th, th2 in it:
            intersect = ((xb >= pos) & (xb <= pos + panel_wd) & 
                    (yb >= 0) & (yb <= panel_len))
            if (blk == 0):
                print(sum(intersect))
            phits[...] += sum(intersect)
            # Keep track of the average angle and stdev
            th[...] += sum(theta[intersect])
            th2[...] += sum(theta[intersect]**2)

    return {'counts':hits, 'theta_avg':th_acc/n,
            'theta_stdev':sqrt(th2_acc/n - (th_acc/n)**2)}


