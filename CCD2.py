import numpy as np
# D is the boundary of the target test area (suppose D is the unit hypercube [0, 1]^q)
# q is the test area dimension
# x is the input point, and the vector of length q represents the coordinates of that point in the test area D
# points: A set of n points, each of which is a vector of length q
# Npoints: The region D is usually discretized by Npoints (N is much larger than n), each of which is a vector of length q
def get_subregion(x, D, q):
 # Calculates the index of the subregion Dt(x) in which the point is located
 index = 0
 for i in range(q):
        if x[i] >= 0.5:
            index += 2**i  # If x_i > 0.5, set the bit of the corresponding dimension to 1
 return index

def ND(points, D, q):
    ND = np.zeros(2**q, dtype=int)
    # For each point, calculate the subarea it belongs to and count it
    for point in points:
        subregion_index = get_subregion(point, D, q)
        ND[subregion_index] += 1
        
    return ND


def CCD(N, n, q, x, D, points, Npoints):
    sum_total = 0
    for point in Npoints:
        for t in range(2**q):
            N1 = ND(points, D, q)[t]
            N2 = ND(Npoints, D, q)[t]

            diff = abs(N1 / n - N2 / N)
            sum_total += diff**2
    CCD = (sum_total / (N * 2**q))**0.5
    return CCD