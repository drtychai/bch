#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#######################################################################
#
#       Author: Justin Angra
#       Last Modified: April 30, 2017
#
#       Overview:  
#           The code below calculates the solution of 
#           Z = log(e^X e^Y) for non-commutative X,Y in the 
#           P. Hall basis of Lie algebra. We have optimized the 
#           explicit Baker-Campbell-Hausdorff formula by
#           structuring it's computational alorgithm as a rooted
#           tree. With this structure, we hope to achieve a 
#           reduction in CPU time by 6 significant digits (with
#           units of seconds). 
#
#           The algorithm was not created, but only implemented by 
#           me. All credit to it's creation goes to Fernando Casas 
#           of Departament de Matematiques, Universitat Jaume I, 
#           Castell´on, Spain and Ander Murua of Konputazio 
#           Zientziak eta A.A. saila, Informatika Fakultatea, 
#           EHU/UPV, Spain
#
#           Please see https://arxiv.org/pdf/0810.2656.pdf
#
#
#       To do (in the follow order):
#           ** Fix pair generator to work past order 5
#           ** Create function to calculate coefficients
#           * Modify algorithm for general basis
#           * Create function to determine start of next order
#           * Add parser for client to pick order of calculation
#
#
#       Usage:
#           python bch.py -x X_MATRIX -y Y_MATRIX
#
#       Notes:
#           * Input matrices can be of any size
#           * Input matrices must be formatted as such:
#               1 2 1
#               3 5 6
#               -4 2 6
#             i.e. each element is separated by a space and each
#               row has it's own line (please see example matrices)
#           * Running against mathematica in simple cases, we see that 
#               by only going up to the 5th order causes a significant 
#               error. Most work needs to be done for this to be an
#               applicable program
#######################################################################

import argparse
import sys
import numpy as np
from time import clock  # timeing
import bch_coeff

def init_dict():
    """Returns two dictionary. Only coded up to 9th order.
    
    USAGE:
        init_dict()

    OUTPUT:
        order_dict   - dictionary of each iteration and if that iteration
                        is the first of the next order (True if so)
        len_dict     - dictionary of each iteration and it's 'length', i.e.
                        of which respective order it is
    """
    first_of_order = [2,3,4,6,9,15,24,42,72,128]
    order_dict = {}
    len_dict = {1:1}
    j=0
    for i in xrange(2,180):
        if i in first_of_order:
            order_dict.update({i:True})
            j+=1
        else:
            order_dict.update({i:False})
        len_dict.update({i:j})

    return order_dict, len_dict

def myComm(op1, op2):
    """Returns the commutator relation for given input.
    For commuting values function returns 0.

    USAGE:
        myComm(op1, op2)

    INPUT:
        op1    - matrix of the same dimentions as op2

        op2    - matrix of the same dimentions as op1
    OUTPUT:
             - matrix representing the commutator relation between op1 and op2
    """
    return np.dot(op1, op2)-np.dot(op2, op1)
    
def main():
    is_first_in_order, get_len = init_dict() # initialize dictionaries
    coeffs = bch_coeff.get_coeff()           # load coefficients
    order=5                                  # set order
    pairs = [] 
    
    # parser to get arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-x', '--arg1', dest='X_mat', required=True)
    parser.add_argument('-y', '--arg2', dest='Y_mat', required=True)
    globals().update(vars(parser.parse_args()))

    try:
        X = np.loadtxt(X_mat, dtype=np.float128)  # Read X matrix file
        Y = np.loadtxt(Y_mat, dtype=np.float128)  # Read Y matrix file
    except Exception as e:
        print e
        sys.exit(2)
    
    # Check if both matrices are of the same size
    assert(len(X) == len(Y)), "Matrices do not have same dimensions" 

    # Get matrix size
    mat_size = len(X)

    # generate ordered pairs 
    # *(current implementation is not accurate above 5th order)*
    i=3
    for n in xrange (2,order+1):
        for j in xrange(1,i-1):
            for k in xrange(j+1,i):
                if (get_len[j]+get_len[k] == n):
                    if (j==1 and is_first_in_order[k]):
                        pairs.append([k,j])
                        i+=1
                    elif j!=1:
                        pairs.append([k,j])
                        i+=1  
               
    #set initial values
    soln = [] #each elements will be a calculated commutator for a pair
    soln.append(X)
    soln.append(Y)
    
    #calculate commutator for each pair
    for i in xrange(len(pairs)):
        A = pairs[i][0]
        B = pairs[i][1]
        soln.append(coeffs[i]*(myComm(soln[A-1],soln[B-1])))
        
    #sum each single commutator solution
    Z = np.zeros((mat_size, mat_size), dtype=np.float)
    for i in xrange(len(soln)):
        Z+=soln[i]
    
    print "Z = " , Z 
    return 
    
if __name__ == "__main__":
    start = clock()
    main()
    print "Total run time: {0} seconds".format(clock() - start)