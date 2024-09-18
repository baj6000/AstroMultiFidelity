#!/usr/local/bin/python
'''

 <<Description>>


 <<Summary>>

'''

__author__ = 'Brandon A. Jones'
__version__ = '$Revision$'[11:-2]
__date__ = '$Date$'[7:26]

################################################################################
#                     I M P O R T     L I B R A R I E S
################################################################################

import numpy as np
import scipy.linalg as spl

import StochCollocation.stochCollocation as SC

################################################################################
#                    E X P O R T E D     C L A S S E S:
################################################################################

#-------------------------------------------------------------------------------
class BiFidelityModel( SC.AbstractSCModel ) :

  def __init__( self ) :

    self.haveSampleGen_    = False
    self.haveLowFidModel_  = False
    self.haveHighFidModel_ = False

    self.haveLowFidelitySamples_ = False
    self.haveRandomInputs_       = False

    self.numHFSamples_ = None
    self.lowRankTol_   = 1e300

    self.haveCVSamples_  = False
    self.haveErrorBound_ = False

    return


  def _getOrderedNodes( self, snapshot, numNodes ) : # KEEP
    '''
    Implements pivoted Cholesky decomposition Algorithm 1 of 

    Zhu, Narayan, and Xiu, ``Computational Aspects of Stochastic
    Collocation with Multifidelity Models", J. Uncert. Quant., Vol 2.,
    pp. 444-463, 2014.
    '''

    V = snapshot.copy()
    M = snapshot.shape[1]

    #  Initialize the ensemble for each parameter z_m;
    weights = np.zeros( (M,) )
    for ii in range( M ) :
      weights[ii] = np.dot( snapshot[:,ii].T, snapshot[:,ii] )

    #  Initialize the nodal selection vector and Cholesky factor;
    P = []
    L = np.zeros( (numNodes, M) )
    L = np.zeros( (M, numNodes) )

    indexMap = list(range( len(weights)))

    for n in range( numNodes ) :

      #  Find the next interpolation point (the next pivot);
      p = np.argmax( weights[n:] )
      p += n

      e = weights[p]

      #  Avoid ill-conditioning
      if e < 1e-10 :
        n = n-1
        break

      #  Update P and exchange column n and p in V;
      P.append( indexMap[p] )
      V[:,[n,p]] = V[:,[p,n]]
      weights[[n,p]] = weights[[p,n]]
      L[[n,p],:] = L[[p,n],:]

      #  Since we are flipping around columns, we need to make sure that we
      #  include a correct map of elements of weights to samples in snapshot
      tempInd = indexMap[p]
      indexMap[p] = indexMap[n]
      indexMap[n] = tempInd

      #  Update L;
      L[n,n] = np.sqrt( weights[n] )
      for t in range( n+1, M ) :
        r = np.dot(  V[:,t].T, V[:,n] ) \
                   - np.sum( L[t,0:numNodes-1]*L[n,0:numNodes-1]  )
        weights[t] = weights[t] - r*r/weights[n]
        L[t,n] = r / L[n,n]

    L = L[0:n+1,0:n+1]
    gramian = np.dot( L, L.T )

    return P, L
    #return P, gramian

  def setRandomSampleGen( self, fcn ) :
    '''
    Set the function used to generate random inputs to the solvers
    '''

    self.sampleGenerator_ = fcn
    self.haveSampleGen_   = True

    return

  def setLowFidelityModel( self, fcn, context ) :
    '''
    Set the function that takes random inputs and generates low-fidelity 
    samples
    '''

    self.lowFidelityModel_   = fcn
    self.lowFidelityContext_ = context
    self.haveLowFidModel_    = True

    return

  def getNumHighFidelitySamples( self ) :
    '''
    Get the number of high-fidelity samples to use.  This will override
    selecting the number based on a convergence tolerance
    '''
    return self.numHFSamples_

  def setNumHighFidelitySamples( self, numSamples ) :
    '''
    Set the number of high-fidelity samples to use
    '''
    self.numHFSamples_ = numSamples
    return

  def setLowFidelityConvergenceTol( self, tol ) :
    '''
    Set the convergence tolerance for the low-rank approximation
    '''
    self.lowRankTol_ = tol
    return

  def setHighFidelityModel( self, fcn, context ) :
    '''
    Set the function that takes random inputs and generates high-fidelity
    samples
    '''

    self.highFidelityModel_   = fcn
    self.highFidelityContext_ = context
    self.haveHighFidModel_    = True

    return

  def getHighFidelitySamples( self, newInputs=False ) :
    '''
    If available, get the high-fidelity samples
    '''
    if self.haveHighFidelitySamples_  and newInputs is False :
      projVector = np.dot( self.lfBasis_, \
                           self.lowFidelitySamples_ )
      coeffs = spl.cho_solve( self.cho_solveInputs_, projVector )
      return np.dot( self.highFidelitySamples_, coeffs )

    elif newInputs is not False :

      lfSamples = self.lowFidelityModel_( newInputs, \
                                          self.lowFidelityContext_)
      projVector = np.dot( self.lfBasis_, lfSamples )
      coeffs = spl.cho_solve( self.cho_solveInputs_, projVector )
      return np.dot( self.highFidelitySamples_, coeffs )

    else :
      return None

  def generateApproximation( self, nSamples ) :
    '''
    Generate the bi-fidelity approximation
    '''

    #  Make sure we have all of the inputs required to run the generation
    if   ( not self.haveSampleGen_ and not self.haveRandomInputs_ ) \
      or not self.haveLowFidModel_ \
      or not self.haveHighFidModel_  :
      raise TypeError( "class not configured before generating approximation" )

    #  If we don't have the random inputs, then generate them
    if not self.haveRandomInputs_ :
      self.randInputs_ = self.sampleGenerator_( nSamples )
      self.haveRandomInputs_ = True

    #  Generate the low-fidelity samples to be used later
    self.lowFidelitySamples_ = self.lowFidelityModel_( self.randInputs_, \
                                                       self.lowFidelityContext_)
    self.haveLowFidelitySamples_ = True

    #  Get the dimension of each sample and store for future reference
    uDim = self.lowFidelitySamples_.shape[0]


    #  Check to see if we should define limits on the number of high-fidelity 
    #  samples to consider based on user inputs
    if self.numHFSamples_ is not None :
      minRank = self.numHFSamples_
      maxRank = self.numHFSamples_

    else :
      minRank = 1
      maxRank = uDim

    #  Now, let's generate the list of important samples and prepare for future
    #  solutions via the collocation method
    for testRank in range( minRank, maxRank+1 ) :

      #  Get the list of important samples and the lower-trangular Cholesky
      #  decomposition of the grammian matrix
      self.importantSamples_, self.grammian_ = self._getOrderedNodes( \
                                                     self.lowFidelitySamples_, \
                                                     testRank )

      #  Setup the Cholesky solver inputs for future calculation of corrected
      #  samples
      self.cho_solveInputs_ = ( self.grammian_, True )
      self.lfBasis_ = self.lowFidelitySamples_[:,self.importantSamples_].T

      projVector = np.dot( self.lfBasis_, \
                           self.lowFidelitySamples_ )
      coeffs = spl.cho_solve( self.cho_solveInputs_, projVector )

      self.approxLFSol_ = np.dot( self.lfBasis_.T, coeffs )
      errMat =  self.approxLFSol_ - self.lowFidelitySamples_

      self.infError_ = np.max( np.abs( errMat ) )
      self.l2Error_  = np.linalg.norm( errMat, 2 )

      #if self.l2Error_ < self.lowRankTol_ :
      if self.infError_ < self.lowRankTol_ :
        break

    #  Store the low-fidelity surrogate coefficients in case we want to compute
    #  error bounds later
    self.coeffs_ = coeffs.copy()

    #  We have the rank of the solution, now let's save it for future use
    self.numHFSamples_ = testRank

    #  Now, we need to generate the important high-fidelity samples.
    self.highFidelitySamples_ = self.highFidelityModel_( \
                                  self.randInputs_[:,self.importantSamples_], \
                                  self.highFidelityContext_ )
    self.haveHighFidelitySamples_ = True

    return
