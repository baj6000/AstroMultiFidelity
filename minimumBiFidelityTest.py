#!/usr/local/bin/python
'''

 Notional use of the biFidelityModel class


 Sample software that includes a notional use of the bi-fidelity software.

'''

__author__ = 'Brandon A. Jones'
__version__ = '$Revision$'[11:-2]
__date__ = '$Date$'[7:26]

################################################################################
#                     I M P O R T     L I B R A R I E S
################################################################################

import numpy as np
import matplotlib.pyplot as plt

import StochCollocation.biFidelityModel as BFM

################################################################################
#                    E X P O R T E D     C L A S S E S:
################################################################################

#-------------------------------------------------------------------------------


################################################################################
#                  E X P O R T E D     F U N C T I O N S:
################################################################################

#-------------------------------------------------------------------------------


################################################################################
#             U N I T     T E S T     C A S E     F U N C T I O N:
################################################################################

#-------------------------------------------------------------------------------


################################################################################
#                       M A I N     F U N C T I O N:
################################################################################

def main():
  '''
  <Description>

  Parameters
  ----------
  inputs : description

  Returns
  -------
  outputs : description

  Examples
  --------
  Examples, if needed
  '''

  #  Here is a simple low-fidelity model.  Note that the dimensions of the
  #  output needs to be number of quantitites of interest x number of samples
  def lowFidFunction( randSamples, context ) :
    return np.cos( context['tspan']*randSamples )

  #  Here is a simple high-fidelity model.  Note that the dimensions of the
  #  output needs to be number of quantitites of interest x number of samples
  def highFidFunction( randSamples, context ) :
    return np.cos(context['tspan']*randSamples) \
               +0.1*np.sin(context['tspan']*randSamples)

  #  Here is our random input generator.  Note that the output needs to be the
  #  number of random inputs x the number of samples
  def genRandomSamples( nSamples ) :
    return np.random.rand( 1, nSamples )*0.1 + 0.1

  #  Set the seed for repeatability
  np.random.seed( 328 )

  #  Generate a time array.  We need to output at least more times than the
  #  rank of our approximation
  tspan = np.linspace( 100., 110., 31 ).reshape((-1,1))

  #  Set the "context" of our models.  These are just parameters that are
  #  passed through to the low-fidelity and high-fidelity models.
  lfContext = {}
  lfContext['tspan'] = tspan.copy()

  hfContext = {}
  hfContext['tspan'] = tspan.copy()

  #  Instantiate our class that will solve for the bi-fidelity model
  lrApprox = BFM.BiFidelityModel()

  #  Set the random input generator
  lrApprox.setRandomSampleGen( genRandomSamples )

  #  Set the low-fidelity model function.  Can instead be a callable class if
  #  desired
  lrApprox.setLowFidelityModel(   lowFidFunction, lfContext )

  #  Set the high-fidelity model function.  Can instead be a callable class if
  #  desired
  lrApprox.setHighFidelityModel( highFidFunction, hfContext )

  #  Set the tolerance for the low-fidelity model, which determines how many
  #  important samples to select.  This one will require some tuning for a
  #  specific application.
  lrApprox.setLowFidelityConvergenceTol( 1e-5 )

  #  Generate the bi-fidelity surrogate using 1000 low-fidelity samples
  lrApprox.generateApproximation( 1000 )

  #  Print out how many important samples were used, i.e., the rank of the
  #  approximation
  print('Generated %d rank solution'%(lrApprox.getNumHighFidelitySamples()))

  #  Now, generate 2000 independent samples.  We will use them to demonstrate
  #  the improved accuracy of our surrogate.
  randInputs = genRandomSamples( 2000 )

  #  Evaluate the high-fidelity model to get truth
  newOutputs = highFidFunction( randInputs, hfContext )

  #  Get the bi-fidelity determined solution
  newHFSamples = lrApprox.getHighFidelitySamples( randInputs )

  #  For comparison, get the low-fidelity samples
  newLFSamples = lowFidFunction( randInputs, lfContext )

  #  Now, we will plot the histogram of the low-fidelity and bi-fidelity
  #  errors.  We will see that the bi-fidelity approximation reduces the error
  #  in the surrogate.

  bfErrors = newHFSamples-newOutputs
  lfErrors = newHFSamples-newLFSamples

  plt.figure()
  plt.subplots_adjust(hspace=0.5)

  plt.subplot(2,1,1)
  plt.hist( lfErrors.reshape((-1,)), 30 )
  plt.xlabel( 'Error in Low-Fidelity Model', size=14 )

  plt.subplot(2,1,2)
  plt.hist( bfErrors.reshape((-1,)), 30 )
  plt.xlabel( 'Error in Bi-Fidelity Model', size=14 )

  plt.suptitle( 'Histogram of Errors', size=16 )

  plt.show()

  return

if __name__ == "__main__":
  main()
