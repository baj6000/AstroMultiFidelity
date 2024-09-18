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

################################################################################
#                    E X P O R T E D     C L A S S E S:
################################################################################

#-------------------------------------------------------------------------------
class AbstractSCModel( object ) :

  def __init__( self ) :

    return


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

  SC = StochasticCollocation()

  #  Now, we test the node selection routine
  Nodes = np.zeros( (2, 7) )
  Nodes[0,0] = -3.0
  Nodes[0,1] = -2.0
  Nodes[0,2] = -1.0
  Nodes[0,3] =  0.0
  Nodes[0,4] =  1.0
  Nodes[0,5] =  2.0
  Nodes[0,6] =  3.0

  np.random.seed(1)
  Nodes = np.random.randn(2,7)
  #Nodes = np.random.randn(6,20)

  nodes, gram = SC._getOrderedNodes( Nodes, 10 )

  lkjlkj

  orderedNodes = []
  allNodes = list(range( Nodes.shape[1]))
  includedNodes = list(allNodes)
  print(includedNodes)
  while len(includedNodes) > 0 :
    if len(includedNodes) == 1 :
      orderedNodes.append( includedNodes[0] )
      break
    nodes, gram = SC._getOrderedNodes( Nodes[:,includedNodes], 1 )
    #print nodes
    for jj in range( len(nodes) ) :
      curNode = includedNodes[nodes[jj]]
      print(curNode)
      orderedNodes.append( curNode )
    includedNodes = list( set(allNodes) - set(orderedNodes) )
    #print includedNodes
    
  print(orderedNodes)
  for ii in range( len(orderedNodes) ) :
    print(np.sqrt( np.sum( Nodes[:,orderedNodes[ii]]**2. ) ))
  nodes = orderedNodes[0]
  #print gram 

  import matplotlib.pyplot as plt
  #print Nodes
  #print Nodes[:,nodes]

  plt.figure()
  plt.plot( Nodes[0,nodes], Nodes[1,nodes], 'o', ms=10, alpha=.5,  )
  plt.scatter( Nodes[0,:], Nodes[1,:] )
  #plt.show()

  return

if __name__ == "__main__":
  main()
