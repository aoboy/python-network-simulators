#!/usr/bin/env python
# -*- coding: utf-8 -*-
# make a horizontal bar chart
"""
Dept. of Automatic Control
    School of Electrical Engineering
        KTH Kungliga Tekniska Högskolan (Royal Institute of Technology)
Created on Sat Dec 01 10:53:52 PM 2014

@author: gonga

New ideas of accerelating epidemic discovery.
A frame consists of: (1) a discovery slot, where nodes use multichannel epidemic
discovery, and (2) listening(announcement) where each node picks a slot to 
listen for a specific neighbor it has not yet heard 
"""

import matplotlib as mpl
import matplotlib.pyplot  as plt
from pylab import *

import math
import numpy as np
from numpy import mean
from scipy import stats
import scipy as sp


from matplotlib import rc
plt.rcParams['ps.useafm'] = True
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#plt.rcParams['pdf.fonttype'] = 3 #[42 or 3]
prop = matplotlib.font_manager.FontProperties(size=8)


simSeed= 997
snapShotCounter=0
rsObj = np.random.mtrand.RandomState(simSeed)

class Node:
    def __init__(self, id, avgNodeDegree):
      self.id = id  
      self.neighLen = avgNodeDegree
      self.chId = 0
    
      self.ownSlot = None
      self.ownChannel =None
     
      self.state = False
      self.networkVec   = []
      self.notKnownNeigh=[]
      self.neighbors  = []
      self.neighbors1Hop = []
      
      self.xc0 = 0
      self.yc0 = 0
      self.hopCount    = 0
      self.radioRadius = 30
      self.Ptrans = 0.1
      '''Const related when a node starts discovery..'''
      self.active = False 
      '''Generate a random number between 0 and 200.'''
      self.startTime = 0
      
      self.Phase_i = 1
      self.adptTxP = 1./2.
      self.adaptNumNeighbors = 2**(self.Phase_i+1)
      self.phaseLenght =8*self.Phase_i*2**self.Phase_i
      
      self.adptTxP = 1./self.adaptNumNeighbors

      self.networkVec.append(self) 
      
    def updateParameters(self, simTime):
        
        if simTime != 0 and simTime%self.phaseLenght == 0:
            
            self.Phase_i = self.Phase_i + 1
            
            newWindow = 8*self.Phase_i*2**self.Phase_i
            
            self.phaseLenght = self.phaseLenght + newWindow
                        
            self.adaptNumNeighbors = 2**(self.Phase_i+0)
            
            self.adptTxP = 1./self.adaptNumNeighbors
        
    
    def addCoordinates(self, x, y):
       self.xc0 = x
       self.yc0 = y
       
    def isNeighbor(self, nodeXY):
       
       if nodeXY.id != self.id:
         d = math.sqrt( (nodeXY.xc0-self.xc0)**2 + (nodeXY.yc0 -self.yc0)**2 )
       
         if d < self.radioRadius:
               return True
               
       return False
       
    def isNeighborByID(self, nID):
        if self.id != nID and nID in self.neighbors:
            return True
        return False
    
    def flushNeighbors(self):
       ''''''
       self.networkVec    = []    
       self.notKnownNeigh = []
       self.neighbors1Hop = []
       
       self.networkVec.append(self) 

    def addReceived(self, node):
        if node not in self.networkVec:
            self.networkVec.append(node)  
        else:
            idx = self.networkVec.index(node) 
            hop = self.networkVec[idx].hopCount
            
            if node.hopCount < hop:
                self.networkVec[idx] = node

       
    def addNeighbor(self, neighId):
        if neighId not in self.neighbors1Hop:
            self.neighbors1Hop.append(neighId)
            
            '''We have discovered the node.. remove it from not Known nodes'''
            if neighId in self.notKnownNeigh:
                self.notKnownNeigh.remove(neighId)
            
    def addNotKnown(self, notID):
        if notID not in self.notKnownNeigh:
            self.notKnownNeigh.append(notID)
                  
    def allNetworkDiscovered(self):
        return len(self.networkVec) == self.neighLen
                  
    def allNeighborsDiscovered(self):
	 return len(self.neighbors) == len(self.neighbors1Hop)
  
    def addChannel(self, chId):
        self.chId = chId
        
    def addOwn(self, channel, ownSlot):
        
        self.ownSlot = ownSlot
        self.ownChannel = channel
        
    def isOwnSlot(self, slotOffset):
        if self.ownSlot == slotOffset:
            return True
        return False        

    def selectOwn(self, pOwn):
        retV= False
        p = rsObj.uniform(0,1)
        if p < pOwn:
            retV = True
        return retV
        
    def transmitOwn(self, ptOwn):
        
        self.state = False        
        self.chId = self.ownChannel
                
        p = rsObj.uniform(0,1)
        
        pt = 0.5
        if len(self.neighbors1Hop) != 0:
             pt = 1./(len(self.neighbors1Hop) + 1)  
        else:
            pt = 0.5
        
        if p < pt:#Own:
            self.state = True
            
    def transmitAdapt(self):
        self.state = False 
        self.chId = self.ownChannel
        
        p = rsObj.uniform(0,1)
        
        if p < self.adptTxP:
            self.state = True
               
    def disseminate(self, node):
        
        node.hopCount = 1
        
        if node.id != self.id:
            self.addNeighbor(node)
        
            '''Add overheard nodes to list'''
            for nodej in node.networkVec:            
                
                if nodej.id != self.id:
                    nodej.hopCount = nodej.hopCount + 1
                    hop = nodej.hopCount                
    
                    if hop == 2:
                        self.addNotKnown(nodej)
                    
                    '''add to the network table..'''
                    #print 'GGG', self.id, nodej.id
                    self.addReceived(nodej)

                
    def transmit(self, pt):
        self.state = False
        p = rsObj.uniform(0,1)
        if p < pt:
            self.state = True  
    
    def setReceiver(self, channel):
        self.state = False
        self.chId = channel
    
    def receive(self, pl):
        p = rsObj.uniform(0,1) 
        
        if p > pl:

            return True
        else:
            return False
            
    def funcProbOf(C, N):
        a = N
        b = 2*C+N-1
        c = C
        f2 = (b - np.sqrt(math.pow(b,2)-4*a*c))/(1.0*2*a)
        return f2    

def confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), sp.stats.stderr(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m-h, m, m+h
#---------------------------------------------------------------------------
def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), sp.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m-h, m, m+h
#---------------------------------------------------------------------------
def funcProbOf(C, N):
    a = N
    b = 2*C+N-1
    c = C
    f2 = (b - np.sqrt(math.pow(b,2)-4*a*c))/(1.0*2*a)
    return f2
#--------------------------------------------------------------------------
def loss_p(p, n):
    a = 1.0
    for i in range(n):
        a = a*(1-p)

    return a
#---------------------------------------------------------------------------
'''Create all nodes in the Network'''
def createNodes(n, channelLen, T_period):
    lnodes={}
    for i in range(1, n+1):
        node = Node(i,n)
        #channel = rsObj.randint(channelLen) #i%channelLen
        channel = i%channelLen
        ownSlot = 1 + i%(T_period-1)
        node.addOwn(channel, ownSlot)
                
        lnodes[i]= node
    return lnodes
#---------------------------------------------------------------------------
def flushIds(Lnodes):
    for i, node in enumerate(Lnodes.values()):
       Lnodes[node.id].flushNeighbors()
       #Lnodes[node.id].networkVec    = []    
       #Lnodes[node.id].notKnownNeigh = []
       #Lnodes[node.id].neighbors1Hop = []  
       
       #Lnodes[node.id].networkVec.append(node)
#---------------------------------------------------------------------------
def printNeighbors(Lnodes):
    for node in Lnodes.values():
        print 'NodeID:', node.id, 'Neighbors: ',
        #for n in node.neighbors1Hop:
        if 1:
            #print n.id,
            print len(node.neighbors1Hop), len(node.neighbors), len(node.networkVec), node.Phase_i
           
    #print '\n'
#---------------------------------------------------------------------------
'''returns True if every node have discovered all its neighbors '''
def all2allComplete(Lnodes):
    for node in Lnodes.values():
        if node.allNetworkDiscovered() == False:
            return False
    return True
#---------------------------------------------------------------------------
def allNeighborsDiscovered(Lnodes):
    vec90pComplete=[]
    #for node in Lnodes.values():
    #    if len(node.neighbors1Hop) != len(node.neighbors):
    #        return False
    #return True
    for node in Lnodes.values():
        if len(node.neighbors1Hop) >= 0.8*len(node.neighbors):
            vec90pComplete.append(1)
    
    #print 'Percentage Nodes:', 1.*len(vec90pComplete)/len(Lnodes)       
    if len(vec90pComplete) >= 0.8*len(Lnodes):
        return True
    return False
#---------------------------------------------------------------------------    
def addTopology(Xmax, Ymax, Nodes):
    lnodes={}
    xco = np.arange(0, Xmax+0.25, 0.25)
    yco = np.arange(0, Ymax+0.25, 0.25)    
    '''add coordinates'''
    for node in Nodes.values():
        x = np.random.choice(xco)
        y = np.random.choice(yco)
        
        Nodes[node.id].addCoordinates(x,y)
        
        lnodes[node.id] = Nodes[node.id]
        
    '''add neighbors based on X and Y'''
    for nodei in lnodes.values():
        for nodej in lnodes.values():
            if nodei.isNeighbor(nodej):
                lnodes[nodei.id].neighbors.append(nodej)
    '''Print Coordinates'''           
    #for node in lnodes.values():
    #    print [node.id, node.xc0, node.yc0],'->',node.neighbors
    
    #saveTopology('./NodePositions2.dat', lnodes)
            
    return lnodes
    
def plot_cdf(a):
    a.sort()
    yvals=np.arange(len(a))/float(len(a))
    plt.plot(a, yvals )  
    
    grid(True)
    #plt.close()
    ylabel('CDF')
    xlabel('Worst case Latency (slots)')
    plt.show()
#---------------------------------------------------------------------------   
#---------------------------------------------------------------------------
def write_file(fname, a):
    with open(fname,'w') as fw:
        for i, val in enumerate(a):
            if i+1 < len(a):
                fw.write(str(val)+',')
            else:
                fw.write(str(val))
#---------------------------------------------------------------------------
def read_file(fname):
    with open(fname,'r') as fr:    
        dt = fr.read().split(',')
        xx =[int(x) for x in dt]
        return xx    
#--------------------------------------------------------------------------- 
def saveTopology(fname, listNodes):
    with open(fname, 'w') as fw:
        for ii, node in enumerate(listNodes.values()):
            if ii+1 < len(listNodes):
                fw.write(str(node.xc0)+',')
                fw.write(str(node.yc0)+':')
            else:
                fw.write(str(node.xc0)+',')
                fw.write(str(node.yc0))                                          
#---------------------------------------------------------------------------         
def readTopology(fname):
    with open(fname,'r') as fr:
        dt = fr.read().split(':')
        #print dt
        xx = [(float(x.split(',')[0]),float(x.split(',')[1])) for x in dt]

        return xx
#---------------------------------------------------------------------------                  
def nodesAddCoordinates(network, lnodes):
    Lnodes={}
    if len(network) < len(lnodes):
        print 'WARNING:::: COORDINATES VECTOR SMALLER THAN SIZE OF THE NETWORK'
        exit(0)
    print '''==>downloading network topology.....'''
    
    
    #plt.axes()
    
    for node, coord in zip(lnodes.values(), network):
            x, y = coord
            lnodes[node.id].addCoordinates(x,y) 
            Lnodes[node.id] = lnodes[node.id]
            
    '''add neighbors based on X and Y'''
    print '==>computing neighborhood...'
    
    for nodei in Lnodes.values():
        
        for nodej in lnodes.values():
            if nodei.id != nodej.id and nodei.isNeighbor(nodej):
                #nodei.neighbors.append(nodej)
                Lnodes[nodei.id].neighbors.append(nodej)
        
        '''PLOT NODE LOCATION'''
        x, y = nodei.xc0, nodei.yc0
        
        #circxy = plt.Circle((x, y), radius = 2, fc='r')#,

        #plt.gca().add_patch(circxy)
        
        #plt.text(x-0.5*2, y-0.5*2, str(nodei.id), fontsize=10, color='k')
        
        '''PLOT LINKS'''
        for neighNode in nodei.neighbors:
            x = (nodei.xc0, neighNode.xc0)
            y = (nodei.yc0, neighNode.yc0)
            #plot(x,y, color='b', linestyle='-', linewidth=0.2)
    

    
    #xlim([0,100])
    #ylim([0,100])
    

    
    #plt.grid(True)
    #title('[100m x 100m]::::::RANDOM NETWORK TOPOLOGY::::::')
    #xlabel('X-$axis$ (m)')
    #ylabel('Y-$axis$ (m)')            
    #plt.show()

    return Lnodes

#--------------------------------------------------------------------------- 
def networkUpdate(lnodes):
    plt.clf()
    
    plt.axes()        
        
    '''PLOT NODE LOCATION'''
    for nodei in lnodes.values():
        x, y = nodei.xc0, nodei.yc0
        #mpl.patches.Circle(x,y, radius=3, color='red')
        circxy = plt.Circle((x, y), radius = 2, fc='k')#,
        plt.gca().add_patch(circxy)
        plt.text(x-0.5*2, y-0.5*2, str(nodei.id), fontsize=10, color='white')
        
        '''PLOT LINKS'''
        for neighNode in nodei.neighbors1Hop:
            
            x = (nodei.xc0, neighNode.xc0)
            y = (nodei.yc0, neighNode.yc0)
            plot(x,y, color='b', linestyle='-', linewidth=0.2)
    
    '''RESTORE GRID AND AXIS NAMES'''
    xlim([0,100])
    ylim([0,100])
    
    plt.grid(True)
    title('[100m x 100m]:::::: [Snapshot] NETWORK TOPOLOGY::::::')
    xlabel('X-$axis$ (m)')
    ylabel('Y-$axis$ (m)')            
    plt.show()            
#---------------------------------------------------------------------------            
def runSim():
  
  global csvFd, CSize, KSize  
  
  global maxHopCount, T_period
  
  T_period = 3
  snapShotCounter = 0
  
  superFrames=[2,3,4,5,6,10]
  ptOwn    = 0.50
  
  numRepeat   = 10
  #nodes = [2, 4]+range(5, 51,5)
  channels=[2,4,8,10,12,14,16]  
  
  nodes = [100]#range(60, 101,10)
  channels=[8]#[4,8]
  
  X_MAX, Y_Max, NumNodes = (100, 100, 100)
  
  
  #networkCoord=readTopology('./NodePositions.dat')
  
  #print networkCoord
  
  #exit()
  K2=[129, 130, 147, 149, 152, 152, 155, 156, 163, 167]
  K4=[166, 190, 198, 199, 202, 206, 207, 207, 214, 217]
  K8=[494, 510, 517, 522, 569, 584, 586, 587, 607, 649]
  
  phase_count = 1
  
  for frameLen in [3]: #superFrames:
      T_period = frameLen
      
      print '----------------------------------------------'
      
      for num_nodes in nodes: 
          for channel in channels:
              #pt = funcProbOf(channel, num_nodes)
              pt = funcProbOf(channel, 10)
    
              vec_times = []
              ##t_slot = 1
              print '==> [T=',frameLen,' N=',num_nodes,' K=',channel,']'
              
              #Lnodes = createNodes(num_nodes, channel, T_period)
              
              Lnodes = {}
              
              '''Add topology'''
              #addTopology(100, 100, Lnodes)
              #return
              
              '''Add Coordinates to nodes..'''
              #Lnodes = nodesAddCoordinates(networkCoord, Lnodes)
              
              
              #return
              
              for t  in range(numRepeat):
                  #A=np.eye(n)
                  t_slot = 0
                  
                  '''clear neighbors list'''
                  flushIds(Lnodes)
                  
                  Lnodes = createNodes(num_nodes, channel, T_period)                  
                  
                  Lnodes = addTopology(X_MAX, Y_Max, Lnodes)
                  
                  #printNeighbors(Lnodes)
                  
                  #while all2allComplete(Lnodes) == False:
                  while allNeighborsDiscovered(Lnodes) == False:
                        t_slot = t_slot + 1 
                                   
                        if t_slot%20 == 0:
                           print '[round: ', t , t_slot,']'  
                           #networkUpdate(Lnodes)
                           
                           if t_slot%100 == 0:
                               ''''''
                               #printNeighbors(Lnodes)
                          
                        '''THIS IS BEACONING SLOT'''                     
                        #'''THIS IS A DIFFERENT SLOT-NON BEACON'''
                        for jj, node in enumerate(Lnodes.values()):
                            
                            '''Make sure we update the new Tx Probability'''
                            Lnodes[node.id].updateParameters(t_slot)
                            
                            #randChannel = rsObj.choice(range(channel))
                            #Lnodes[node.id].setReceiver(randChannel)
                            #if node.selectOwn(1.0):

                            if node.selectOwn(ptOwn):
                                '''it is my slot.. transmit with probability: pOwn'''
                                #Lnodes[node.id].transmitOwn(pt)
                                Lnodes[node.id].transmitAdapt()
                                '''@note HERE MAYBE try to CHECK status here to determine if we...'''
                            else:
                                '''first I select one channel of not heard nodes yet'''
                                vCh=[chRand for chRand in range(channel) if chRand != node.ownChannel]                                    
                                if  len(vCh):
                                    '''choose channel UNIFORMLY based on the not YET heard nodes'''
                                    randChannel=rsObj.choice(vCh)
                                    #randChannel = Lnodes[rsObj.choice(vCh)].chId
                                    Lnodes[node.id].setReceiver(randChannel)
                                else:
                                    '''OTHERWISE choose any random channel COLUMN'''
                                    randChannel = rsObj.choice(range(channel))
                                    Lnodes[node.id].setReceiver(randChannel)
                        '''Now select who succeeded'''
                        for kk, node in enumerate(Lnodes.values()):
                            '''....I am a receiver....'''
                            if node.state == False:
                                txList=[]
                                '''check for transmitters'''
                                for nodeTx in Lnodes.values():
                                    '''Not Myself, and Transmitter is on the same channel I am..'''
                                    if nodeTx.state and node.isNeighbor(nodeTx) and node.chId == nodeTx.chId:
                                        '''This is a transmitter'''
                                        txList.append(nodeTx)
                                if len(txList) == 1:
                                    '''ONE NODE TRANSMITTED..SUCCESS'''                                        
                                    Lnodes[node.id].disseminate(txList[0]) 
                  #end While
                  vec_times.append(t_slot)
              #end for
              '''PLOT CDF...'''
              plot_cdf(vec_times)   
              
              vec_times.sort()
              print vec_times
              
              dfname='./data/AdptV1'+'K'+str(channel)+'.csv'
              
              #write_file(dfname, vec_times)              
            
              discTime = mean(vec_times)#/(1-0.2)
              
              print 'Mean for',' [T=', T_period,'N =', num_nodes, 'K=',channel, round(pt, 3), '->', round(discTime, 2),']'                                    
      #end for
  ###
  '''comment'''
  ## 
#---------------------------------------------------------------------------
'''======================================================================'''
if __name__ == '__main__':
    
    runSim()      
    
 
'''======================================================================'''
