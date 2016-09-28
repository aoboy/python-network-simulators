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
rsObj = np.random.mtrand.RandomState(simSeed)

class Node:
    def __init__(self, id, avgNodeDegree):
      self.id = id  
      self.neighLen = avgNodeDegree
      self.chId = 0
    
      self.ownSlot = None
      self.ownChannel =None
     
      self.state = False
      self.neighVec   = []
      self.notKnownNeigh=[]

      self.neighVec.append(self.id)  
    
    def initNotKnown(self):
        self.notKnownNeigh = []
        for i in range(1, self.neighLen + 1):
            if i not in self.neighVec:
                self.notKnownNeigh.append(i)
    
    def flushNeighbors(self):
       self.neighVec = []    
       self.notKnownNeigh=[]
       
       self.neighVec.append(self.id) 
       
       #self.initNotKnown()
           
    def addNeighbor(self, neighId):
        if neighId not in self.neighVec:
            self.neighVec.append(neighId)
            ''''''
            if neighId in self.notKnownNeigh:
                self.notKnownNeigh.remove(neighId)
            
    def addNotKnown(self, notID):
        if notID not in self.notKnownNeigh:
            self.notKnownNeigh.append(notID)
                  
    def neighVecFull(self):
	 return len(self.neighVec) == self.neighLen
  
    def addChannel(self, chId):
        self.chId = chId
        
    def addOwn(self, channel, ownSlot):
        
        self.ownSlot = ownSlot
        self.ownChannel = channel
        
    def isOwnSlot(self, slotOffset):
        if self.ownSlot == slotOffset:
            return True
        return False        
        
    def transmitOwn(self, ptOwn):
        
        self.state = False        
        self.chId = self.ownChannel
                
        p = rsObj.uniform(0,1)
        
        if p < ptOwn:
            self.state = True
               
    def disseminate(self, node):
        
        self.addNeighbor(node.id)
        
        for idx, nid in enumerate(node.neighVec):
            #self.addNeighbor(nid)
            if nid != self.id:
                self.addNotKnown(nid)

                
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
        channel = i%channelLen
        ownSlot = 1 + i%(T_period-1)
        node.addOwn(channel, ownSlot)
        lnodes[i]= node
    return lnodes
#---------------------------------------------------------------------------
def flushIds(Lnodes):
    for i, node in enumerate(Lnodes.values()):
        Lnodes[node.id].flushNeighbors()
#---------------------------------------------------------------------------
'''returns True if every node have discovered all its neighbors '''
def all2allComplete(Lnodes):
    for node in Lnodes.values():
        if node.neighVecFull() == False:
            return False
    return True
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
def plot_cdf(a):
    sorted_v=np.sort( a )
    yvals=np.arange(len(sorted_v))/float(len(sorted_v))
    plt.plot(sorted_v, yvals)   
    plt.show()
    #plt.close()
    print 'gonga....'
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
def runSim():
  
  
      
  global csvFd, CSize, KSize  
  
  global maxHopCount, T_period
  
  T_period = 3
  
  superFrames=[2,3,4,5,6,10]
  ptOwn    = 0.2
  
  numRepeat   = 20
    
  #nodes = [2, 4]+range(5, 51,5)
  channels=[2,4,8,10,12,14,16]  
  
  nodes = [100]#range(60, 101,10)
  channels=[4]#[4,8]
  
  
  for frameLen in [2]: #superFrames:
      T_period = frameLen
      
      print '----------------------------------------------'
      
      for num_nodes in nodes: 
          for channel in channels:
              pt = funcProbOf(channel, num_nodes)
    
              vec_times = []
              ##t_slot = 1
              print '==> [T=',frameLen,' N=',num_nodes,' K=',channel,']'
              
              Lnodes = createNodes(num_nodes, channel, T_period)
              
              for t  in range(numRepeat):
                  #A=np.eye(n)
                  t_slot = 0
                  '''clear neighbors list'''
                  flushIds(Lnodes)
                  
                  while all2allComplete(Lnodes) == False:
                      t_slot = t_slot + 1 
                      
                      '''THIS IS BEACONING SLOT'''
                      if t_slot%T_period == 0:
                      
                          '''each node uniformly at random selected a channel'''
                          rand_ch = rsObj.randint(0, channel, num_nodes).tolist()
                          
                          '''Transmission vector'''
                          tx_vec = rsObj.uniform(0,1, num_nodes).tolist()
                          
                          for ii, node in enumerate(Lnodes.values()):
                              '''Check if node is receiver'''
                              if tx_vec[ii] > pt:
                                  txList = []
                                  for jj, nodeX in enumerate(Lnodes.values()):
                                      if node.id != nodeX.id and rand_ch[ii] == rand_ch[jj]:
                                          if tx_vec[jj] < pt:
                                              txList.append(nodeX)
                                  if len(txList) == 1:
                                      Lnodes[node.id].disseminate(txList[0])
                          '''END of TX'''   
                      else:
                            '''THIS IS An ANNOUNCEMENT SLOT-NON BEACON'''
                            for jj, node in enumerate(Lnodes.values()):
                                '''determine which slot we are'''
                                slotOffset= 1 + t_slot%(T_period-1)
                                
                                if node.isOwnSlot(slotOffset):
                                    '''it is my slot.. transmit with probability: pOwn'''
                                    #Lnodes[node.id].transmitOwn(ptOwn)
                                    Lnodes[node.id].transmitOwn(ptOwn)
                                    '''@note HERE MAYBE try to CHECK status here to determine if we...'''
                                else:
                                    '''first I select one channel of not heard nodes yet'''
                                    vCh=[nodeC.chId for nodeC in Lnodes.values() if nodeC.isOwnSlot(slotOffset) == True and nodeC.id != node.id and nodeC.id not in node.neighVec]
                                    #vCh=[nodeC.id for nodeC in Lnodes.values() if nodeC.isOwnSlot(slotOffset) and nodeC.id not in node.neighVec and nodeC.id != node.id]
                                    if len(vCh):
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
                                '''I am a receiver'''
                                if node.state == False:
                                    txList=[]
                                    '''check for transmitters'''
                                    for nodeTx in Lnodes.values():
                                        '''Not Myself, and Transmitter is on the same channel I am..'''
                                        if nodeTx.state == True and node.id != nodeTx.id and node.chId == nodeTx.chId:
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
