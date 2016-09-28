  # -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 13:53:52 2014

@author: gonga
"""

#!/usr/bin/env python
# make a horizontal bar chart
import matplotlib
matplotlib.use('Agg')


import numpy as np
import math
import time
from time import gmtime, strftime, localtime
import Queue
import threading
from threading import Thread

from scipy import stats
import scipy as sp
import scipy.stats
#---------------------------------------------------------------------------


primeNums=[2063,2069,2081,2083,2087,2089,2099,2111,2113,2129]
simSeed = 997#43
rsObj = np.random.mtrand.RandomState(simSeed)

def mean_conf_int(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return (round(m-h,2), round(m,3), round(m+h,2))
#---------------------------------------------------------------------------
class Node:
    def __init__(self, id, avgNodeDegree):
      self.id = id  
      self.neighLen = avgNodeDegree
      self.chId = 0
      self.hc   = 0;
      self.state = False
      self.neighVec   = []
      self.adjMat     = []
      
      self.neighHc    = [0]
      
      self.adjMat.append(self.id)
      self.neighVec.append(self.id)      
    
    def flushNeighbors(self):
	self.neighVec = []
        self.neighVec.append(self.id)
        
    def txVecOnes(self):  
      for v in self.txVec:
	  if v == 0:
	    return False
      return True
      
    def updateTxVecIdx(self, idx):
      self.txVec[idx] = 1   
    
    def fillAdjMat(self, nodes):
        for nid in nodes:
            if nid not in self.adjMat:
                self.adjMat.append(nid)
    
    def addNeighbor(self, neighId):
        if neighId not in self.neighVec:
            self.neighVec.append(neighId)
            #if self.id == 2: print neighId     
    def addNeighborHc(self, neighId, hopc):
        if neighId not in self.neighVec:
            self.neighVec.append(neighId)
            self.neighHc.append(hopc)
            
        if neighId in self.neighVec:
            idx = self.neighVec.index(neighId)
            hc = self.neighHc[idx]
            if hc > hopc:
                self.neighHc[idx] = hopc
                
    
    def neighVecFull(self, n):
	 return len(self.neighVec) == n

    def addChannel(self, chId):
        self.chId = chId
        
    def contaminate(self, node):
       '''direct transmission'''
       self.addNeighborHc(node.id, 1) 
       
       '''epidemic dissemination with hop-count filtering'''
       for idx, nid in enumerate(node.neighVec): 
           hopc = node.neighHc[idx] + 1
           if hopc <= maxHopCount:
               self.addNeighborHc(nid, hopc)
    def disseminate(self, node):
        self.addNeighbor(node.id)
        for idx, nid in enumerate(node.neighVec):
            self.addNeighbor(nid)
                
    def aloha(self, pt):
        self.state = False
        p = rsObj.uniform(0,1)
        if p < pt:
            self.state = True        
#---------------------------------------------------------------------------    
#---------------------------------------------------------------------------
class SimNetwork:
    def __init__(self, numChannels, numNodes, numRepeat, tx_prob):
        self.numChannels = numChannels
        self.numNodes = numNodes
        self.nRuns = numRepeat
        self.tx_prob = tx_prob
        self.runTime = 0
        self.fracNodes=[]
        
    def createNodes(self, num_nodes):
        nodes = {}   
        for i in range(0, num_nodes):
            nodes[i] = Node(i, num_nodes)            
        return nodes  
         
    def addAdjMat(self, nodesL):
        n = len(nodesL)
        for node in nodesL.values():
            node.fillAdjMat(range(n)) #example of a clique...
    
    def allNodesFound(self,nodes):
        num_nodes = len(nodes)
        for node in nodes.values():
            
            if node.neighVecFull(num_nodes) == False:
                  return False
        '''Compute Fraction of nodes'''        
        return True
    def optTxProb(self, n, k):    
        a = n
        b = 2*k+n-1
        c = 1.*k
        f2 = (b - np.sqrt(b**2-4*a*c))/(2*a)
        return f2        
    def runSim(self):
        sum_time = 0.0
        time_vals = []
        nodesL = self.createNodes(self.numNodes)
        #self.addAdjMat(nodesL)  
 
        for runs in range(0, self.nRuns):
            
            time_slots = 0
            while self.allNodesFound(nodesL) == False:
                time_slots =  time_slots + 1
                
                channelIds =rsObj.randint(0, self.numChannels, len(nodesL))
                channelIds = channelIds.tolist() 
                for idx, node in enumerate(nodesL.values()):
                    channel = channelIds[idx]
                    node.addChannel(channel)
                    node.aloha(self.tx_prob)                                      
                ###end FOR
                for node in nodesL.values():
                    if node.state == False: #receiver
                        txList = []
                        for m in nodesL.values():
                            if m.id != node.id and node.chId == m.chId:
                                if m.state:
                                    txList.append(m)
                                #END IF
                            #END IF
                        ##END FOR
                        if len(txList) == 1:
                            node.disseminate(txList[0])
                    #end IF
            ###end While
            time_vals.append(time_slots)
            sum_time = sum_time + time_slots
            for node in nodesL.values():
                node.flushNeighbors()
        ### end for

        conf_int  = mean_conf_int(time_vals)
        std_discT = np.std(time_vals)
        '''add std value..'''
        
        #print 'mean:', conf_int, sum_time/self.nRuns,std_discT
        self.runTime = conf_int, round(std_discT,2)
#---------------------------------------------------------------------------
class Simulate(threading.Thread):
    def __init__(self, queue):
        threading.Thread.__init__(self)    
        self.queue = queue        
    
    def run(self):
        while True:
            objCoupon = self.queue.get()            
            objCoupon.runSim()
            self.queue.task_done()    
#---------------------------------------------------------------------------
class CouponThread(threading.Thread):
    def __init__(self, queue):
        threading.Thread.__init__(self)    
        self.queue = queue  
        
    def run(self):
        while True:
            objCoupon = self.queue.get()            
            objCoupon.aloha()
            
            self.queue.task_done()     
#---------------------------------------------------------------------------
numRepeat   = 500
#---------------------------------------------------------------------------
def funcProbOf(C, N):

    a = N
    b = 2*C+N-1
    c = C
    f2 = (b - np.sqrt(b**2-4*a*c))/(1.0*2*a)
    return f2
#--------------------------------------------------------------------------
def create_file(fname):
    csv_file = open(fname, 'wb')        
    return csv_file
#---------------------------------------------------------------------------
def create_filemode(fname, fmode):
    csv_file = open(fname, fmode)        
    return csv_file
#-------------------------------------------------------------------------
def runSim():
  
  myQueue       = Queue.Queue()
    
  for n in range(1,51):  
      myTh = Simulate(myQueue)      
      myTh.setDaemon(True)
      myTh.start()

  global couponQueue
  couponQueue = Queue.Queue()
  for cTh in range(0, 80):#vecChannels:
      myTh = CouponThread(couponQueue)
      myTh.setDaemon(True)
      myTh.start() 

  fileName = './MedalTh_Mar13_14.csv'      
  global csvFd, CSize, KSize
  csvFd = create_filemode(fileName, 'a')
  
  global maxHopCount
  
  maxHopCount = 15 
  
  numNodes = [2]+range(5,51,5)
  numChannels = [1,2,4,8,3,5,6,7]
  for channelID in [3,5,6,7,10,12,16]:#[8, 4, 2]:#vecChannels: 
      print 'Channel:',channelID, 'Time:', strftime("%a, %d %b %Y %H:%M:%S", localtime()), ' NT:', threading.active_count() 
      #for n, netSize in enumerate([10]):                    
      for n, netSize in enumerate(numNodes):      
             KSize = channelID 
             CSize = netSize                                    
                                    
             couponObj = {}
             pt = funcProbOf(channelID, netSize)
             couponObj[n] = SimNetwork(channelID, netSize, numRepeat, pt)             
                         
             for objC in couponObj.values():
                 myQueue.put(objC)
                 
             myQueue.join()

             for objSim in couponObj.values():
                 (dt_l, dt_m, dt_h), std_val = objSim.runTime
                 print '>(N, k, hop -> E[T], std)',
                 print  (netSize, channelID, maxHopCount, '->',dt_m, std_val)

                 csvFd.write(str(maxHopCount)+',')
                 csvFd.write(str(objSim.tx_prob)+',')
                 csvFd.write(str(channelID)+',')
                 csvFd.write(str(netSize)+',')
                 csvFd.write(str(dt_m)+',')
                 csvFd.write(str(dt_l)+',') 
                 csvFd.write(str(dt_h)+',')                 
                 csvFd.write(str(std_val)+'\n')
                 csvFd.flush()
             '''end of repetitions'''       
  ##--------------------------------------------------------------------
  csvFd.close()
#---------------------------------------------------------------------------
'''======================================================================'''
if __name__ == '__main__':
    
    runSim()
    
    print 'Execution finished.. : ', threading.active_count(),' at ', strftime("%a, %d %b %Y %H:%M:%S", localtime())
    
 
'''======================================================================'''
