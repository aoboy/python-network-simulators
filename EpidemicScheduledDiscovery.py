#!/usr/bin/env python
# -*- coding: utf-8 -*-
# make a horizontal bar chart
"""
Created on Fri Nov 30 10:53:52 2012

@author: gonga
"""


import matplotlib
matplotlib.use('Agg')

import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
#pdf = PdfPages('./pdfs/comp_channels_wlosses.pdf')

from pylab import *
import matplotlib.pyplot  as pyplot

import numpy as np
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
     
      self.state = False
      self.neighVec   = []

      self.neighVec.append(self.id)  
      
    
    def flushNeighbors(self):
       self.neighVec = []       
       self.neighVec.append(self.id)   
           
    def addNeighbor(self, neighId):
        if neighId not in self.neighVec:
            self.neighVec.append(neighId)
            #if self.id == 2: print neighId    
                  
    def neighVecFull(self):
	 return len(self.neighVec) == self.neighLen
  
    def addChannel(self, chId):
        self.chId = chId
               
    def disseminate(self, node):
        
        self.addNeighbor(node.id)
        
        for idx, nid in enumerate(node.neighVec):
            self.addNeighbor(nid)

                
    def transmit(self, pt):
        self.state = False
        p = rsObj.uniform(0,1)
        if p < pt:
            self.state = True  
       
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
#--------------------------------------------------------------------------
def create_file(fname):
    csv_file = open(fname, 'wb')        
    return csv_file
#---------------------------------------------------------------------------
def create_filemode(fname, fmode):
    csv_file = open(fname, fmode)        
    return csv_file
#---------------------------------------------------------------------------
def createNodes(n, T_period):
    lnodes={}
    for i in range(n):
        node = Node(i,n)
        lnodes[i]= node
    return lnodes
#---------------------------------------------------------------------------
def flushIds(Lnodes):
    for i, node in enumerate(Lnodes.values()):
        Lnodes[i].flushNeighbors()
#---------------------------------------------------------------------------
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
    plt.plot(sorted_v, yvals )   
    #plt.close()
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
  
  T_period = 20
  
  numRepeat   = 50
  
  fileName = './EpidemicScheduledDiscovery.csv'  
  csvFd = create_filemode(fileName, 'a')
  
  #nodes = [2, 4]+range(5, 51,5)
  channels=[2,4,8]  
  
  nodes = [10]#range(60, 101,10)
  channels=[3]#[4,8]
  
  
  if 1:
      print '----------------------------------------------'
      
      for num_nodes in nodes: 
          for channel in channels:
              pt = funcProbOf(channel, num_nodes)
    
              vec_times = []
              ##t_slot = 1
              Lnodes = createNodes(num_nodes, T_period)
              
              for t  in range(numRepeat):
                  #A=np.eye(n)
                  t_slot = 0
                  '''clear neighbors list'''
                  flushIds(Lnodes)
                  
                  while all2allComplete(Lnodes) == False:
                      t_slot = t_slot + 1 
                      
                      if t_slot%T_period == 0:
                      
                          '''each node uniformly at random selected a channel'''
                          rand_ch = rsObj.randint(0, channel, num_nodes).tolist()
                          
                          '''Transmission vector'''
                          tx_vec = rsObj.uniform(0,1, num_nodes).tolist()
                          
                          for ii, node in enumerate(Lnodes.values()):
                              '''Check if node is receiver'''
                              if tx_vec[ii] > pr:
                                  txList = []
                                  for jj, nodeX in enumerate(Lnodes.values()):
                                      if node.id != nodeX.id and rand_ch[ii] == rand_ch[jj]:
                                          if tx_vec[jj] < pt:
                                              txList.append(nodeX)
                                  if len(txList) == 1:
                                      Lnodes[ii].disseminate(txList[0])
                          '''END of TX'''   
                        else:
                            '''CHECK which type of slot it is..'''
                            
                                          
                          '''for each channel perform '''
                          for ii in range(channel):                          
                              nodes_ch = [jj for jj, kk in enumerate(rand_ch) if kk == ii]
                              #if t_slot%10 == 0:
                                  #print ii, nodes_ch
                              '''if more than one node that selected channel ii'''
                              if len(nodes_ch) > 1:
                                  '''who is a transmitter and who is a receiver...'''
                                  tv = rsObj.uniform(0,1,  len(nodes_ch)).tolist()
                                  '''index of transmitters'''
                                  sv = [jj for jj, pr in enumerate(tv) if pr < pt]
                                  '''guarantee that we have at least a Rx'''
                                  if len(sv) == 1: 
                                      for jj in nodes_ch:
                                          k      = sv[0]
                                          t_idx  = nodes_ch[k]
                                          t_node = Lnodes[t_idx]
                                          Lnodes[jj].disseminate(t_node)
                      else:
                          ''''''
                          for node in Lnodes.values():
                              ''''''
                              
                  #end While
                  vec_times.append(t_slot)
              #end for
              '''PLOT CDF...'''
              plot_cdf(vec_times)            
            
              discTime = mean(vec_times)#/(1-0.2)
              max_discTime = max(vec_times)
              std_discT= np.std(vec_times)
              discTime_l, discTime, discTime_h = mean_confidence_interval(vec_times)
              
              csvFd.write(str(n)+str(','))
              csvFd.write(str(channel)+str(','))
              csvFd.write(str(round(pt,3))+str(','))
              csvFd.write(str(round(discTime,2))+str(','))
              csvFd.write(str(round(discTime_l,2))+str(','))
              csvFd.write(str(round(discTime_h,2))+str(','))
              
              csvFd.write(str(round(max_discTime,2))+str(','))
              csvFd.write(str(round(std_discT,2))+str('\n'))
              csvFd.flush()
              
              print 'N =', num_nodes, channel, round(pt, 3), '->', round(discTime, 2), round(std_discT, 2)                                     
      #end for
  ###
  '''comment'''
  csvFd.close()
  ## 
#---------------------------------------------------------------------------
'''======================================================================'''
if __name__ == '__main__':
    
    runSim()      
    
 
'''======================================================================'''
