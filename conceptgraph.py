# gensim modules
import logging
from operator import itemgetter
import sys
import networkx as nx
from datetime import datetime, timedelta
import matplotlib
import pdb
from collections import Counter
import math

import rake
import operator
import os.path

import collections

log = logging.getLogger()
log.setLevel(logging.DEBUG)

ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter(
             '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
log.addHandler(ch)

def compute_psis(N,t):
    psis = {}
    psis[N] = 1.
    for i in xrange(N-1,0,-1):
        psis[i] = psis[i+1]*t/(float(i+1.))+1.
    return psis


# the graph of events
class EventGraph():
    def __init__(self, filename):
        if os.path.isfile("events_propogation.gpickle"):
            log.info('Load graph after the propogation of words.')
            self.G = nx.read_gpickle("events_propogation.gpickle")
            self.run()
        elif os.path.isfile("events.gpickle"):
            log.info('Load graph.')
            self.G = nx.read_gpickle("events.gpickle")
            log.info('Propogation.')
            self.wordpropogation(self.G)
            log.info('Process graph.')
            self.run()
        else:
            self.construct(filename)
            self.wordpropogation(self.G)
            self.run()

    def construct(self, filename):
        log.info('Read raw GEXF file.')
        self.G = nx.read_gpickle(filename + ".gpickle")
        self.conceptDict = self.G.graph["conceptDict"]
        self.catgDict = self.G.graph["catgDict"]

        ndays = 10
        log.info("Connect the events containing same keywords published within %d days..." % ndays)
        for ccpt, eventlist in self.conceptDict.items():
            window = 50
            mylist = self.listOfDates(eventlist)
            for i in range(len(mylist)):
                for j in range(i+1, min(i+window, len(mylist))):
                    if datetime.strptime(mylist[j][1], '%Y-%m-%d') - \
                       datetime.strptime(mylist[i][1], '%Y-%m-%d') < timedelta(days=ndays):
                        self.G.add_edge(mylist[i][0], mylist[j][0])
                    else:
                        break;
        errornodes = []
        with open("./RAKE-tutorial/SmartStoplist.txt", "r") as stoppath:
            stopwords = stoppath.read().split() 
        for n in self.G:
            try:
                #w = self.G.node[n]["data"]["title"].split()\
                #    + self.G.node[n]["data"]["summary"].split()
                w = self.G.node[n]["data"]["title"].split()
                w = [_ for _ in w if not _ in stopwords and len(_) > 1]
                self.G.node[n]["words"] = Counter(w)
                self.G.node[n]["mixwords"] = self.G.node[n]["words"]
            except:
                errornodes.append(n)
        self.G.remove_nodes_from(errornodes)
        print "-" * 20

	nx.write_gpickle(self.G, "events.gpickle")
        log.info("Completed. Graph information:")
        print(nx.info(self.G))

    def subG(self):
        subnodes = set() 
        for ccpt in ["http://en.wikipedia.org/wiki/Donald_Trump"\
                    ,"http://en.wikipedia.org/wiki/Police"\
                    ,"http://en.wikipedia.org/wiki/Hillary_Rodham_Clinton"]: 
            for event, _ in self.conceptDict[ccpt]:
                subnodes.add(event)
        subG = self.G.subgraph(list(subnodes))                    
        print(nx.info(subG))
	nx.write_gpickle(subG, "subG.gpickle")

    def wordpropogation(self, subG):
        log.info("Words Propogation starts:")
        print(nx.info(subG))
        for _ in range(4):
            log.info("Round %d:" % _)
            #propagation
            log.info("propogation..")
            for n in subG: 
                subG.node[n]["mixwords"] = subG.node[n]["words"]
            for n in subG:
                for ch in subG.out_edges(n):
                    subG.node[ch[1]]["mixwords"] = \
                        self.mix(subG.node[ch[1]]["mixwords"], subG.node[n]["words"])
            for n in subG: 
                subG.node[n]["words"] = subG.node[n]["mixwords"]
    
            #inflation
            log.info("inflation..")
            for n in subG: 
                C = 0
                for word, count in subG.node[n]["words"].items():
                    subG.node[n]["words"][word] = math.pow( count, 1.2 )
                    C += math.pow( count, 1.2 )
                for word in subG.node[n]["words"].keys():
                    subG.node[n]["words"][word] /= C 
    
            #cutoff
            log.info("cutoff..")
            for n in subG: 
                for word, count in subG.node[n]["words"].items():
                    if count < 0.01: 
                        del subG.node[n]["words"][word]
	nx.write_gpickle(subG, "events_propogation.gpickle")

    def expand(self, n):
        ret = set([(n,0.1)]) 
        path = nx.single_source_shortest_path(self.G,n,cutoff = 3)
        for nn in path.keys():
            d = self.distance(self.G.node[n]["words"], self.G.node[nn]["words"])
            if d > 0.25:
                ret.add((nn,d))  
        return ret

    def investigate(self):
        stre = "4562441"
        print self.G.node[stre]["words"] 
        print self.G.node[stre]["data"]["title"] 
        print "================="
        for n in self.G[stre]:
            print self.G.node[n]["words"] 
            print self.G.node[n]["data"]["title"] 
            print "================="
        sys.exit(1)

    def difussion(self):
        node = "4562441"
        score = self.heatdiffusion([node])
        heatnodes = self.sweep(score)
        print heatnodes
        sys.exit(1)

    def run(self):
        #self.difussion()
        #self.investigate()

        degree_sequence = sorted(self.G.degree().items(), key=operator.itemgetter(1), reverse=True)
        eventlist = []
        blacklist = []
        log.info("Compute rankings by social scores..")
        #for n in self.G:
        for n, deg in degree_sequence[:250]:
            myscore = 0 
            if float(self.G.node[n]["data"]["score"]) > 1000 and not n in blacklist:
                for nn,d in self.expand(n):
                    myscore += self.G.node[nn]["data"]["score"] * d
                    blacklist += nn 
                myscore += 1000 * self.G.degree(n)
                eventlist.append((n,myscore))
        eventlist.sort(key=lambda tup: tup[1], reverse = True)
        log.info("Top 50 results...")
        for event in eventlist[:50]:
            print [_ for _,value in sorted(self.G.node[event[0]]["words"].items(), key=operator.itemgetter(1))]
            print "ID", n
            print "own:", float(self.G.node[event[0]]["data"]["score"]) 
            print "total:", event[1] 
            print "date:", self.G.node[event[0]]["data"]["date"]
            print "="*20

    def view(self, n):
        print self.G.node[n]["words"]
        print self.G.node[n]["data"]["summary"] 
        print "="*20
        for ch in self.G.out_edges(n):
            print self.G.node[ch[1]]["words"] 
            print self.G.node[ch[1]]["data"]["summary"] 
            print "-"*20

    def distance(self, d1, d2):
        ret = 0
        for key in set(d1.keys()).intersection(set(d2.keys())):
            ret += d1[key] * d2[key]
        return ret

    def mix(self, d1, d2):
        ret = {}
        d = self.distance(d1, d2)
        if d < 0.2: return d1
        for key in d1.keys() + d2.keys():
            if key in d1:
                if key in d2:
                    ret[key] = d2[key] + d1[key] 
                else:
                    ret[key] = d1[key] 
            else:
                ret[key] = d2[key]
        return ret

    def printTopConcept(self):
        conceptCount = []
        for ccpt, eventlist in self.conceptDict.items():
            conceptCount.append((ccpt,
                            sum([float(self.G.node[event]["data"]["score"])
                                for event,date in eventlist])))
        ccptList = sorted(conceptCount, key=itemgetter(1))[-50:]
        print(ccptList)

    # return a list of events, sorted according to the post date
    def listOfDates(self, eventlist):
        datelist = [(event, date)
                    for event, date in eventlist
                    if int(date.split('-')[0]) > 2013]
        datelist.sort(key=itemgetter(1))
        return datelist

    def heatdiffusion(self,seed): 
        ## Setup parameters that can be computed automatically
        N = 11 # see paper for how to set this automatically
        t = 5.
        eps = 0.01
        psis = compute_psis(N,t)
            
        ## Estimate hkpr vector 
        # G is graph as dictionary-of-sets, 
        # t, tol, N, psis are precomputed
        x = {} # Store x, r as dictionaries
        r = {} # initialize residual
        Q = collections.deque() # initialize queue
        for s in seed: 
          r[(s,0)] = 1./len(seed)
          Q.append((s,0))
        while len(Q) > 0:
          (v,j) = Q.popleft() # v has r[(v,j)] ...
          rvj = r[(v,j)]
          # perform the hk-relax step
          if v not in x: x[v] = 0.
          x[v] += rvj 
          r[(v,j)] = 0.
          mass = (t*rvj/(float(j)+1.))/(self.G.degree(v))
          for u in self.G[v]:   # for neighbors of v
            next = (u,j+1) # in the next block
            if j+1 == N: 
              x[u] += rvj/(self.G.degree(v))
              continue
            if next not in r: r[next] = 0.
            thresh = math.exp(t)*eps*(self.G.degree(u))
            thresh = thresh/(N*psis[j+1])
            if r[next] < thresh and \
               r[next] + mass >= thresh:
               Q.append(next) # add u to queue
            r[next] = r[next] + mass
        
        # Find cluster, first normalize by degree
        for v in x:
            x[v] = x[v]/self.G.degree(v)
            #print "hk[%2i] = %.16lf"%(v, x[v])
        return x

    def sweep(self, x):
        ## Step 2 do a sweep cut based on this vector 
        # now sort x's keys by value, decreasing
        sv = sorted(x.iteritems(), key=lambda x: x[1], reverse=True)
        S = set()
        volS = 0.
        cutS = 0.
        bestcond = 1.
        bestset = sv[0]
        for p in sv:
          s = p[0] # get the vertex
          volS += self.G.degree(s) # add degree to volume
          for v in self.G[s]:
            if v in S:
              cutS -= 1
            else:
              cutS += 1
          #print "v: %4i  cut: %4f  vol: %4f   conduct: %4f"%(s, cutS,volS, cutS/volS)
          S.add(s)
          if cutS/volS < bestcond:
            bestcond = cutS/volS
            bestset = set(S) # make a copy
        #print "Best set conductance: %f"%(bestcond)
        #print "  set = ", str(bestset)
        return bestset


if __name__ == "__main__":
    myG = EventGraph("mygraph")
    #myG.run()
