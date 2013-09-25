# graph.py
# James Mithen
# j.mithen@surrey.ac.uk

"""
Graph class to replace networkx module.  Note this class is modelled
on the graph implementation in the popular third party Python module
'networkx'. In particular, the function 'connected_nodes' is
essentially identical to the networkx function
'single_source_shortest_path_length'. See networkx module for more
information

CLASSES:
Graph - A simple graph class
"""

class Graph(object):
    """
    Very bare graph class.  Note this is only used to compute the
    connected components.  Be careful using it for anything else.
    """
    def __init__(self,nodes=None):
        """Initialise graph with a list of nodes."""
        
        self.node = {}
        self.adj = {}
        for n in nodes:
            self.node[n] = {}
            self.adj[n] = {}

    def __iter__(self):
        """Iterate over nodes when 'for i in G' is used."""
        
        return iter(self.node)

    def __getitem__(self,n):
        """Return dict of neighbours of node n when G[n] is used."""
        
        return self.adj[n]

    def add_edge(self,u,v):
        """Add an edge between two nodes u and v."""
        
        # add nodes u and/or v if do not already exist
        if u not in self.node:
            self.node[u] = {}
            self.adj[u] = {}
        if v not in self.node:
            self.node[v] = {}
            self.adj[v] = {}
        # add the edge
        ddict = self.adj[u].get(v,{})
        self.adj[u][v] = ddict
        self.adj[v][u] = ddict
        return

def connected_nodes(G,n):
    """
    Return dictionary of all nodes of graph G connected to node n
    (including node n itself).  This uses BFS.
    """
    
    seen = {}
    level = 0
    nextlevel = {n:1}
    while nextlevel:
        thislevel = nextlevel
        nextlevel = {}
        for v in thislevel:
            if v not in seen:
                seen[v] = level
                nextlevel.update(G[v])
        level = level + 1
    return seen

def connected_comps(G):
    """Return list of connected components of graph G, ordered by size."""
    
    seen = {}
    comps = []
    for n in G:
        if n not in seen:
            # get all nodes connected to n
            d = connected_nodes(G,n)
            comps.append(d.keys())
            seen.update(d)
    comps.sort(key=len,reverse=True)
    return comps
