
# coding: utf-8

# # BMI/CS 576 Fall 2018 - HW1
# The objectives of this homework are to practice
# 
# * with the basic algorithms for sequence assembly
# * reasoning about graphs and paths for the sequence assembly task
# 
# ## HW policies
# Before starting this homework, please read over the [homework polices](https://canvas.wisc.edu/courses/119014/pages/hw-policies) for this course.  In particular, note that homeworks are to be completed *individually*.

# ## PROBLEM 1: The greedy algorithm for fragment assembly (60 points)
# Write a function, `greedy_assemble`, that takes as input a list of read strings and uses the greedy fragment assembly algorithm to output a *single* superstring that contains all reads as substrings. You must use the graph-based (Hamiltonian path) version of the greedy algorithm. We will assume that:
# 1. we are assembling a single-stranded sequence and
# 2. that no read is a substring of any other read.
# 
# To keep things simple, for this homeowork, we will allow overlaps of any length (including zero).  In practice, for sequence assembly we would typically require some minimum overlap length.
# 
# For the purpose of making this algorithm deterministic, we must establish tiebreaking criteria for edges in the overlap graph that have the same weight. For two edges with the same weight, we will first choose the edge whose source vertex read is first in lexicographical order. If the source vertices are identical, then we choose the edge whose target vertex read is first in lexicographical order. For example, if e1 = ATCGGA → GGAT and e2 = ATCGGA → GGAA, we will attempt to use edge e2 first because GGAA < GGAT according to lexicographical order.

# In[933]:


class Vertex:
    def __init__(self, node):
        self.id = node
        self.invertex = ""
        self.outvertex = ""
        self.indegree = 0
        self.outdegree = 0

    def __str__(self):
        instring = str(self.id) + ' incoming: ' + str(self.invertex)
        outstring = str(self.id) + ' outcoming: ' + str(self.outvertex)
        return instring + " " + outstring
    
    def add_invertex(self, invertex, weight=0):
        self.invertex = invertex
        self.indegree +=1

    def add_outvertex(self, outvertex, weight=0):
        self.outvertex = outvertex
        self.outdegree +=1
        
    def get_invertex(self):
        return self.invertex
    
    def get_outvertex(self):
        return self.outvertex

    def get_id(self):
        return self.id

    def get_inweight(self, vertex):
        if vertex in self.invertex:
            return self.invertex[1]
        else:
            return None
    def get_outweight(self, vertex):
        if vertex in self.outvertex:
            return self.outvertex[1]
        else:
            return None
        
    def has_outvertex(self):
        return self.outvertex != ""
        
    def has_invertex(self):
        return self.invertex != ""
    
    def get_indegree(self):
        return self.indegree
    
    def get_outdegree(self):
        return self.outdegree

class Graph:
    def __init__(self):
        self.vert_dict = {}
        self.num_vertices = 0
        self.num_edges = 0
        self.begin_vertex = ""

    def __iter__(self):
        return iter(self.vert_dict.values())

    def add_vertex(self, node):
        self.num_vertices += 1
        new_vertex = Vertex(node)
        self.vert_dict[node] = new_vertex
        #return new_vertex

    def get_vertex(self, vertex):
        if vertex in self.vert_dict:
            return self.vert_dict[vertex]
        else:
            return None

    def has_vertex(self, vertex):
        if vertex in self.vert_dict:
            return True
        else:
            return False
    
    def add_edge(self, frm, to, weight = 0):
        if frm not in self.vert_dict:
            self.add_vertex(frm)
        if to not in self.vert_dict:
            self.add_vertex(to)

        self.vert_dict[frm].add_outvertex(to, weight)        
        self.vert_dict[to].add_invertex(frm, weight)
        
        if self.num_edges == 0:
            self.begin_vertex = frm
        
        self.num_edges += 1
            
    def get_begin_vertex(self):
        return self.begin_vertex
    
    def set_begin_vertex(self, frm):
        self.begin_vertex = frm
    
    def get_num_edges(self):
        return self.num_edges

    def get_vertices(self):
        return self.vert_dict.keys()


# In[934]:


def reads_to_list(reads_combination):
    read_list = []
    for comb in reads_combination:
        overlap = overlap_length(comb[0],comb[1])
        read_list.append((comb[0],comb[1],-overlap,))
    return read_list


# In[935]:


def overlap_length(left, right):
    """Returns the length of the longest suffix of left that is a prefix of right"""
    i = - len(left)
    while i < 0:
        if right.startswith(left[i:]):
            break
        else:
            i = i+1
    return -i


# In[936]:


def merge_ordered_reads(reads):
    """Returns the shortest superstring resulting from merging 
    the elements of reads, an ordered list of strings"""
    merged_read = []
    length = len(reads)
    if length == 0:
        return ""
    merged_read.append(reads[0])
    for i in range(length):
        if i+1 >= length:
            break
        else:
            reads_length = len(reads[i+1])
            overlap = overlap_length(reads[i],reads[i+1])
            merged_read.append(str(reads[i+1])[overlap:])
    return "".join(merged_read)


# In[937]:


def check_circle(g, vertex, frm):
    while vertex.has_outvertex():
        out_vertex = vertex.get_outvertex()
        if out_vertex == frm:
            return True
        else:
            vertex = g.get_vertex(out_vertex)
    return False


def add_edge(reads_list , g, reads_len):
    num_edge = 0
    while num_edge < reads_len-1:
        tup = reads_list.pop()
        frm = tup[0]
        to = tup[1]
        weight = tup[2]
        
        if g.get_vertex(frm).get_outdegree() == 0 and g.get_vertex(to).get_indegree() == 0:
            if g.get_num_edges() !=0 and check_circle(g, g.get_vertex(to), frm):
                continue    
            g.add_edge(frm, to, weight)
            
            num_edge += 1

    return g



# Code for PROBLEM 1
# You are welcome to develop your code as a separate Python module
# and import it here if that is more convenient for you.
import itertools
import queue
#from collections import defaultdict
        

def greedy_assemble(reads):
    """Returns a string that is a superstring of the input reads, which are given as a list of strings.
    The superstring computed is determined by the greedy algorithm as described in HW1, with specific tie-breaking
    criteria.
    """
    reads_len = len(reads)
    g = Graph()
    
    for read in reads:
        g.add_vertex(read)
        
    reads_combination = list(itertools.permutations(reads,2))
    reads_list = reads_to_list(reads_combination) # find the overlap between two reads and add to a list
    reads_list.sort(key=lambda tup: tup[1],reverse=True )
    reads_list.sort(key=lambda tup: tup[0],reverse=True )
    reads_list.sort(key=lambda tup: tup[2],reverse=True )
    
    #print(reads_list)
    edged_g = add_edge(reads_list , g, reads_len)
    super_list = []
    one_vertex = g.get_vertex(edged_g.get_begin_vertex())

    while one_vertex.has_invertex():
        #print(vertex)
        invertex = one_vertex.get_invertex()
        
        one_vertex = g.get_vertex(invertex)
    super_list.append(one_vertex.get_id())
    
    while one_vertex.has_outvertex():
        outvertex = one_vertex.get_outvertex()
        super_list.append(outvertex)
        one_vertex = g.get_vertex(outvertex)
    
    #print(super_list[0])
    final_list = merge_ordered_reads(super_list)
    print(final_list)
    return final_list
