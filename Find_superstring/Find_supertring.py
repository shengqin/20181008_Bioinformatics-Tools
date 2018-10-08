
# coding: utf-8

# ## PROBLEM : The greedy algorithm for fragment assembly 
# Write a function, `greedy_assemble`, that takes as input a list of read strings and uses the greedy fragment assembly algorithm to output a *single* superstring that contains all reads as substrings. You must use the graph-based (Hamiltonian path) version of the greedy algorithm. We will assume that:
# 1. we are assembling a single-stranded sequence and
# 2. that no read is a substring of any other read.
# 
# To keep things simple, for this homeowork, we will allow overlaps of any length (including zero).  In practice, for sequence assembly we would typically require some minimum overlap length.
# 
# For the purpose of making this algorithm deterministic, we must establish tiebreaking criteria for edges in the overlap graph that have the same weight. For two edges with the same weight, we will first choose the edge whose source vertex read is first in lexicographical order. If the source vertices are identical, then we choose the edge whose target vertex read is first in lexicographical order. For example, if e1 = ATCGGA → GGAT and e2 = ATCGGA → GGAA, we will attempt to use edge e2 first because GGAA < GGAT according to lexicographical order.

# In[1]:


class Vertex:
    """A vertex in directed graph that allows to add incoming and outcoming vertex."""
    
    def __init__(self, node):
        """Initiate the vertex with default values."""
        
        self.id = node
        self.invertex = ""
        self.outvertex = ""
        self.indegree = 0
        self.outdegree = 0

    def __str__(self):
        """Print the vertex with incoming and outcoming vertex."""
        
        instring = str(self.id) + ' incoming: ' + str(self.invertex)
        outstring = str(self.id) + ' outcoming: ' + str(self.outvertex)
        return instring + " " + outstring
    
    def add_invertex(self, invertex):
        """Record the incoming vertex."""
        
        self.invertex = invertex
        self.indegree +=1

    def add_outvertex(self, outvertex):
        """Record the outcoming vertex."""
        
        self.outvertex = outvertex
        self.outdegree +=1
        
    def get_invertex(self):
        """Get incoming vertex."""
        
        return self.invertex
    
    def get_outvertex(self):
        """Get outcoming vertex."""
        
        return self.outvertex

    def get_id(self):
        """Get the name of this vertex."""
        
        return self.id
        
    def has_outvertex(self):
        """Check if has outcoming vertex."""
        
        return self.outvertex != ""
        
    def has_invertex(self):
        """Check if has incoming vertex."""
        
        return self.invertex != ""
    
    def get_indegree(self):
        """Get the incoming degree."""
        return self.indegree
    
    def get_outdegree(self):
        """Get the outcoming degree."""
        return self.outdegree

class Graph:
    """A graph to store vertex and edges."""
    
    def __init__(self):
        """Initiate a graph with default fields."""
        
        self.vert_dict = {}
        self.num_vertices = 0
        self.num_edges = 0
        self.begin_vertex = ""

    def add_vertex(self, node):
        """Great a new vertex in the graph."""
        
        self.num_vertices += 1
        new_vertex = Vertex(node)
        self.vert_dict[node] = new_vertex

    def get_vertex(self, vertex):
        """Get a vertex."""
        
        if vertex in self.vert_dict:
            return self.vert_dict[vertex]
        else:
            return None

    def has_vertex(self, vertex):
        """Check if the graph has a vertex."""
        
        if vertex in self.vert_dict:
            return True
        else:
            return False
    
    def add_edge(self, frm, to):
        """Add an adge to a graph."""

        self.vert_dict[frm].add_outvertex(to)        
        self.vert_dict[to].add_invertex(frm)

        
        # record the first added edge for the convenience of tracking
        if self.num_edges == 0:
            self.begin_vertex = frm
        self.num_edges += 1
               
    def get_begin_vertex(self):
        """Get the recorded vertex."""
        return self.begin_vertex
    
    def get_num_edges(self):
        """Get the number of edges in the graph."""
        return self.num_edges

    def get_vertices(self):
        """Get a vertex in the graph."""
        return self.vert_dict.keys()


# In[2]:


def reads_to_list(reads_combination):
    """Find the overlap between two reads"""
    
    read_list = []
    for comb in reads_combination:
        overlap = overlap_length(comb[0],comb[1])
        read_list.append((comb[0],comb[1],-overlap,))
    return read_list


# In[3]:


def overlap_length(left, right):
    """Returns the length of the longest suffix of left that is a prefix of right"""
    
    i = - len(left)
    while i < 0:
        if right.startswith(left[i:]):
            break
        else:
            i = i+1
    return -i


# In[4]:


def check_circle(g, vertex, frm):
    """Check if adding the edge creates a circle"""
    
    while vertex.has_outvertex():
        out_vertex = vertex.get_outvertex()
        if out_vertex == frm:
            return True
        else:
            vertex = g.get_vertex(out_vertex)
    return False


# In[5]:


def add_edge(reads_list , g, reads_len):
    """Check if edge can be added to graph, if so, add the edge"""
    
    num_edge = 0
    while num_edge < reads_len-1:
        tup = reads_list.pop()
        frm = tup[0]
        to = tup[1]
        weight = tup[2]
        
        if g.get_vertex(frm).get_outdegree() == 0 and g.get_vertex(to).get_indegree() == 0:
            if g.get_num_edges() !=0 and check_circle(g, g.get_vertex(to), frm):
                continue    
            g.add_edge(frm, to)
            
            num_edge += 1

    return g


# In[6]:


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


# In[7]:


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
    
    # initiate
    reads_len = len(reads)
    g = Graph()
    super_list = [] 
    
    # add every read as vertex
    for read in reads: 
        g.add_vertex(read)
    
    # find possible edges and sort them
    reads_combination = list(itertools.permutations(reads,2))
    reads_list = reads_to_list(reads_combination) 
    reads_list.sort(key=lambda tup: tup[1],reverse=True )
    reads_list.sort(key=lambda tup: tup[0],reverse=True )
    reads_list.sort(key=lambda tup: tup[2],reverse=True )
    
    # add edge to list
    edged_g = add_edge(reads_list , g, reads_len) 

    # find the starting point of the edges
    one_vertex = g.get_vertex(edged_g.get_begin_vertex()) 
    while one_vertex.has_invertex():
        invertex = one_vertex.get_invertex()
        one_vertex = g.get_vertex(invertex)
    super_list.append(one_vertex.get_id())
    
    # add reads to superstring following the edges
    while one_vertex.has_outvertex():
        outvertex = one_vertex.get_outvertex()
        super_list.append(outvertex)
        one_vertex = g.get_vertex(outvertex)
    
    # merge and return
    final_list = merge_ordered_reads(super_list)
    return final_list


# Tests for `greedy_assemble` are provided at the bottom of this notebook.

# ## Practice: Assembling a small subset of Ebola virus reads

# Included with this notebook is the file `ebola_reads.txt` which is small subset of the Illumina reads used to assemble the genome of an isolate of the Ebola virus, which caused a major epidemic in West Africa. 
# 
# Use your greedy assemble function to assemble these reads. Once correctly assembled, these reads form a short segment of the genome of this virus. To allow your assembler to succeed, the reads have been cleaned of errors and have have been oriented so that they all come from the same strand of the genome.  You may find the following function below of use, which produces a list of reads from the contents of a file.

# In[8]:


def read_strings_from_file(filename):
    return [line.rstrip() for line in open(filename)]


# Once you have assembled the genomic segment, use the [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) web service to search the NCBI database of proteins with your assembled sequence. You should use BLASTX with its default settings. Based on the results of your BLASTX search, which gene is contained within this genomic segment?

# In[ ]:


greedy_assemble(read_strings_from_file("ebola_reads.txt"))


# #
# # *ANSWER *
# #
# - Answer: it is polymerase coded by AIE11805.1
# ![image.png](attachment:image.png)

# ### BONUS CHALLENGE
# A subset of Illumina reads from an Ebola virus genome sequencing experiment that cover the entire genome are included in the file `ebola_full_genome_reads.txt`. Can you get your code to assemble these reads in under 2 minutes? The code takes 59s to run. Wow!

# In[ ]:


greedy_assemble(read_strings_from_file("ebola_full_genome_reads.txt"))


# ### Tests for problem 1

# In[22]:


def test_greedy_assemble_with_files(reads_filename, superstring_filename):
    reads = read_strings_from_file(reads_filename)
    [superstring] = read_strings_from_file(superstring_filename)
    assert greedy_assemble(reads) == superstring 


# In[ ]:


# TEST: greedy_assemble returns a string
sanity_test_reads = read_strings_from_file("tests/test_reads.txt")
assert isinstance(greedy_assemble(sanity_test_reads), str)
print("SUCCESS: greedy_assemble returns a string passed!")


# In[ ]:


# TEST: greedy_assemble returns a superstring
def is_superstring(s, reads):
    return all(read in s for read in reads)
assert is_superstring(greedy_assemble(sanity_test_reads), sanity_test_reads)
print("SUCCESS: greedy_assemble returns a superstring passed!")


# In[ ]:


# TEST: greedy_assemble small test 1
small_test1_reads = ["GTT", "ATCTC", "CTCAA"]
assert greedy_assemble(small_test1_reads) == "ATCTCAAGTT"
print("SUCCESS: greedy_assemble small test 1 passed!")


# In[ ]:


# TEST: greedy_assemble small test 2
small_test2_reads = ["CGAAG", "ATCGA", "AGAG", "GGG"]
assert greedy_assemble(small_test2_reads) == "ATCGAAGAGGG"
print("SUCCESS: greedy_assemble small test 2 passed!")


# In[ ]:


# TEST: greedy_assemble small test 3
small_test3_reads = ["C", "T", "G", "A"]
assert greedy_assemble(small_test3_reads) == 'ACGT'
print("SUCCESS: greedy_assemble small test 3 passed!")


# In[ ]:


# TEST: greedy_assemble large test 1
test_greedy_assemble_with_files("tests/large_test1_reads.txt", "tests/large_test1_superstring.txt")
print("SUCCESS: greedy_assemble large test 1 passed!")

