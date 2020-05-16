
# coding: utf-8

# # BIOINF 1 #

# In[1]:



'''PatternToNumber(Pattern)
Transform a pattern of DNA letters into an index number;
Convert k-mers into the 4^k different integers between 0 and 4^k − 1'''

#V1
def PatternToNumber1(Pattern):
    '''
    Transform a pattern of DNA letters into an index number;
Convert k-mers into the 4^k different integers between 0 and 4^k − 1
    
    Args:
        Pattern: nucleotide pattern to be turned into an index number.

    '''
    d={'A':0,'C':1,'G':2,'T':3}
    number=0
    for i in range(0,len(Pattern)):
        number += d[Pattern[i]]*(4**(len(Pattern)-1-i))
    return(number);
#--------------------------------------
# V2
def PatternToNumber(Pattern):
    '''
    Transform a pattern of DNA letters into an index number;
Convert k-mers into the 4^k different integers between 0 and 4^k − 1
    
    Args:
        Pattern: nucleotide pattern to be turned into an index number.

    '''
    d={'A':0,'C':1,'G':2,'T':3}
    if Pattern == '':
        return 0
    symbol = Pattern[-1]
    Pattern= Pattern[:-1]
    return(4 * PatternToNumber(Pattern) + d[symbol])

#---
#PatternToNumber('ATC')


# In[6]:

''' NumberToPattern(Number, k)'''

#V2
def NumberToPattern(Number, k):
    ''' 
    Takes a Number and return a k -length nucleotide (ATGC) pattern, very usefull for indexing
    
    Args:
        Number: the number to turn into a nucleotide pattern.
        k: length of nucleotide pattern

    '''
    d={'A':0,'C':1,'G':2,'T':3}
    pattern=[]
    if k==1:
        #return d.keys()[d.values().index(Number)]
        return (list(d.keys())[list(d.values()).index(int(Number))])
    prefixIndex = Number//4
    r = Number%4
    PrefixPattern = NumberToPattern(prefixIndex, (k-1))
    symbol = (list(d.keys())[list(d.values()).index(r)])
    return (PrefixPattern+symbol)

#print(NumberToPattern(5437, 8))


# In[1]:

'''PatternCount(Text, Pattern)

Count the frequency on a pattern(word) in a text(sequence)'''

def PatternCount(Text, Pattern):
    '''
    Count the frequency on a Pattern(word) in a Text(sequence)
    
     Args:
        Text: nucleotide sequence.
        Pattern: nucleotide pattern to be counted in the sequence
    
    '''
    count=0
    for n in range(0,(len(Text)-len(Pattern)+1)):
                if Text[n:n+len(Pattern)]== Pattern:
                    count += 1
    print(count)                
    return              
#PatternCount('AACTGGTCACACATC','AC')
#PatternCount('CGCGATACGTTACATACATGATAGACCGCGCGCGATCATATCGCGATTATC', 'CGCG')


# In[1]:

'''FrequentWords(Text, k)

Find all most frquent words of size k in a sequence'''

def FrequentWords(Text, k):
    '''
    Find all most frquent words of size k in a Text (DNA sequence)

    Args:
        Text: nucleotide/DNA sequence.
        k: size of most frequent words to be found in Text
        
    '''
    
    FrequentPattern=list()
    d={}
    for n in range(0,(len(Text)-(k-1))):
        word=Text[n:n+k]
        if word in d:
            d[word] += 1
        else:
            d[word] = 1     
    maxCount= max(d.values())

    max_keys = [k for k in d if d[k] == maxCount]  
    return(max_keys)     
          
#FrequentWords('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4) 
#print ('The script took '+str(time.time()-start_time)+' s')
#FrequentWords('TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT', 3) 


# In[3]:

'''FindingFrequentWordsBySorting(Text , k)'''


def FindingFrequentWordsBySorting(Text , k):
    '''
    Find all most frquent words of size k in a Text (DNA sequence) by sorting index lists

    Args:
        Text: nucleotide/DNA sequence.
        k: size of most frequent words to be found in Text
        
    '''
    
    from Bioinf_functions import NumberToPattern
    from Bioinf_functions import PatternToNumber
    
    FrequentPatterns = list()
    Index={}
    Count={}
    for n in range(len(Text)-(k-1)):
        Pattern = Text[n:n+k]
        Index[n] = int(PatternToNumber(Pattern))
        Count[n] = 1
    SortedIndex = sorted(Index.values())
    for i in range(1,len(Text)-(k-1)):
        if SortedIndex[i] == SortedIndex[i-1]:
            Count[i] = (Count[i-1]) + 1
    maxCount= max(Count.values())
    for i in range(0,len(Text)-(k-1)):
        if Count[i] == maxCount:
            Pattern = NumberToPattern(SortedIndex[(i)], k)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns

#FindingFrequentWordsBySorting('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4)  
#print ('The script took '+str(time.time()-start_time)+' s')
#this script was slower, but might become faster when dealing with long sequences


# In[6]:

'''ReverseComplement(sequence)
Find the reverse complement of a sequence'''

def ReverseComplement(sequence):
    '''
    Find the reverse complement of a sequence
    
    Args:
        sequence: Nucleotide sequence to return the reverse complement of.
        
    '''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(sequence))
    return(reverse_complement);


#ReverseComplement('ACGTTGCA')
#ReverseComplement('TTGTGTC')


# In[49]:

'''PatternMatch(Text,Pattern)
Finds all the positions of a pattern in a text'''

def PatternMatch(Text,Pattern):
    '''
    Finds all the positions of a Pattern in a Text
    
    Args:
        Text: nucleotide sequence.
        Pattern: Nucleotide pattern which positions in Text has to be found.
        
    '''
    output=list()
    for n in range(0,(len(Text)-(len(Pattern)-1))):
        word=Text[n:n+len(Pattern)]
        if word == Pattern:
            output.append(n)

    print(' '.join(map(str,output)))
    return;

#PatternMatch('GATATATGCATATACTT','ATAT')


# In[50]:

'''ClumpFind(Text,k,L,t)
Finds clump of at least t patterns of length k(kmers) within a window of length L from a Text(sequence)'''

def ClumpFind(Text,k,L,t):
    '''
    Finds clump of at least t patterns of length k(kmers) within a window of length L from a Text(sequence)
    
    Args:
        Text: nucleotide sequence.
        k: k-mer length.
        L: window length.
        t: min number of patterns in window.
        
    '''
    from collections import defaultdict
    kmer_dic = defaultdict(list)
    suc_keys = []    
    for n in range(len(Text)-(k-1)):
        kmer_dic[Text[n:n+k]].append(n)

    for key in kmer_dic.keys():
        if len(kmer_dic[key])>=t:
            index_bank = []
            for index in kmer_dic[key]:
                index_bank.append(index)

            for list_index in range(0,len(index_bank)-(t-1)):
                if index_bank[list_index+t-1]+k-index_bank[list_index]<=L:
                    suc_keys.append(key) 
                    break
                else:
                    continue        
        else:
            continue


    print (' '.join(map(str,suc_keys)))
    return;

#ClumpFind('CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA',5,50,4)


# In[17]:

'''FrequencyArray(Text , k)'''

def FrequencyArray(Text,k):
    '''
    Frequency array ....
    '''
    from Bioinf_functions import PatternToNumber
    
    FrequencyArray={}
    for i in range(4**k):
        FrequencyArray[i]= 0   
    for n in range(0,len(Text)-(k-1)):
        Pattern = Text[n:n+k]
        j=PatternToNumber(Pattern)
        FrequencyArray[j] += 1
    return(FrequencyArray.values());

#print(' '.join(map(str,(FrequencyArray('ATATATAC',2)))))


# In[1]:

'''Hamming Distance(p,q)

We say that position i in k-mers p1 … pk and q1 … qk is a mismatch if pi ≠ qi.
For example, CGAAT and CGGAC have two mismatches. 
The number of mismatches between strings p and q is called the Hamming distance between these strings
and is denoted HammingDistance(p, q).'''

def HammingDistance(p,q):
    '''
    We say that position i in k-mers p1 … pk and q1 … qk is a mismatch if pi ≠ qi.
    For example, CGAAT and CGGAC have two mismatches. 
    The number of mismatches between strings p and q is called the Hamming distance between these strings
    and is denoted HammingDistance(p, q).
    
    Args:
        p: DNA sequence 1.
        q: DNA sequence 2.
        
    '''
    HamDi=int(0)
    for i in range(len(p)):
        if p[i]==q[i]:
            continue
        else:
            HamDi+=+1
    return(HamDi)


# In[2]:

'''ApproximatePatternCount(Text,Pattern,d)

is a modify version of PatternMatch that account for a max of d mismatch in patterns'''

'''Our goal now is to modify our previous algorithm for the Frequent Words Problem in order to find DnaA boxes
by identifying frequent k-mers, possibly with mismatches. Given strings Text and Pattern as well as an integer d,
we define Countd(Text, Pattern) as the total number of occurrences of Pattern in Text with at most d mismatches.
For example, Count1(AACAAGCTGATAAACATTTAAAGAG, AAAAA) = 4 because AAAAA appears four times in this string with at
most one mismatch: AACAA, ATAAA, AAACA, and AAAGA. Note that two of these occurrences overlap.'''

#--
def ApproximatePatternCount(Text,Pattern,d):
    '''
    is a modify version of PatternMatch that account for a max of d mismatch in patterns.

    From Coursera course: Our goal now is to modify our previous algorithm for the Frequent Words Problem in order to find DnaA boxes
    by identifying frequent k-mers, possibly with mismatches. Given strings Text and Pattern as well as an integer d,
    we define Countd(Text, Pattern) as the total number of occurrences of Pattern in Text with at most d mismatches.
    For example, Count1(AACAAGCTGATAAACATTTAAAGAG, AAAAA) = 4 because AAAAA appears four times in this string with at
    most one mismatch: AACAA, ATAAA, AAACA, and AAAGA. Note that two of these occurrences overlap.
    
    Args:
        Text: DNA sequence
        Pattern: DNA pattern.
        d: Max number of mismatches.
        
    '''
    
    from Bioinf_functions import HammingDistance
    count=0
    for n in range(0,(len(Text)-(len(Pattern)-1))):
        word=Text[n:n+len(Pattern)]
        if HammingDistance(Pattern,word)<=d:
            count +=+1
    return(count)


# In[3]:

'''ImmediateNeighbors(Pattern)

generate the 1-neigborhood of Pattern'''

def ImmediateNeighbors(Pattern):
    '''
    Generate a list of all patterns that have a 1 nucletodide difference to Pattern (Neighbors)
    
    Args:
        Pattern: DNA pattern
    
    '''
    Neighborhood =list()
    nucl='ATCG'
    for i in range(0,len(Pattern)):
        symbol = Pattern[i]
        for n in nucl:               
            if Pattern[i]== n:
                       continue
            else:
                Neighbor = Pattern.replace(Pattern[i], n)
                Neighborhood.append(Neighbor)
    return(Neighborhood)


# In[4]:

'''Neighbors(Pattern, d)

Find all the near-patterns for which the haming distance with pattern is <= d
'''


def Neighbors(Pattern, d):
    '''
       Find all the near-patterns for which the haming distance with pattern is <= d
    
    Args:
        Pattern: DNA pattern.
        d: max number of mismatches allowed (max hamming distance).
        
    '''
    from Bioinf_functions import HammingDistance
    
    if d == 0:
        Neighborhood =(str(Pattern))
        return Neighborhood
    if len(Pattern) == 1: 
        Neighborhood=['A','C','G','T']
        return Neighborhood
    nucl=['A','C','G','T']
    Neighborhood=list()
    SuffixNeighbors = Neighbors((Pattern[1:]), d)
    for Text in SuffixNeighbors:
        if HammingDistance((Pattern[1:]), Text) < d:
            for n in nucl:
                Neighborhood.append(str(n + Text)) 
        else:
            Neighborhood.append(str(Pattern[0]+ Text))
    return Neighborhood


# In[5]:

'''FrequentWordsWithMismatches(Text, k, d)'''
#--------------------------------------


#--
def FrequentWordsWithMismatches(Text, k, d):
    '''
    Finds most frequent words of k nucleotides in Text, allowing d mismatches
    
    Args:
        Text: DNA sequence
        k: k-mer size
        d: max number of mismatches
    '''
    from Bioinf_functions import HammingDistance
    from Bioinf_functions import Neighbors

    dico={}
    for n in range(0,(len(Text)-(k-1))):
        Neighborhood = Neighbors(Text[n:n+k], d)
        for word in Neighborhood:
            if word in dico.keys():
                dico[word]+=1
            else:
                dico[word]= 1    

    max_value = max(dico.values())
    max_keys = [k for k in dico if dico[k] == max_value]
    return(max_keys)


# In[6]:


#--
def ApproximatePatternMatch(Text,Pattern,d):
    '''
    Missing description...
    '''
    
    from Bioinf_functions import HammingDistance
    
    posi=list()
    for n in range(0,(len(Text)-(len(Pattern)-1))):
        word=Text[n:n+len(Pattern)]
        if HammingDistance(Pattern,word)<=d:
            posi.append(n)
    return(posi)


# In[7]:

'''MOTIFENUMERATION(Dna, k, d)
A brute force algorithm for motif finding'''



#--

def MOTIFENUMERATION(Dna, k, d):
    '''
    A brute force algorithm for motif finding...
    '''
    from Bioinf_functions import HammingDistance
    from Bioinf_functions import Neighbors

    i=1
    p=list()
    while i<len(lines):
        list_i=list()
        for n in range(len(lines[i])-k+1):
            Neigh=Neighbors(lines[i][n:n+k],d)
            for word in Neigh:
                list_i.append(word)
        p.append(list_i)
        i=i+1  
    return(set.intersection(*map(set, p)))


# In[8]:

'''d(Pattern, Dna)
the sum of distances between Pattern and all strings in Dna'''


def d(Pattern, Dna):
   '''
   the sum of distances between Pattern and all strings in Dna...
   '''
   
   from Bioinf_functions import HammingDistance
   
   score=0
   k=len(Pattern)
   for i in range(0,len(Dna)):
       Hamdist={}
       minDist={}
       seq=Dna[i]
       for n in range(0,len(seq)-(k-1)):
           Hamdist[seq[n:n+k]]=HammingDistance(Pattern,seq[n:n+k])
       minimum=int(min(Hamdist.values()))
       score+=minimum
   return(score)


# In[9]:

'''MEDIANSTRING(Dna, k)
A k-mer Pattern that minimizes d(Pattern, Dna) among all k-mers Pattern.'''

def MEDIANSTRING(Dna, k):
    '''
    A k-mer Pattern that minimizes d(Pattern, Dna) among all k-mers Pattern.
    
    '''
    
    from Bioinf_functions import HammingDistance
    from Bioinf_functions import NumberToPattern
    from Bioinf_functions import d

    distance= float('inf')
    kmer_pattern=list()
    Median=list()
    for n in range(4**k):
        kmer_pattern.append(NumberToPattern(n, k))

    for Pattern in kmer_pattern:
        #if distance > d(Pattern, Dna):
        if 0 == d(Pattern, Dna):    
            
            print distance
            Median.append(Pattern)
    return(Median)


# In[10]:

'''Profile_most(Text,k,Profile)

Profile-most Probable k-mer:
a k-mer that was most likely to have been generated by Profile among all k-mers in Text'''

def Profile_most(Text,k,Profile):
    '''
    Profile-most Probable k-mer:
    a k-mer that was most likely to have been generated by Profile among all k-mers in Text

    '''
    
    from numpy import *
    from Bioinf_functions import PatternToNumber
    
    max_score=0
    most_probable_kmer=Text[0:0+k]
    for n in range(0,len(Text)-(k-1)):
        kmer=Text[n:n+k]
        score=1
        for i in range(0,len(kmer)):
            probaN=Profile[PatternToNumber(kmer[i])][i]
            return probaN
            score *= probaN
            return score
        if score>max_score:
            max_score=score
            most_probable_kmer=kmer

    return(most_probable_kmer)


# In[11]:

'''MakeProfile(Motif_list)
Make profile matrix from a list of motifs'''

#V1: Without pseudo counts, 0 probability possible

def MakeProfileRaw(Motif_list):
    '''
    Make profile matrix from a list of motifs without pseudo counts; 0 probability possible
    
    '''
    
    from numpy import *
    import numpy as np
    from Bioinf_functions import PatternToNumber
    
    k=len(Motif_list[0])
    matrix_occ = np.zeros((4,k))
    for i in range(len(Motif_list)):
        for n in range(len(Motif_list[i])):
            matrix_occ[(PatternToNumber(Motif_list[i][n])),n]+=1   
    Profile_matrix= matrix_occ / matrix_occ.sum(axis=0, dtype='float')       
    return Profile_matrix  

##V2: With pseudo counts, get rid of 0 probability

def MakeProfile(Motif_list):
    '''
    Make profile matrix from a list of motifs with pseudo counts; gets rind of 0 probability
    
    '''
    
    from numpy import *
    import numpy as np
    from Bioinf_functions import PatternToNumber

    k=len(Motif_list[0])
    matrix_occ = np.zeros((4,k))
    matrix_occ +=1
    for i in range(len(Motif_list)):
        for n in range(len(Motif_list[i])):
            matrix_occ[(PatternToNumber(Motif_list[i][n])),n]+=1   
    Profile_matrix= matrix_occ / matrix_occ.sum(axis=0, dtype='float')
    #return matrix_occ.sum(axis=0, dtype='float')
    return Profile_matrix  


# In[12]:

'''Consensus(Profile)

output consensus sequence from a profile matrix'''

def Consensus(Profile):
    '''
    output consensus sequence from a profile matrix...
    
    '''
    
    import numpy as np
    
    maxarg= list(Profile.argmax(0))
    Cons=''
    for posi in maxarg:
        Cons+=NumberToPattern(int(posi), 1)
    return Cons


# In[16]:

from Bioinf_functions import HammingDistance
from Bioinf_functions import NumberToPattern
from Bioinf_functions import PatternToNumber
from Bioinf_functions import Consensus
from Bioinf_functions import MakeProfile
from Bioinf_functions import d
from Bioinf_functions import Profile_most


def score(motifs):
    from collections import Counter
    # change to motifs as list of columns instead of list of rows
    matrix = [[] for k in range(len(motifs[0]))]
    for i in range(len(motifs[0])):
        for j in range(len(motifs)):
            matrix[i].append(motifs[j][i])
    score = 0
    for col in range(len(matrix)):
        # count letters in a column, delete most common letter, sum remaining
        c = Counter(matrix[col])
        del c[c.most_common(1)[0][0]]
        score += sum(c.values())
    return score

def Pr(kmer,Profile):
    '''
    return probability of kmer matching a profile
    
    '''
    
    score=1
    for n in range(len(kmer)):
        score *= Profile[n][kmer[n]]
    return score


def GREEDYMOTIFSEARCH(Dna,k,t):
    '''
    Description missing...
    
    '''
    
    from Bioinf_functions import MakeProfile
    from Bioinf_functions import Profile_most
    from Bioinf_functions import score
    
    
    #starting conditions:
    BestMotifs_list=list()
    for string in range(len(Dna)):
        BestMotifs_list.append(Dna[string][0:0+k])

    #list all kmers in first string
    kmer_1st_string=list()
    for n in range(len(Dna[0])-(k-1)):
        Motif=Dna[0][n:n+k]
        kmer_1st_string.append(Motif)
    #-     
    for motif in kmer_1st_string:
        Motif_list=list()    
        Motif_list.append(motif)
        for i in range(1,len(Dna)):
            Text=Dna[i]
            Profile=array(MakeProfile(Motif_list))
            new_motif=Profile_most(Text,k,Profile) 
            Motif_list.append(new_motif)
        Motif_score= score(Motif_list)
        if Motif_score < score(BestMotifs_list):
            BestMotifs_list=Motif_list
    return(BestMotifs_list)


# In[17]:

def RANDOMIZEDMOTIFSEARCH(Dna, k, t):
    '''
    missing description...
    
    '''
 
    from Bioinf_functions import MakeProfile
    from Bioinf_functions import Profile_most
    from Bioinf_functions import Consensus
    
    Motifs=list()
    for n in range(len(Dna)):
        rdm_posi=random.randint(0, (len(Dna[n])-(k)))
        Motifs.append(Dna[n][rdm_posi:rdm_posi+k])
    BestMotifs=Motifs
    #print BestMotifs

    while True:
        Profile=MakeProfile(Motifs)
        Motifs=list()
        for n in range(len(Dna)):
            Motifs.append(Profile_most(Dna[n],k,Profile))
        Motif_score= d(Consensus(MakeProfile(Motifs)), Motifs)
        if Motif_score < d(Consensus(MakeProfile(BestMotifs)), BestMotifs):
                BestMotifs=Motifs
        else:
            return BestMotifs,Motif_score


# In[18]:

'''Random(array)
this random number generator, denoted Random(array), where array=(p1, …, pn), 
models an n-sided biased die and returns integer i with probability pi'''
from Bioinf_functions import Pr
from Bioinf_functions import MakeProfile
from Bioinf_functions import PatternToNumber


def Random(Motif_list,Profile):
    '''
    this random number generator, denoted Random(array), where array=(p1, …, pn), 
    models an n-sided biased die and returns integer i with probability p
    '''
    
    from Bioinf_functions import Pr

    Array=list()
    for motif in Motif_list:
        Array.append(Pr(motif,Profile))
    PArray=[v/sum(Array) for v in Array]
    i = np.random.choice( len(PArray), p=PArray)
    return i


# In[2]:


def GIBBSSAMPLER(Dna, k, t, N):
    '''
    GIBBs Sampler...
    
    '''
    import random
    from Bioinf_functions import MakeProfile
    from Bioinf_functions import Random
    
    #dic for faster indexing of kmers
    dic={}
    for t in range(len(Dna)):
        dic[t]=list()
        for n in range(0,len(Dna[t])-(k-1)):
            dic[t].append(Dna[t][n:n+k])

    #get first set of Motifs 
    Motifs=list()
    for n in range(len(Dna)):
        rdm_posi=random.randint(0, (len(Dna[n])-(k)))
        Motifs.append(Dna[n][rdm_posi:rdm_posi+k])
        BestMotifs=Motifs[:]
        BestScore=score(BestMotifs)
    ###        
    x=0
    while x < N:
        i=(Random(Motifs,MakeProfile(Motifs)))
        Motifs.remove(Motifs[i])
        Profile_sans=MakeProfile(Motifs)
        ii=Random(dic[i],Profile_sans)    
        Motifs_i= dic[i][ii]  
        Motifs.insert(i,Motifs_i)
        sM=score(Motifs)
        if sM < BestScore:
            BestScore=sM
            BestMotifs=Motifs[:]
            x+=1
        else: 
            x+=1
    return BestMotifs, score(BestMotifs)


# # BIOINF 2 #

# In[ ]:

def KmerList(word,k):
    i=0
    liste=[]
    while i < len(word)-(k-1):
        kmer=word[i:i+k]
        liste.append(kmer)
        i+=1
    liste.sort()
    return (liste)

def DeBruijn_k(Patterns,k):
    from collections import defaultdict
    from collections import OrderedDict
    from Bioinf_functions import KmerList
    
    def Prefix(word):
        return word[:-1]

    def Suffix(word):
        return word[1:]
    
    pref={}
    suff={}
    link=defaultdict(list)
    liste=[]
    for word in Patterns:
        pref[word]=Prefix(word)
        suff[word]=Suffix(word)
    for word in Patterns:
        if word in sorted(pref):
            link[pref[word]].append(suff[word])      
    return [str(key)+' -> '+str(','.join(value)) for key,value in sorted(link.items())]


# In[3]:

def GenomePath2String(path):
    def Prefix(word):
        return word[:-1]

    def Suffix(word):
        return word[1:]
    seq=path[0]
    for n in range(1,len(path)):   
        if Prefix(path[n]) ==  Suffix(path[n-1]):
                seq+=path[n][-1]
    return seq  


# In[2]:

def eulerianCycle(graph):
    '''Find an Eulerian Cycle in a given graph
    #
    # Graph has to be strongly connected and balanced.
    #
    # @param graph:  A dict of lists, sources are dict keys, targets are in the lists'''
    
    # init cycle with any node
    cycle = [graph.keys()[0]]
    # run until all edges are moved from graph to cycle
    while len(graph) > 0:
        # whenever cycle closes rotate to a node with remaining targets
        if cycle[0] == cycle[-1]:
            while not cycle[0] in graph:
                cycle.pop(0)
                cycle.append(cycle[0])
        # the last node of cycle is the new source
        source = cycle[-1]
        # move one target at a time from graph to the end of cycle
        cycle.append(graph[source].pop())
        # clean up empty dict entries of graph
        if len(graph[source]) == 0: del graph[source]
    return cycle


# In[4]:

def EulerianPath(Graph):
    from collections import Counter
    from collections import defaultdict
    from Bioinf_functions import eulerianCycle
    # Find ends of graph
    outdegrees=Counter(dict([(k,(len(Graph[k]))) for k in Graph]))
    indegrees=Counter(value for values in Graph.itervalues() for value in values)
    for k in outdegrees+indegrees:
        if outdegrees[k]<indegrees[k]:
            pre=k
        if outdegrees[k]>indegrees[k]:
            post=k
    # Link ends together to balance graph        
    if pre in Graph:
        Graph[pre].append(post)
    else:
        Graph[pre]=[post]
        
    # run EulerianCycle    
    eul_cycle=eulerianCycle(Graph)[1:]
    
    # Rotate cycle and remove ends link
    while eul_cycle[0] != post or eul_cycle[-1] != pre:  
        eul_cycle.append(eul_cycle.pop(0))
    return eul_cycle


# In[6]:

def StringSpelledByGappedPatterns(GappedPatterns, k, d):
    ''' description '''
    from Bioinf_functions import GenomePath2String 
    FirstPatterns=[]
    for n in range(len(GappedPatterns)): # FirstPatterns long version
        FirstPatterns.append(GappedPatterns[n][0])
    SecondPatterns=[]
    for n in range(len(GappedPatterns)):
        SecondPatterns.append(GappedPatterns[n][1])
        
    PrefixString=GenomePath2String(FirstPatterns)
    SuffixString=GenomePath2String(SecondPatterns)
    for i in range(int(k+d+1), len(PrefixString)):
        if PrefixString[i] != SuffixString[i - k - d]:
            return "there is no string spelled by the gapped patterns"
    return PrefixString+SuffixString[-(k+d):] 


# In[1]:

def DeBruijn_k_pairs(GappedPatterns,k):
    ''' Construct a deBruijn graph for a collection of paired reads'''
    
    from collections import defaultdict
    from collections import OrderedDict
    
    def Prefix_pair(word):
        return '|'.join([read[:-1] for read in word.split('|')])

    def Suffix_pair(word):
        return '|'.join([read[1:] for read in word.split('|')])
    
    pref={}
    suff={}
    link=defaultdict(list)
    liste=[]
    for word in GappedPatterns:
        pref[word]=Prefix_pair(word)
        suff[word]=Suffix_pair(word)
    for word in GappedPatterns:
        if word in sorted(pref):
            link[pref[word]].append(suff[word])      
    return [str(key)+' -> '+str(','.join(value)) for key,value in sorted(link.items())]


# In[1]:

def CircularSpectrum_string(Peptide, AminoAcid, AminoAcidMass):
    ''' generate circular spectrum of peptide
    when peptide input is under the string form:
    EXAMPLE: peptide NQEL = [114,128,129,113] or 'NQEL' '''
        
    PrefixMass={}
    PrefixMass[0]=0
    for i in range(1,len(Peptide)+1):
        for j in range(len(AminoAcid)):
            if AminoAcid[j] == Peptide[i-1]:
                PrefixMass[i] = PrefixMass[i-1] + AminoAcidMass[j]
    peptideMass = PrefixMass[len(Peptide)]

    CircularSpectrum=[0]
    for i in range(len(Peptide)):
        for j in range(i+1,len(Peptide)+1):
            CircularSpectrum.append(PrefixMass[j]-PrefixMass[i])
            if i > 0 and j < len(PrefixMass)-1:
                CircularSpectrum.append(peptideMass -(PrefixMass[j]-PrefixMass[i]))
    return sorted(CircularSpectrum)


def CircularSpectrum_list(Peptide, AminoAcid, AminoAcidMass):
    ''' generate circular spectrum of peptide
    when peptide input is under the list form:
    EXAMPLE: peptide NQEL = [0,114,128,129,113] '''

    from collections import defaultdict
    PrefixMass={}
    PrefixMass[0]=0
    for i in range(1,len(Peptide)):
        for j in range(len(AminoAcid)):
            if AminoAcid[j] == Peptide[i]:
                PrefixMass[i] = PrefixMass[i-1] + AminoAcidMass[j]
    peptideMass = PrefixMass[len(Peptide)-1]

    CircularSpectrum=[0]
    for i in range(len(Peptide)-1):
        for j in range(i+1,len(Peptide)):
            CircularSpectrum.append(PrefixMass[j]-PrefixMass[i])
            if i > 0 and j < len(PrefixMass)-1:
                CircularSpectrum.append(peptideMass -(PrefixMass[j]-PrefixMass[i]))
    return sorted(CircularSpectrum)
    


# In[ ]:


def LinearSpectrum_string(Peptide, AminoAcid, AminoAcidMass):
    ''' generate circular spectrum of peptide
    when peptide input is under the string form:
    EXAMPLE: peptide NQEL = [114,128,129,113] or 'NQEL' '''
    PrefixMass={}
    PrefixMass[0]=0
    for i in range(len(Peptide)):
        for j in range(len(AminoAcid)):
            if AminoAcid[j] == Peptide[i]:
                PrefixMass[i+1] = PrefixMass[i] + AminoAcidMass[j]
    LinearSpectrum=[0]
    for i in range(len(Peptide)):
        for j in range(i+1,len(Peptide)+1):
            LinearSpectrum.append(PrefixMass[j]-PrefixMass[i])
    return sorted(LinearSpectrum)

def LinearSpectrum_list(Peptide, AminoAcid, AminoAcidMass):
    ''' generate circular spectrum of peptide
    when peptide input is under the list form:
    EXAMPLE: peptide NQEL = [0,114,128,129,113] '''
    PrefixMass={}
    PrefixMass[0]=0
    for i in range(1,len(Peptide)):
        for j in range(len(AminoAcid)):
            if AminoAcid[j] == Peptide[i]:
                PrefixMass[i] = PrefixMass[i-1] + AminoAcidMass[j]
    LinearSpectrum=[0]
    for i in range(len(Peptide)-1):
        for j in range(i+1,len(Peptide)):
            LinearSpectrum.append(PrefixMass[j]-PrefixMass[i])
    return sorted(LinearSpectrum)


# In[2]:

class StringPeptide:

    def __init__(self, x):
        self.string = x
        self.description = "Peptide represented by a string of letters for each aa"
        self.AminoAcid='GASPVTCILNDKQEMHFRYW'
        self.AminoAcidMass=[57,71,87,97,99,101,103,113,113,114,115,128,128,129,131,137,147,156,163,186]
        
    def CircularSpectrum(self):
        ''' generate circular spectrum of peptide
        when peptide input is under the string form:
        EXAMPLE: peptide NQEL = [114,128,129,113] or 'NQEL' '''
             
        PrefixMass={}
        PrefixMass[0]=0
        for i in range(1,len(self.string)+1):
            for j in range(len(self.AminoAcid)):
                if self.AminoAcid[j] == self.string[i-1]:
                    PrefixMass[i] = PrefixMass[i-1] + self.AminoAcidMass[j]
        peptideMass = PrefixMass[len(self.string)]
        CircularSpectrum=[0]
        for i in range(len(self.string)):
            for j in range(i+1,len(self.string)+1):
                CircularSpectrum.append(PrefixMass[j]-PrefixMass[i])
                if i > 0 and j < len(PrefixMass)-1:
                    CircularSpectrum.append(peptideMass -(PrefixMass[j]-PrefixMass[i]))
        return sorted(CircularSpectrum)
    
    def LinearSpectrum(self):
        ''' generate circular spectrum of peptide
        when peptide input is under the string form:
        EXAMPLE: peptide NQEL = [114,128,129,113] or 'NQEL' '''
        PrefixMass={}
        PrefixMass[0]=0
        for i in range(len(self.string)):
            for j in range(len(self.AminoAcid)):
                if self.AminoAcid[j] == self.string[i]:
                    PrefixMass[i+1] = PrefixMass[i] + self.AminoAcidMass[j]
        LinearSpectrum=[0]
        for i in range(len(self.string)):
            for j in range(i+1,len(self.string)+1):
                LinearSpectrum.append(PrefixMass[j]-PrefixMass[i])
        return sorted(LinearSpectrum)

    def CircularPeptideScore(self,Spectrum):
        ''' score of cyclic peptide against a spectrum'''
        floatSpectrum=Spectrum[:]
        score=0
        for spike in self.CircularSpectrum():
            if spike in floatSpectrum:
                floatSpectrum.remove(spike)
                score+=1
        return score        

    def LinearPeptideScore(self,Spectrum):
        ''' score of cyclic peptide against a spectrum'''
        floatSpectrum=Spectrum[:]
        score=0
        for spike in self.LinearSpectrum():
            if spike in floatSpectrum:
                floatSpectrum.remove(spike)
                score+=1
        return score


# In[1]:

class ListPeptide:
        
    def __init__(self, x):
        self.string = x
        self.description = "Peptide represented by a string of letters for each aa"
        self.AminoAcid=[57,71,87,97,99,101,103,113,114,115,128,129,131,137,147,156,163,186]
        self.AminoAcidMass=[57,71,87,97,99,101,103,113,114,115,128,129,131,137,147,156,163,186]

        
    def CircularSpectrum(self):
        ''' generate circular spectrum of peptide
        when peptide input is under the list form:
        EXAMPLE: peptide NQEL = [0,114,128,129,113] '''

        from collections import defaultdict
        PrefixMass={}
        PrefixMass[0]=0
        for i in range(1,len(self.string)):
            for j in range(len(self.AminoAcid)):
                if self.AminoAcid[j] == self.string[i]:
                    PrefixMass[i] = PrefixMass[i-1] + self.AminoAcidMass[j]
        peptideMass = PrefixMass[len(self.string)-1]
        CircularSpectrum=[0]
        for i in range(len(self.string)-1):
            for j in range(i+1,len(self.string)):
                CircularSpectrum.append(PrefixMass[j]-PrefixMass[i])
                if i > 0 and j < len(PrefixMass)-1:
                    CircularSpectrum.append(peptideMass -(PrefixMass[j]-PrefixMass[i]))
        return sorted(CircularSpectrum)
    
    def LinearSpectrum(self):
        ''' generate circular spectrum of peptide
        when peptide input is under the list form:
        EXAMPLE: peptide NQEL = [0,114,128,129,113] '''
        PrefixMass={}
        PrefixMass[0]=0
        for i in range(1,len(self.string)):
            for j in range(len(self.AminoAcid)):
                if self.AminoAcid[j] == self.string[i]:
                    PrefixMass[i] = PrefixMass[i-1] + self.AminoAcidMass[j]
        LinearSpectrum=[0]
        for i in range(len(self.string)-1):
            for j in range(i+1,len(self.string)):
                LinearSpectrum.append(PrefixMass[j]-PrefixMass[i])
        return sorted(LinearSpectrum)
    
    def CircularPeptideScore(self,Spectrum):
            ''' score of cyclic peptide against a spectrum'''
            floatSpectrum=Spectrum[:]
            score=0
            for spike in self.CircularSpectrum():
                if spike in floatSpectrum:
                    floatSpectrum.remove(spike)
                    score+=1
            return score        

    def LinearPeptideScore(self,Spectrum):
        ''' score of cyclic peptide against a spectrum'''
        floatSpectrum=Spectrum[:]
        score=0
        for spike in self.LinearSpectrum():
            if spike in floatSpectrum:
                floatSpectrum.remove(spike)
                score+=1
        return score 


# In[3]:

def Expand(Peptides,AminoAcidNames):
       '''
       Expand each peptides by adding one aa
       '''
       NewPeptides=[]
       for n in range(len(Peptides)):
           peptide=Peptides[n]
           for aa in AminoAcidNames:
               newpeptide=peptide[:]
               newpeptide.append(aa)
               NewPeptides.append(newpeptide)
       return NewPeptides


# In[4]:

def Trim(Leaderboard, Spectrum, N,AminoAcid, AminoAcidMass):
    import pandas as pd
    from Bioinf_functions import ListPeptide
    PepScoredic = {}
    
    for j in range(len(Leaderboard)):
        Peptide = ListPeptide(Leaderboard[j])
        Peptide.AminoAcid=AminoAcid[:]
        Peptide.AminoAcidMass=AminoAcidMass[:]
        PepScoredic['-'.join(map(str,Peptide.string))]=Peptide.LinearPeptideScore(Spectrum)
    df=pd.DataFrame(list(PepScoredic.iteritems()),columns=['Peptide','Score'])
    df=df.sort(columns='Score',ascending=False)
    df.index = range(len(df)) #reset index
    for j in range (N,len(Leaderboard)):
        if df.Score[j]< df.Score[N-1]:
            df2=df[df.index < j]
            break
            
    return [map(int,item.split('-')) for item in list(df[df.index < j].T.itertuples())[0][1:]]


# In[1]:

def LEADERBOARDCYCLOPEPTIDESEQUENCING(N,Spectrum,AminoAcidNames,AminoAcidMass):
    from Bioinf_functions import ListPeptide,Expand,Trim

    Leaderboard=[[0]]
    LeaderPeptide=[[0]]
    LeaderPeptideScore=0

    while len(Leaderboard)>0:
        Leaderboard=Expand(Leaderboard,AminoAcidNames)
        for peptide in Leaderboard[:]:
            
            Peptide= ListPeptide(peptide)
            Peptide.AminoAcid=AminoAcidNames[:]
            Peptide.AminoAcidMass=AminoAcidMass[:]
            
            if Peptide.LinearSpectrum()[-1] == Spectrum[-1]:
                    if Peptide.CircularPeptideScore(Spectrum)>LeaderPeptideScore:
                        LeaderPeptide=Peptide.string
                        LeaderPeptideScore=Peptide.CircularPeptideScore(Spectrum)
                        Leaderboard.remove(LeaderPeptide)
            elif sum(Peptide.string) > Spectrum[-1]:
                Leaderboard.remove(Peptide.string)
        try:
            Leaderboard = Trim(Leaderboard, Spectrum, N,AminoAcidNames, AminoAcidMass)
        except ValueError:
            #print 'End of pipeline'
            return '-'.join(map(str,LeaderPeptide[1:]))


# In[1]:

### BIOINF 3 ###




# # BIOINF 4 #

# In[7]:

def delta(i,k):
    if i==k:
        return 0
    else:
        return 1


# In[6]:

def SmallParsimony_rooted(T,seqs):
    '''Small Parsimony Problem. Find the most parsimonious labeling of the internal nodes of a rooted tree.
     Input: A rooted binary tree with each leaf labeled by a string of length m.
     T={0:[4],1:[4],2:[5],3:[5],4:[0,1,6],5:[2,3,6],6:[4,5]}
     seqs={0:'CAAATCCC',1:'ATTGCGAC',2:'CTGCGCTG',3:'ATGGACGA',4:'',5:'',6:''}
     Output: A labeling of all other nodes of the tree by strings of length m that minimizes the parsimony score of the tree
     tree_score=16
     Character={0: 'CAAATCCC', 1: 'ATTGCGAC', 2: 'CTGCGCTG', 3: 'ATGGACGA', 4: 'ATAGCCAC', 5: 'ATGGACTA', 6: 'ATAGACAA'}
     '''
    from Bioinf_functions import delta
    ### find sequence length ###
    for k,v in seqs.items():
        if len(v)>0:
            sequence_len=len(v)
            break
    ### Init sequence dict, and tree score
    Character=seqs.copy()
    Character = Character.fromkeys(Character, '')
    tree_score=0
    ### proceed nucleotibe by nucleotide ###
    for n in range(sequence_len):
        
        leaves=[k for k,v in T.items() if len(T[k])==1]
        Tag={} #marks which nodes were visited
        s={'A':{},'C':{},'G':{},'T':{}}
        minscore_nucl={}
        #### initiate ####
        for v in T.keys():
            Tag[v]=0
            if v in leaves:
                Tag[v]=1 #mark node as visited
                for k in ['A','C','G','T']:
                    if k==seqs[v][n]:
                        s[k].update({v:0})
                        minscore_nucl[v]=k
                    else:
                        s[k].update({v:float('Inf')})
        #### calculate parsimony score  ####
        ripe_nodes=[k for k,v in Tag.items() if v==0] #find nodes not visited yet
        while len(ripe_nodes)>0:
            for v in ripe_nodes:
                Tag[v]=1
                for k in ['A','C','G','T']:
                    s[k].update({v:sum([min([s[i][item]+delta(i,k) for i in s.keys() ]) for item in T[v] if item<=v])})
                minscore=min([s[k][v] for k in s.keys()])
                minscore_nucl[v]=[key for key,val in s.items() if s[key][v]==minscore ]
            ripe_nodes=[k for k,v in Tag.items() if v==0]
        #### backtrack ####
        root=max(T.keys()) #last common ancester/root
        minval=min([s[k][root] for k in s.keys()]) #lowest score for last common ancestor
        for key in s.keys():
            if s[key][root] ==minval:
                Ancestral_char=key # Ancestral character
                break
        tree_score+=minval
        Character[root] += Ancestral_char
        queue = []
        queue.append(root)
        while queue:
            parent_node = queue.pop()
            son = T[parent_node][0]
            daughter = T[parent_node][1]
            parent_Char = Character[parent_node][n]
            if parent_Char in minscore_nucl[son]:
                Character[son]+=parent_Char
            else: 
                Character[son]+=minscore_nucl[son][0]
            if parent_Char in minscore_nucl[daughter]:
                Character[daughter]+=parent_Char
            else: 
                Character[daughter]+=minscore_nucl[daughter][0]
            if len(T[son])>1:
                queue.append(son)
            if len(T[daughter])>1:
                queue.append(daughter)
    return tree_score,Character


# In[4]:

def SmallParsimony_unrooted(T,seqs):
    '''
    INPUT
    T =  {0: [4], 1: [4], 2: [5], 3: [5], 4: [0, 1, 5], 5: [2, 3, 4]}
    seqs =  {0: 'TCGGCCAA', 1: 'CCTGGCTG', 2: 'CACAGGAT', 3: 'TGAGTACC', 4: '', 5: ''}
    OUTPUT
    score=17
    sequence={0: 'TCGGCCAA', 1: 'CCTGGCTG', 2: 'CACAGGAT', 3: 'TGAGTACC', 4: 'CCTGGCAA', 5: 'CAAGGAAC'}
    '''
    from Bioinf_functions import delta,SmallParsimony_rooted
    #### find two internal nodes
    for k,v in T.items():
        if len(v)>1:
            node0=k
            for item in v:
                if len(T[item])>1:
                    node1=item
                    break
            
    #### turn unrooted tree to rooted tree by adding root between the two internal nodes
    new_node=len(T)
    node0_newlinks=[item for item in T[node0] if item!= node1]+[new_node]
    node1_newlinks=[item for item in T[node1] if item!= node0]+[new_node]

    Tree=T.copy()
    seqs.update({new_node:''})
    Tree.update({node0:node0_newlinks,node1:node1_newlinks,new_node:[node0,node1]})

    #### run SmallParsimony

    score,sequence=SmallParsimony_rooted(Tree,seqs)
    sequence.pop(new_node)
    return score,sequence


# In[9]:

def NearestNeighborsofTree(node0,node1,T):
    '''
    INPUT
    node0=5
    node1=4
    T={0: [4], 1: [4], 2: [5], 3: [5], 4: [0, 1, 5], 5: [2, 3, 4]}
    OUTPUT
    T2={0: [5], 1: [4], 2: [4], 3: [5], 4: [1, 5, 2], 5: [3, 4, 0]}
    T3={0: [4], 1: [5], 2: [4], 3: [5], 4: [0, 5, 2], 5: [3, 4, 1]}
    '''
    
    #### find subtrees and swap them 
    a,b=[item for item in T[node0] if item!= node1]
    c,d=[item for item in T[node1] if item!= node0]
    T2=T.copy()
    T2.update({node0:[item for item in T2[node0] if item!= a]+[c],node1:[item for item in T2[node1] if item!= c]+[a],a:[item for item in T2[a] if item!= node0]+[node1],c:[item for item in T2[c] if item!= node1]+[node0]})
    T3=T.copy()
    T3.update({node0:[item for item in T3[node0] if item!= a]+[d],node1:[item for item in T3[node1] if item!= d]+[a],a:[item for item in T3[a] if item!= node0]+[node1],d:[item for item in T3[d] if item!= node1]+[node0]})
    return T2,T3


# # BIOINF 5 #

# # BIOINF 6 #

# In[1]:

def TrieConstruction(Patterns):
    ''' Solve the Trie Construction Problem.
    INPUT:
    Patterns=['ATA','ATAGA','ATC','GAT']
    OUTPUT:
    {0: {1: 'A', 7: 'G'}, 1: {2: 'T'}, 2: {3: 'A', 6: 'C'}, 3: {4: 'G'}, 4: {5: 'A'}, 5: {}, 6: {}, 7: {8: 'A'}, 8: {9: 'T'}, 9: {}}
    '''
    Trie={}
    root=0
    Trie[root]={}
    newNode=1
    for pat in Patterns:
        currentNode=root
        for i in range(len(pat)):
            currentSymbol= pat[i]
            if currentSymbol in Trie[currentNode].values():
                endingNode=[k for k,v in Trie[currentNode].items() if v==currentSymbol][0]
                currentNode=endingNode
                continue
            else:
                Trie[currentNode].update({newNode:currentSymbol})
                Trie[newNode]={}
                currentNode=newNode
                newNode+=1
    return Trie


# In[2]:

def TrieMatching(Text, Trie):
    ''' finds whether any strings in Patterns match a prefix of Text.
    INPUT:
    Text='AATCGGGTTCAATCGGGGT'
    Patterns=['ATCG','GGGT']
    OUTPUT= [1, 4, 11, 15] # indices of all starting positions in Text where a string from Patterns appears as a substring.
    '''
    
    def PrefixTrieMatching(Text, Trie):
        v=0
        n=0
        symbol = Text[0]
        path=[0]
        while True:
            if len(Trie[v])==0:
                return path
            elif symbol in Trie[v].values():
                w=[key for key,val in Trie[v].items() if val==symbol][0]
                path+=[w]
                v = w
                n+=1
                try:
                    symbol = Text[n]
                except IndexError:
                    pass
            else:
                return "no matches found"
    #####################
    
    matches=[]
    posi=0
    while len(Text)>0:
        match=PrefixTrieMatching(Text, Trie)
        if match !='no matches found':
            matches+=[posi]
        posi+=1
        Text=Text[1:]
    return matches


# In[3]:

def ModifiedSuffixTreeConstruction(Text):
    '''
    Construct the Suffix Tree of a String
    INPUT:
    Text='ATAAATG$'
    OUTPUT= {0: {1: [0, 1], 9: [1, 1], 29: [6, 2], 30: [7, 1]}, 1: {16: [3, 1], 2: [1, 1]},
    2: {8: [2, 6], 25: [6, 2]}, 8: 0, 9: {27: [6, 2], 15: [2, 6]}, 15: 1, 16: {20: [4, 4], 
    23: [5, 3]}, 20: 2, 23: 3, 25: 4, 27: 5, 29: 6, 30: 7}
    # suffix tree in dictionary with the following format:
    #{nodein: {nodeout:[<initial position of the substring to which it corresponds in Text>, <length of this substring>]}}
    '''
    
    def ModifiedSuffixTrieConstruction(Text):
        root=0
        Trie = {}
        Trie[root]={}
        newNode=root
        for i in range(len(Text)):
            currentNode = root
            for j in range(i,len(Text)):
                currentSymbol = Text[j]
                if currentSymbol in [a for a,b in Trie[currentNode].values()]:
                    currentNode = [k for k,v in Trie[currentNode].items() if v[0]==currentSymbol][0] #ending node of this edge
                else:
                    newNode+=1
                    Trie[currentNode].update({newNode:[currentSymbol,j]})
                    Trie[newNode]={}
                    currentNode = newNode
            if currentSymbol == '$': 
                Trie[newNode]=i
        return Trie

    def MaximalNonBranchingPaths(graph):
        sets=[]
        for k,v in graph.items():
            if type(v)==dict:
                sets+=[[k,[out for out in v.keys()]]]
        d = {} 
        for node,outs in sets:
            d[node] = d.get(node, { 'in': [], 'out': [] })
            for out in outs:
                d[node]['out'].append(out)
                d[out] = d.get(out, { 'in': [], 'out': [] })
                d[out]['in'].append(node)
        paths = []
        verticesUsed = set()  # this set is used to find isolated cycles
        # Find standard non-branching paths
        for vertex,inout in d.items():
            if len(inout['in']) != 1 or len(inout['out']) != 1:
                for outgoingEdge in inout['out']:
                    nonBranchingPath = [vertex, outgoingEdge]
                    verticesUsed.add(vertex)
                    verticesUsed.add(outgoingEdge)
                    while len(d[outgoingEdge]['in']) == 1 and len(d[outgoingEdge]['out']) == 1:
                        outgoingEdge = d[outgoingEdge]['out'][0]
                        nonBranchingPath.append(outgoingEdge)
                        verticesUsed.add(outgoingEdge)
                    paths.append(nonBranchingPath)
        return paths
    
        #############################################################################################
        
    Trie = ModifiedSuffixTrieConstruction(Text)
    Paths = MaximalNonBranchingPaths(Trie)
    for path in Paths:
        start=Trie[path[0]][path[1]][1]
        lenpath=len(path)-1
        for n in range(lenpath):
            if type(Trie[path[n]])==int:
                break
            elif len(Trie[path[n]])==1:
                del Trie[path[n]]
            else:
                del Trie[path[n]][path[n+1]]
        Trie[path[0]].update({path[-1]:[start,lenpath]})
    return Trie


# In[1]:

def LongestRepeat(Text):
    
    '''
    Find the longest repeat in a string
    INPUT:
    Text='TTATATCGTTTTATCGTT$' # a string with an extra symbol '$' at the end
    OUTPUT='TATCGTT'
    '''
    
    from Bioinf_functions import ModifiedSuffixTreeConstruction
    
    graph=ModifiedSuffixTreeConstruction(Text)

    ### trim tree to keep only inner nodes ###        
    leaves=[k for k,v in graph.items() if type(v) == int]
    for k,v in graph.items():
        if k in leaves:
            del graph[k]
        else:
            for node in v.keys():
                if node in leaves:
                    del graph[k][node]

    ### get max path ###
    sets=[]
    for k,v in graph.items():
        if type(v)==dict:
            sets+=[[k,[out for out in v.keys()]]]

    d = {} 
    for node,outs in sets:
        d[node] = d.get(node, { 'in': [], 'out': [] })
        for out in outs:
            d[node]['out'].append(out)
            d[out] = d.get(out, { 'in': [], 'out': [] })
            d[out]['in'].append(node)

    maxreplen=0
    leaves=[k for k,v in graph.items() if len(v)==0]
    backtrack={}
    while len(leaves)>0:
        nodeend=leaves.pop(-1)
        nodeout=nodeend
        replen=0
        path=[nodeout]
        backtrack[str(path)]=''
        while True:
            nodein=d[nodeout]['in'][0]
            replen+=graph[nodein][nodeout][1]
            repstart=graph[nodein][nodeout][0]
            path+=[nodein]
            backtrack[str(path)]=Text[repstart:repstart+graph[nodein][nodeout][1]]+backtrack[str(path[:-1])]
            if nodein==0:
                break
            else:
                nodeout=nodein
        if replen>maxreplen:
            maxreplen=replen
            repstart=graph[nodein][nodeout][0]
            maxrep=Text[repstart:repstart+maxreplen]
            goodpath=path
    return backtrack[str(goodpath)]


# In[ ]:



