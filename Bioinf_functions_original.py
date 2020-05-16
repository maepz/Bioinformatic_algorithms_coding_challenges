
# coding: utf-8

# In[1]:

'''PatternCount(Text, Pattern)

Count the frequency on a pattern(word) in a text(sequence)'''

def PatternCount(Text, Pattern):
    count=0
    for n in range(0,(len(Text)-len(Pattern)+1)):
                if Text[n:n+len(Pattern)]== Pattern:
                    count += 1
    print(count)                
    return              
#PatternCount('AACTGGTCACACATC','AC')
#PatternCount('CGCGATACGTTACATACATGATAGACCGCGCGCGATCATATCGCGATTATC', 'CGCG')


# In[3]:

'''FrequentWords(Text, k)

Find all most frquent words of size k in a sequence'''

def FrequentWords(Text, k):
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


# In[52]:

'''FindingFrequentWordsBySorting(Text , k)'''

#------------------------------------------------
def PatternToNumber(Pattern):
    d={'A':0,'C':1,'G':2,'T':3}
    if Pattern == '':
        return 0
    symbol = Pattern[-1]
    Pattern= Pattern[:-1]
    return(4 * PatternToNumber(Pattern) + d[symbol])
#--

def NumberToPattern(Number, k):
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
#--

def FindingFrequentWordsBySorting(Text , k):
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
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(sequence))
    return(reverse_complement);


#ReverseComplement('ACGTTGCA')
#ReverseComplement('TTGTGTC')


# In[49]:

'''PatternMatch(Text,Pattern)
Finds all the positions of a pattern in a text'''

def PatternMatch(Text,Pattern):
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


# In[1]:

'''PatternToNumber(Pattern)
Transform a pattern of DNA letters into an index number;
Convert k-mers into the 4^k different integers between 0 and 4^k − 1'''

#V1
def PatternToNumber1(Pattern):
    d={'A':0,'C':1,'G':2,'T':3}
    number=0
    for i in range(0,len(Pattern)):
        number += d[Pattern[i]]*(4**(len(Pattern)-1-i))
    return(number);
#--------------------------------------
# V2
def PatternToNumber(Pattern):
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


# In[17]:

'''FrequencyArray(Text , k)'''

def PatternToNumber(Pattern):
    d={'A':0,'C':1,'G':2,'T':3}
    number=0
    for i in range(0,len(Pattern)):
        number = number + d[Pattern[i]]*(4**(len(Pattern)-1-i))   
    return(number);
#---------------------------
def FrequencyArray(Text,k):
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

def HammingDistance(p,q):
    HamDi=int(0)
    for i in range(len(p)):
        if p[i]==q[i]:
            continue
        else:
            HamDi+=+1
    return(HamDi)
#--
def ApproximatePatternCount(Text,Pattern,d):
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


#--------------------------------------------------
def HammingDistance(p,q):
    HamDi=int(0)
    for i in range(len(p)):
        if p[i]==q[i]:
            continue
        else:
            HamDi+=+1
    return(HamDi)
#--
def Neighbors(Pattern, d):
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
def HammingDistance(p,q):
    HamDi=int(0)
    for i in range(len(p)):
        if p[i]==q[i]:
            continue
        else:
            HamDi+=+1
    return(HamDi)
#--
def Neighbors(Pattern, d):
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
#--
def FrequentWordsWithMismatches(Text, k, d):
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

def HammingDistance(p,q):
    HamDi=int(0)
    for i in range(len(p)):
        if p[i]==q[i]:
            continue
        else:
            HamDi+=+1
    return(HamDi)
#--
def ApproximatePatternMatch(Text,Pattern,d):
    posi=list()
    for n in range(0,(len(Text)-(len(Pattern)-1))):
        word=Text[n:n+len(Pattern)]
        if HammingDistance(Pattern,word)<=d:
            posi.append(n)
    return(posi)


# In[7]:

'''MOTIFENUMERATION(Dna, k, d)
A brute force algorithm for motif finding'''

def HammingDistance(p,q):
    HamDi=int(0)
    for i in range(len(p)):
        if p[i]==q[i]:
            continue
        else:
            HamDi+=+1
    return(HamDi)
#--
def Neighbors(Pattern, d):
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

#--

def MOTIFENUMERATION(Dna, k, d):
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
def HammingDistance(p,q):
   HamDi=int(0)
   for i in range(len(p)):
       if p[i]==q[i]:
           continue
       else:
           HamDi+=+1
   return(HamDi)

#--
def d(Pattern, Dna):
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
#--
def HammingDistance(p,q):
    HamDi=int(0)
    for i in range(len(p)):
        if p[i]==q[i]:
            continue
        else:
            HamDi+=+1
    return(HamDi)
#--
def NumberToPattern(Number, k):
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

#--
def d(Pattern, Dna):
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

#--
def MEDIANSTRING(Dna, k):
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
from numpy import *

#--

def PatternToNumber(seq):
    d={'A':0,'C':1,'G':2,'T':3}
    number=0
    revseq=seq[::-1]
    for i in range(len(revseq)):
        number += d[revseq[i]]*(4**i)
    return(number)

#--


def Profile_most(Text,k,Profile):
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
from numpy import *
import numpy as np
def MakeProfileRaw(Motif_list):
    k=len(Motif_list[0])
    matrix_occ = np.zeros((4,k))
    for i in range(len(Motif_list)):
        for n in range(len(Motif_list[i])):
            matrix_occ[(PatternToNumber(Motif_list[i][n])),n]+=1   
    Profile_matrix= matrix_occ / matrix_occ.sum(axis=0, dtype='float')       
    return Profile_matrix  

##V2: With pseudo counts, get rid of 0 probability

def MakeProfile(Motif_list):
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
import numpy as np
def Consensus(Profile):
    maxarg= list(Profile.argmax(0))
    Cons=''
    for posi in maxarg:
        Cons+=NumberToPattern(int(posi), 1)
    return Cons


# In[16]:

def HammingDistance(p,q):
    HamDi=int(0)
    for i in range(len(p)):
        if p[i]==q[i]:
            continue
        else:
            HamDi+=+1
    return(HamDi)

#--
def NumberToPattern(Number, k):
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
#--
def PatternToNumber(Pattern):
    d={'A':0,'C':1,'G':2,'T':3}
    if Pattern == '':
        return 0
    symbol = Pattern[-1]
    Pattern= Pattern[:-1]
    return(4 * PatternToNumber(Pattern) + d[symbol])

def score(motifs):
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
    score=1
    for n in range(len(kmer)):
        score *= Profile[n][kmer[n]]
    return score

def MakeProfile(motifs):
    # construct profile-matrix (each column as an dict) from count-matrix
    k=len(motifs[0])
    counter = [{"A":1,"C":1,"G":1,"T":1} for x in range(k)] # start with pseudocounts
    for i in range(k):
        for j in range(len(motifs)):
            letter = motifs[j][i]
            counter[i].update({letter:1+counter[i][letter]})     
    for d in counter:
        for key in d: d[key] /= (len(motifs)+4) # +4 due to pseudocounts
 
    return counter 
#--
def Profile_most(Text,k,Profile):
    max_score=0
    most_probable_kmer=Text[0:0+k]
    for n in range(0,len(Text)-(k-1)):
        kmer=Text[n:n+k]
        score=1
        for i in range(0,len(kmer)):
            probaN=Profile[PatternToNumber(kmer[i])][i]
            score *= probaN
        if score>max_score:
            max_score=score
            most_probable_kmer=kmer
        else:
            continue
    return(most_probable_kmer)
#--
def Consensus(Profile):
    maxarg= list(Profile.argmax(0))
    Cons=''
    for posi in maxarg:
        Cons+=NumberToPattern(int(posi), 1)
    return Cons
#--
def d(Pattern, Dna):
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
#--
def GREEDYMOTIFSEARCH(Dna,k,t):
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


def Pr(kmer,Profile):
    score=1
    for i in range(0,len(kmer)):
        probaN=Profile[PatternToNumber(kmer[i])][i]
        score *= probaN
    return score

def MakeProfile(Motif_list):
    k=len(Motif_list[0])
    matrix_occ = np.zeros((4,k))
    matrix_occ +=1
    for i in range(len(Motif_list)):
        for n in range(len(Motif_list[i])):
            matrix_occ[(PatternToNumber(Motif_list[i][n])),n]+=1   
    Profile_matrix= matrix_occ / matrix_occ.sum(axis=0, dtype='float')
    #return matrix_occ.sum(axis=0, dtype='float')
    return Profile_matrix  

def PatternToNumber(Pattern):
    d={'A':0,'C':1,'G':2,'T':3}
    if Pattern == '':
        return 0
    symbol = Pattern[-1]
    Pattern= Pattern[:-1]
    return(4 * PatternToNumber(Pattern) + d[symbol])

def Random(Motif_list,Profile):
    Array=list()
    for motif in Motif_list:
        Array.append(Pr(motif,Profile))
    PArray=[v/sum(Array) for v in Array]
    i = np.random.choice( len(PArray), p=PArray)
    return i


# In[19]:

def score(motifs):
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
    score=1
    for i in range(0,len(kmer)):
        probaN=Profile[PatternToNumber(kmer[i])][i]
        score *= probaN
    return score

def MakeProfile(Motif_list):
    k=len(Motif_list[0])
    matrix_occ = np.zeros((4,k))
    matrix_occ +=1
    for i in range(len(Motif_list)):
        for n in range(len(Motif_list[i])):
            matrix_occ[(PatternToNumber(Motif_list[i][n])),n]+=1   
    Profile_matrix= matrix_occ / matrix_occ.sum(axis=0, dtype='float')
    #return matrix_occ.sum(axis=0, dtype='float')
    return Profile_matrix  

def PatternToNumber(Pattern):
    d={'A':0,'C':1,'G':2,'T':3}
    if Pattern == '':
        return 0
    symbol = Pattern[-1]
    Pattern= Pattern[:-1]
    return(4 * PatternToNumber(Pattern) + d[symbol])

def Random(Motif_list,Profile):
    Array=list()
    for motif in Motif_list:
        Array.append(Pr(motif,Profile))
    PArray=[v/sum(Array) for v in Array]
    i = np.random.choice( len(PArray), p=PArray)
    return i

def GIBBSSAMPLER(Dna, k, t, N):

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
            #print BestMotifs, sM,BestScore, x
            x+=1
        else: 
            x+=1
            #print BestMotifs, score(BestMotifs),x
    return BestMotifs, score(BestMotifs)


# In[ ]:



