#!/usr/bin/env python3

import Bio
from Bio import SeqIO 
from Bio import pairwise2
from Bio import Align
from alignerParams import alignerInit
import numpy as N

import sys
#sys.stdout = open('trash', 'w')

inputStream=SeqIO.parse("/Users/mbucher/Sinoxolo/wuhan.fasta","fasta") 

printChunkAlignment=False

verbose=False
printSequences=False

i=0
for record in inputStream:
  if i> 1: break
  i+=1

record_wuhan=record

#inputStream=SeqIO.parse("abridged_sequences.fasta","fasta") 

def is_base(x):
  if x=='A': return(True)
  if x=='C': return(True)
  if x=='G': return(True)
  if x=='T': return(True)
  return(False)

def are_bases(x,y):
  if is_base(x) and is_base(y): 
   return(True)
  else:
   return(False)

def get_mutations(string_a, string_b):
   count=0
   for a,b in zip(string_a,string_b):
     if a!=b :
       if count > 60:
         if count < 29000:
           if are_bases(a,b):
             print(count,a,b)
     count+=1
   return(0)

def consolidate_alignments(list_one,list_two):
  list_one_bis=[]
  list_two_bis=[]
  list_one_out=[]
  list_two_out=[]
  for a in list_one:
    for b in a:
      list_one_bis.append(b)
  for a in list_two:
    for b in a:
      list_two_bis.append(b)
  current_pair_one=list_one_bis.pop(0)
  current_pair_two=list_two_bis.pop(0)
  while(1):
    next_pair_one=list_one_bis.pop(0)
    next_pair_two=list_two_bis.pop(0)
    if (( current_pair_one[1]==next_pair_one[0]) and (current_pair_two[1]==next_pair_two[0])):
       current_pair_one=( current_pair_one[0], next_pair_one[1])
       current_pair_two=( current_pair_two[0], next_pair_two[1])
       if len(list_one_bis)==0:
         list_one_out.append(current_pair_one)
         list_two_out.append(current_pair_two)
         break
    else:
       list_one_out.append(current_pair_one)
       list_two_out.append(current_pair_two)
       current_pair_one=next_pair_one
       current_pair_two=next_pair_two
    if len(list_one_bis)==0:
       list_one_out.append(current_pair_one)
       list_two_out.append(current_pair_two)
       break
  result=(list_one_out,list_two_out)
  return(result)

def myPrint(inputString):
  count=0
  chunk_size=100
  chunk_list=[inputString[i:i+chunk_size] for i in range(0,len(inputString),chunk_size)]
  for chunk in chunk_list:
    print(count)
    print(chunk)
    count+=chunk_size

def myPrintBis(inputString,countStart=0):
  inputStringA, inputStringB, inputStringC=( inputString[0], inputString[1], inputString[2])
  count=countStart
  chunk_size=100
  chunk_listA=[inputStringA[i:i+chunk_size] for i in range(0,len(inputStringA),chunk_size)]
  chunk_listB=[inputStringB[i:i+chunk_size] for i in range(0,len(inputStringB),chunk_size)]
  chunk_listC=[inputStringC[i:i+chunk_size] for i in range(0,len(inputStringC),chunk_size)]
  for chunkA,chunkB,chunkC in zip( chunk_listA, chunk_listB, chunk_listC):
   if printChunkAlignment:
    print(count)
    print(chunkA)
    print(chunkB)
    print(chunkC)
    count+=chunk_size

def get_record(value):
  inputStream=SeqIO.parse("cog_2020-09-03_sequences.fasta","fasta") 
  for i,record in enumerate(inputStream):
    if i==value:
      return(record)
  raise Exception("Record not found")

def align_record(record):
  seq1=record_wuhan
  seq2=record
  print(seq1.id)
  print(seq2.id)
  print(i)
  pad =100
  len1=len(str(seq1.seq))
  len2=len(str(seq2.seq))
  alignment_list=[]
  aligned_query_list=[]
  aligned_target_list=[]
  for i_start1 in range(0,30000,5000):
    if i_start1 != 0 : 
      i_start2=i_start1-pad
    else:
      i_start2=0
    i_end1=i_start1+5000
    i_end2=i_end1+pad
    if i_end1 > len1: i_end1=len1
    if i_end2 > len2: i_end2=len2
    seq1Bis=(str(seq1.seq)).upper()[i_start1:i_end1] 
    seq2Bis=str(seq2.seq)[i_start2:i_end2] 
    seq1Bis=seq1Bis.replace('N','X')
    seq2Bis=seq2Bis.replace('N','X')
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    alignerInit(aligner)
    alignments = aligner.align(seq1Bis, seq2Bis)
    if verbose:
      print(len(alignments))
    def myShift(pair,offset):
      a,b=pair
      a+=offset
      b+=offset
      return((a,b))
    for alignment in sorted(alignments):
        aligned_query,aligned_target=alignment.aligned
        aligned_queryShifted =[ myShift(myPair,i_start1) for myPair in aligned_query]
        aligned_targetShifted=[ myShift(myPair,i_start2) for myPair in aligned_target]
        if verbose:
          print(aligned_queryShifted)
          print(aligned_targetShifted)
        aligned_query_list.append(aligned_queryShifted)
        aligned_target_list.append(aligned_targetShifted)
        #alignment_list.append(alignment.aligned)
        #alignment_list
        #print("Score = %.1f:" % alignment.score)
        #print(alignment)
        output_form=str.split(format(alignment))[1]  
        output_formBis=str.split(format(alignment))
        #print(output_form)
        myPrintBis(output_formBis,countStart=i_start2)
        #for i in range(4): print("\n")
  if verbose:
    print(aligned_query_list)
    print(aligned_target_list)
  result=consolidate_alignments(aligned_query_list, aligned_target_list)
  print(result)
  wuhan_length=len(str(seq1.seq))
  uk_aligned=list(wuhan_length*"_")
  for segAlignQuery,seqAlignTarget in zip(aligned_query_list,aligned_target_list):
    #print(segAlignQuery)
    #print(seqAlignTarget)
    for la,lb in zip(segAlignQuery,seqAlignTarget):
      #print(la)
      #print(lb)
      #print(list(zip(la,lb)))
      q1,q2=la
      t1,t2=lb
      #print(q1)
      #print(q2)
      #print(t1)
      #print(t2)
      #print(str(seq2.seq)[t1:t2])
      for q,t in zip(range(q1,q2),range(t1,t2)):
        uk_aligned[q]=str(seq2.seq)[t]
  wuhan_string=str(seq1.seq)
  chunk_size=100
  #print("Trash")
  if printSequences:
    for index in range(0,wuhan_length,chunk_size):
      print(index)
      for c in uk_aligned[index:index+chunk_size]: 
        print(c,end='')
      print("")
      for c in wuhan_string.upper()[index:index+chunk_size]: 
        print(c,end='')
      print("")
      for a,b in zip(wuhan_string.upper()[index:index+chunk_size],uk_aligned[index:index+chunk_size]):
        if a==b:
          print(" ",end='')
        else:
          print("X",end='')
      print("")
      print("")
    print("")
  result=get_mutations(uk_aligned, wuhan_string.upper())
  #print(result)
  return  


#import pickle:
#fstream=open("firstTwoHundred.dat","wb")
#pickle.dump(output_list,fstream)
#fstream.close()

while True:
  try: 
    value=input("Enter genome number (or quit):")
    if value=='quit': break
    value=int(value)
    print(value)
    record=get_record(value)
    print(record)
    align_record(record)
  except Exception:
    print("Must be a number")
  finally:
    pass
    if verbose:
      print("Boo")
