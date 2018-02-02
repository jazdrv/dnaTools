#!/bin/python

import numpy as np

def getVCFcalls(pos, people):
	# Extract call information about a position from the Calls table, comparing call to ancestral allele.
	# Uses:			Calls table
	# Produces:		2D array, defining alleles called for that position
	# Select variants at position <pos> from list
	# <<< variants = SELECT ID,anc,der FROM Variants WHERE buildID="hg38" AND pos=testpos
	# Assign all calls as negative (ruleID=1) until proven otherwise
	# Missing calls (ruleID=0) are handled by getBEDCalls
	vcfcalls=np.zeros[(length(people),length(variants)),dtype=int]
	# Select ancestral variant
	for var in variants[][0]:
		if (variants[var][1]==variants[var][2]) :
			ancID = variants[var][0]
			ancVar = var
	# Get calls
	for var in variants[][0]:
		# <<< varCalls = SELECT pID,ruleID FROM Calls WHERE vID = var
		for person in people:
			for vcall in varCalls:
				if (vcall == person):
					if (var == ancVar):
						vcfcalls[person]=varcalls[vcall,1]
					else:
						vcfcalls[person,var]=varcalls[vcall,1]
	vcfcalls = np.delete(vcfcalls, ancVar, axis=0)
	return vcfcalls
	
def getBEDcalls(pos, people):
	# Extract coverage information about a position from the Coverage table, comparing call to ancestral allele.
	# Uses:			Coverage table
	# Produces:		1D array, defining a coverage mask for that position
	# Get coverage information
	# <<< personIDs = SELECT pID FROM BEDCalls WHERE (pos>coverage1 AND pos<=coverage2)
	# Assign all calls to uncalled (ruleID=0)
	bedcalls=np.zeros(length(people),dtype=int)
	foreach person in people:
		foreach personID in personIDs:
			if (person == personID):
				bedcalls(person) = 1
	return bedcalls
	
def getHaplotreeRows(variants,people):
	# Append calls and coverage for a list of positions to a haplotree matrix
	# Uses:		getVCFcalls, getBEDcalls
	# Produces:	2D array giving call information on all variants at each position
	# Assumes:	If a person is positive for one variant at that position, they are negative for all others
	treerows=np.zeros((length(people),0),dtype=int)
	foreach varpos in variants:
		vcfcalls=getVCFcalls(varpos, people)
		bedcalls=getBEDcalls(bedcalls, people)
		# Mutliply arrays together, so that people who are uncalled get assigned from ruleID=1 to 0
		treerows=np.append(treerows,vcfcalls*bedcalls[None,:],axis=0)
	return treerows
		
def getPersonRange(matrix,variant,pSortIdx):
	# Find the first and last person positive for a variant in a range
	# Can also be used to check if a variant is sorted if:
	# 		last-first+1<npositive
	# Produces:	Array of [first, last, npositive]
	

def getAllPersonRanges(matrix,vSortIdx,pSortIdx):
	# Run getPersonRange over a portion of the haplotree to get range of people, returning statistics on number of passed/failed/uncalled positive/negative/mixed calls
	# Uses:		getPersonRange
	# Produces:	1D table of pRange[variant]=[first, last, 
	# nppos,npneg,npmix,nfpos,nfneg,nfmix,nnc]

def sortVariants(matrix,vSortIdx,pSortIdx):
	# Sort variants, based on quality, frequency, first person and starting position
	# Uses:		findQuality, getAllPersonRanges

def sortPeople(matrix,variant,prange,pSortIdx):
	# Re-sort a subset of people, based on positive and negative alleles
	# Uses:		getPersonRange

def isPhylogenyOK(matrix,variant,vSortIdx,pSortIdx):
	# Binary test of whether a phylogeny is ok: for all positive calls for a variant, search up the variant sort order to find if they map to the same parent
	# Uses:		findPersonRange
	# Returns:	True/false

def insertVariant(matrix,variant,vSortIdx,newpos):
	# Insert a variant into the haplotree by moving it up the calls matrix
	# Uses:		findPersonRange

def isEquivalent(matrix,variant,vSortIdx,pSortIdx):
	# Test whether a variant is equivalent to another in the calls matrix
	# Uses:		findPersonRange
	# Returns:	True/false

def editCall(matrix,variant,person,assignment):
	# Edit a call to assign it positive/negative/junk
	# Assignment (1=junk, 2=neg; 3=pos)
	# Uses:		findPersonRange

def makehaplotree:
# Pseudo-code for making the haplotree

	backMutTolPerc=0.01 # Maximum fractional and absolute numbers
	backMutTolAbs=2 # of negatives to flag as potential back mutations

	# Populate NGS-tested people from People table
	# <<< people = SELECT * FROM people WHERE (ngs=.true.)
	npeople=len(people)
	# Populate sharedVariants list from Calls table
	# <<< variants = SELECT * FROM variants WHERE (npos>1 && npos<npeople)
	callsMatrix=np.zeros((npeople,length(variants),dtype=int) # 2D matrix of calls
	vSortIdx=np.zeros(length(variants),dtype=int) # Variant sort order
	pSortIdx=np.zeros(npeople,dtype=int) # People sort order
	lastGoodVar=0 # Separates passed rows in haplotree

	# Read shared variants into the calls matrix
	callsMatrix=getHaplotreeRows(sharedVariants,ngsPeople)
	# Populate the ranges for each variant
	pstats=getAllPersonRanges(callsMatrix,vSortIdx,pSortIdx)
	# Mask out potential false positives
	for variant in sharedVariants:
		if (pstats.nppos[variant]/npeople >= 1-backMutTolPerc || npeople-pstats.nppos[variant] <= backMutTolAbs):
			# <<< editCalls(callsMatrix,variant,people[variant=negative],1)

	# Sort the variants by quality
	sortVariants(callsMatrix,vSortIdx,pSortIdx)
	# <<< lastGoodVar= last variant where everyone is positive
	# Sort through the variants, checking whether they are phylogenically ok 
	for variant in sharedVariants:
		isPhylogenyOK(callsMatrix,variant,vSortIdx)
		sortPeople(callsMatrix,variant,pSortIdx)
	# etc., etc., etc.