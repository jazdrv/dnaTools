#!/bin/python

# Pseudo-code for running haplotree comparisons

# This code does not work. It's not even written properly in Python.

# INHERITANCE FROM WORK PACKAGE A
# How this process works depends strongly on the setup acquired from Work Package A. We assume we have:
# (1) Parsed VCF calls into an Calls SQL table, Variants SQL table
# (2) Identified and called all shared variants (and perhaps singletons) among all the VCF files and, where no call exists, checked whether callable coverage in the BED file implies they are negative.
# (3) Obtained quality information from the VCF for each of these calls, such that we know which are passed and which have failed.
# (4) Obtained a reference mutation list, correcting reference positives, so we know which alleles are ancestral and which are derived.
# We will assume we are only dealing with sequencing tests for the time being, that reference calls can clearly be identified from uncalled base pairs, and that there are not large orders of magnitude difference in coverage.

# We will assume that the calls are stored with a ruleID as integer, defined as:
#	 7	Positive	Call=PASS, genotype=1/1 (or 0/0 if a ref pos)
#	 6	Mixed	Call=PASS, genotype=0/1
#	 5	Positive?	Call=FAIL, genotype=1/1 (0/0)
#	 4	Mixed?	Call=FAIL, genotype=0/1
#	 3	Negative?	Call=FAIL, genotype=0/0 (1/1)
#	 2	Negative*	No call, but covered in one of the BED file ranges, not a reference positive
#	 1	Negative	Call=PASS, genotype=0/0 (or 1/1 if a ref pos)
#	 0	Uncalled	No call, no BED file coverage
# Additional codes may be needed to encode quality information (e.g. FGC *,**,*** quality thresholds). We can then sort and process this array to form the haplotree.
# Numbers can be added to this call bitwise. We assume:
#	16 = junk;
#	32 = negative;
#	48 = positive; 
# such that VALUE%16 will give the original and int(VALUE/16) will give the overlain assignment.
# We will need a sort index and statistics for both people and variants.
# We can assume that we will need to process ~5000 people with ~30,000 variants needing investigation (out of ~12 million variants), hence ~150 million calls in total.
# Likely extrema for this analysis would be ~50000 people with ~100,000 variants needing investigation (out of ~23 million variants), hence ~5 billion calls in total.


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
	# Produces:	Array of stats for that variant:
	#	[0,1] first positive occurrence, last positive occurence
	#	[2,3] first called non-negative occurence, last -"-
	#	[4-11] count of: uncalled, negative, BED negative, failed negative, failed mixed, failed positive, mixed, positive
	#	[12-15] count of assignments: none, junk, fail, pass
	pfirst = plast = 1	# Passed first and last position
	allfirst = alllast = 1	# Positive/mixed (pass/fail) first and last
	foundFirstPos = foundFirstPPos = 0 # Have we found the first one yet?
	vcount=np.zeros(11,dtype=int) # Count of original calls (0-7)
								  # count of assignments (8-11)
	for person in pSortIdx:
		# Flag start of non-negative called values
		if (matrix[person,variant]%16>3)
			alllast=person
			if (foundFirstPos==0):
				foundFirstPos=1
				allfirst=person
		# Flag start of positive (or mixed) passed called values
		if (matrix[person,variant]%16>=6)
			plast=person
			if (foundFirstPos==0):
				foundFirstPPos=1
				pfirst=person
		vcount[matrix[variant,person]%16]++
		vcount[int(matrix[variant,person]/16)+8]++
	vstats=np.array([pfirst,plast,allfirst,alllast],dtype=int)
	vstats=np.append(vstats,vcount,axis=1)
	return vstats

def getAllPersonRanges(matrix,vSortIdx,pSortIdx):
	# Run getPersonRange over a portion of the haplotree to get range of people, returning statistics on number of passed/failed/uncalled positive/negative/mixed calls
	# Uses:		getPersonRange
	# Produces 1D table
	vstats=np.zeros((0,15),dtype=int)
	foreach variant in matrix[:,None]:
		vstats=np.append(getPersonRange(matrix,variant,PSortIdx))
	return vstats
		
def sortVariants(matrix,vSortIdx,pSortIdx,vstats):
	# Sort variants, based on quality, frequency, first person and starting position
	# Uses:		findQuality, getAllPersonRanges

def sortPeople(matrix,variant,prange,pSortIdx):
	# Re-sort a subset of people, based on positive and negative alleles
	# Uses:		getPersonRange

def isPhylogenyOK(matrix,variant,vSortIdx,pSortIdx,vstats):
	# Binary test of whether a phylogeny is ok: for all positive calls for a variant, search up the variant sort order to find if they map to the same parent
	# Uses:		findPersonRange
	# Returns:	True/false
	parentVariant=0
	isGood=.true.
	for person in pSortIdx[vstats[2]:vstats[3]]: # Loop over all called, non-negative tests
		if (matrix[variant,person]%16>3):
			testVar=vSortIdx[variant]
			while (matrix[testVar,person]%16<0):
				testVar--
		

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

def makeCSV(matrix,vSortIdx,pSortIdx,pstats,vstats):
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
	# Get variant stats
	