### BACKGROUND

This project was borne out of two needs: 

(1) a desire to have a matching service for BigY tests that was better than could be provided by Family Tree DNA itself; and 

(2) a desire to understand the cultural and migratory history of R-U106, by pinning the ages and geographic distributions of clades to known archaeological or historic cultures, by providing ages and geographical distributions for each clade. 

Part (1) was originally run by David Carlisle, who developed a degenerative condition, and was taken over by Andrew Booth, who died suddenly. 

Part (2) was already being run by me. As one of the few people who could, I took on the role on Andrew's death. However, since David's software was bespoke for Macs, I needed to make a hack-and-slash attempt to reproduce the output for myself. 

+ This resulted in the old Build 37 "version 1" pipeline is on my website:
www.jb.man.ac.uk/~mcdonald/genetics.html

+ along with examples of the output. The list of SNP names can be found here:
http://ybrowse.org/gbrowse2/gff/snps_hg19.csv

+ and the repository for U106 input files is on the U106 forum:
https://groups.yahoo.com/neo/groups/R1b1c_U106-S21/files/ Big-Y files/

+ I've been helped in this by both Harald and Jef. Jef has also helped me extend the effort to cover P312, so knows the software and its shortfalls.

### SUMMARY OF PLANS

* The driving aim now is to repurpose this service for Build 38, but also to improve it, extend its reach, make it more flexible and general, and to decrease the amount of manual maintenance going into generating the haplogroup tree.
* While it would be possible to convert the original code to deal with the Build 38 files with a week or two's work, I feel this hack-and-slash effort has run its course. Part of the reasoning behind this is that too much of my time is being spent patching up the haplotrees with manual adjustments that can be automated. Another motivation is that I'm under increasing moral pressure to include other types of test into this analysis, and I think I've working out the mathematical reasoning that allows me to do so. A third motivation is that I need to learn Python for a project I'm doing at work, and SQL will be useful in another project I'm involved in, which motivates the choice of language and database software: Python + SQL, in order that I can learn both.
* At the end of this e-mail, I've attached a more detailed description of what I want to see happen. Different parts of this have different priorities, and the overall aims are fairly ambitious. One of the reasons for this is that I am trying to attract government funding to buy myself research time to work on this. This will require a project with significant scientific impact (and ideally also commercial impact) in order to get any money. This also motivates stricter ethical controls over data access and privacy.

+ The data repository we're using for this is here:
www.haplogroup-r.org/submit_data.php

+ For privacy and accountability reasons, the intention is for this to remain closed. To make progress, I've contacted some related testers in R-DF98 personally, and obtained their Build 38 data for use as test cases for use by anyone who wants to help with the coding. These include Jef and myself. They can be downloaded here:

+ https://www.dropbox.com/sh/dq9fejhcmoolkb4/AAB9nASddMtiB6B3R0AYEbAta?dl=0

+ Currently, only Family Tree DNA has converted to Build 38, but the other companies are expected to follow suit.

### PROGRESS AND DIRECTIONS

* So far, we haven't made much progress, other than me climbing the simultaneous learning curves of Python and SQL, and my hacking together some of our old code to form the Python script attached to this e-mail. Jef and I have had a bit of a think about database structures, but haven't finalised anything yet.
* I have an e-mail from Zak today, saying he has done the following:
	1. repo setup on bitbucket ... I sent you an invite
	2. sent you a way we (+ others) can "collaborate" ... irc #dnatools
* If anyone wants to do anything differently or has other suggestions, let me know. Otherwise, Zak can send out further invites. I'm in your (combined) more-experienced hands when it comes to anything to do with multi-person code writing.

### DETAILED PLAN

### PART A: Assemble the SNP data

* Goal: To capture data on people and get it into an internal database for analysis
* The records being stored on James's haplogroup-r.org datab ase contain VCF+BED zipped files with associated meta-data. A check needs to be made against this data and new files downloaded and unzipped. A data structure "Kits" references data structure "People", as some people have more than one kit at one company. Kits contains information on the kit (number, company, coverage, etc.), while people needs to contain aggregated and meta-data for each person: MDKA latitude; longitude; MDKA date; MDKA date uncertainty; MDKA text string; most-recent known haplogroup; shared testing coverage; shared testing coverage within age BED region; total testing coverage.
* VCF files need parsed for a list of variants.
* This list of variants needs compared against a database (e.g. YBrowse's snps_hg38.csv) to identify ancestral/derived values (replacing reference values).
* The list of variants needs parsed to identify and standardise identical indels which are reported differently between tests.
* Multiple kits (e.g. BigY+YElite) from the same person need identified and merged, and new VCF and BED structures created for them.
* A data matrix (Calls) is needed, encoding kits and variants. For each entry, an array is needed showing: (a) postive/negative assignment (b) positive/negative call, (c) called/uncalled, (d) read quality, (e) mapping quality, (f) read depth, (g) within the age analysis BED regions.
* The VCF files need parsed to recover (b), (d), (e) and (f) in Calls
* The BED files need parsed to recover (c) in Calls.
* Comparison to the age analysis BED file is needed to populate (g) in Calls.
* Non-NGS tests can similarly be added to Calls, including those from YSeq, the literature, and SNP packs. For YSeq, or where corresponding BED information is unavailable, no-calls must be inserted for non-positive variants.

### PART B: Assemble the tree

* Goal: To use the captured data to automatically generate a haplotree
* Create a data structure (Tree) noting clade ID; parent ID; clade name; {list of variants}; {quality scores for phylogenic accuracy}; {one or more flags for variants for later analysis}; {list of child clades}; {ancestral STR alleles}; origin latitude; origin longitude; shared coverage; shared coverage within age BED; oldest/defining MDKA; SNP age; {SNP age 95% confidence interval}; {SNP age probability distribution function [PDF] (e.g. in steps of 10 years)}; {SNP parent clade age PDF}; STR age; {STR age 95% confidence interval}; {STR age PDF}; {STR parent clade age PDF}; combined age; {combined age 95% confidence interval}; {combined age PDF},
* Select mutations with 100% calls and make a tree:
	1. Sort Calls matrix vertically by variant to get the most-common SNPs to the top.
	2. Sort Calls horizontally by person to group people with common SNPs into clades.
	3. Form clades from square blocks of positive calls {Variants,People} in this 2D space.
	4. Assign positive calls as positive and negative calls as negative in Calls (hence assignment and call are redundant, but the distinction is needed later).
	5. Identify members of these clades from People, associate Variants with a clade, and assign a parent and children to each clade.
* Progress individually through progressively worse-called SNPs...
	1. Parse a list of manual implications, which dictate where positives can be implied from comparative clades (e.g. if a kit is P310+ and L11 is uncalled, a positive is implied for L11).
	2. Parse a reject list, which manually removes bad mutations from the tree.
	3. Recognise if a mutation is clearly equivalent to an existing clade or clearly forms a new clade (e.g. if all positive variants are P312+ but some U106+ tests or low-coverage (YSeq/pack) tests are uncalled, the no calls can be assigned as negative).
	4. In ambiguous cases, identify all phylogenically possible locations in the tree, mark the variant as approximately located in the quality score array, and probabilistically chose a location among them (e.g. a variant is most likely located in a clade with more equivalent SNPs). Ambiguous singleton SNPs should be identified as such, not merged up, as they will later not be counted by the age analysis.
* Identify and merge MNPs. Also decide whether they are just poorly recorded complex indels.
* Consult a reference list of named haplogroups (e.g. a text dump of FTDNA's haplotree) to assign a name to each clade in the Tree. Otherwise, automatically assign one.
* Assign oldest or defining MDKAs to nodes in Tree from People (likely with cross-reference to a manually curated list).
* Assign most-recent known haplogroups to People.
* Compute a coverage (BED / callableLoci) for each branch in Tree, and for each tester in People:
	1. The total testing coverage for People should be the coverage of their NGS test. In cases where more than one test is taken (e.g. BigY+YElite, BigY+WGS, BigY+YSeq), the concatenation of tests should be used.
	2. The shared testing coverage for Tree nodes should be the count of the number of bases tested by two or more People under that node. The most effective way to do this is to generate new BED files for each node of the tree.
	3. A *nix-style join should be done against the age analysis BED regions to determine the coverage used for the later age analysis.
* Output a visual representation of the tree and sorted matrix of Calls.

### PART C: Importing STR matches

* Goal: To predict haplogroups for testers without NGS tests, giving probabilities, which will also provide extra data for part D.
* A data structure (STRs) is needed, containing kit number; person; and a standardised array of STR results. Nulls must be representable.
* The People array must be consulted to populate person in the STRs array.
* A list of known ancestral STRs (e.g. U106, P312, P311) can be parsed to populate said array in the Tree structure.
* The tree can be parsed from top to bottom, assigning ancestral STRs to each node. An effective way to do this seems to be to start from the ancestral STRs from the parent clade, and assign mutations from this if two-thirds* of the sub-clades share this mutation (*not sure if this should be >= or > 2/3). Exceptions to this can be managed in step 3.
* Parse list of non-NGS testers, and match each up to a clade in a probabilistic (binomial) sense (cf. SAPP) using the identified mutations in that clade.

### PART D: Geographical analysis

* Goal: To identify the population distributions of each clade, debias them to absolute numbers via correct weighting, and use the clade-averaged position to predict the population's origin.
* Working from the bottom to the top of the Tree, combine the geographic origins of People into their most-recent known haplogroups. Weighting here is important. First, each country/region needs weighted according to the depth of testing (regional weight = population at average MRCA date / testing population). Haplogroups can then be combined, working up the tree:
	1. If a clade has no sub-clades, the nodal location is = {Sum(Latitude * regional weight)/Sum(regional weight) ; Sum(Longitude * regional weight)/Sum(regional weight)}.
	2. If a clade has one or more sub-clades, the weight for each sub-clade should be multiplied by the square root of the number of testers within it, then included in the sum in the same way. (The best choice of weight may need some adjustment from a square root, but that would be a good first guess.)
	3. Update latitude and longitude of origin in Tree.
* Output a series of visualisations allowing migrations to be tracked.

### PART E: Age analysis

* Goal: To provide the age of the most-recent common ancestor for each clade
* For each Person, use Calls and Tree to count the singletons and shared SNPs which are (a) within the BED file of each Person and (b) within the regions of shared coverage for the most-recent known haplogroup.
* For each NGS test, use Poisson statistics to derive an array representing the PDF for the age of that clade, and multiply the PDF arrays together.
* If the MRCA for that node can be limited or directly set by oldest/defining MDKAs in People, or from ancient DNA, the PDFs should be limited accordingly.
* Assign the PDF to Tree and work up to the top of the tree.
* Work back down the tree, generating an array representing the PDF based on the parent of each clade. Store these separately in Tree for now.
* Repeat steps #2 to #5 for Y-STRs, using Walsh's binomial method (best for pairs of people or step-wise measures) and/or Ken Nordtvedt's inter-clade variance method (best for comparing groups of people).
* Multiply the SNP and STR PDFs together; do the same for the parent SNP and STR PDFs.
* Multiply the PDFs together with the parent PDFs to create final ages, 95% c.i.s and PDFs in Tree.
* Output the final ages as lists and visualisations.

The final output will be a time-indexed, geo-coded haplogroup tree, which can be used to track relationships and migrations over the whole of the root haplogroup's history. The aim is to apply this to R-L11 and possibly a wider set of haplogroups.

### PRIORITIES

This is a big project, so there is a need to prioritise how things happen. This prioritisation is driven largely by what the community needs, which I see as the following:

* The immediate need is to get back up and running the systems we had in place before. That means we need to start being able to generate a haplotree from input Build 38 BigY VCF/BED file pairs, and being able to output that haplotree in a meaningful form, cf. the "HTML table of clades" and "full report" in version 1. (PART A and B in the above).

In its simplest form, this means translation of the "version 1" code into Python+SQL. This is something that someone could take on immediately. The hard part is likely to be ensuring that the database architecture reflects the generalisation we want to make. This is probably also the best point to make the more basic changes to how this process operates: better automation of the haplotree generation, better treatment of different types of mutation. I have some ideas about how to improve the automation (see PART B  above), but Jef probably has more experience with finding the exceptions than I do these days. Care should be paid towards keeping the code suitably generic that we can incorporate data from other sources. My expectation here is that we will need different channels for getting data into the database - one routine for BigY, one for FGC data, one for YSeq WGS tests, etc. However, to begin with, we'd prioritise the BigY data to get the system working.

A key saving that must also be made is in the software's memory requirements (see attached "structures.pdf" - but note that these descriptions have advanced since this was written). Ideally, I want it to be able to chug through a few thousand kits on my laptop. The most memory-intensive part will be the list/matrix matching people to variants (e.g. data on 5000 people, each with 3,000,000 variant calls).

* The other parts (C, D and E) are fairly modular and have roughly equal priority in my mind, although people will expect the age analysis to be ready first.

Parts C and D will need to be written from scratch, although part C is a partially solved problem, thanks to Dave Vance's SAPP code and the numerous haplogroup predictors that are already in existence.

Part E, the age analysis, is also a solved problem, but is in need of some substantial revision: the current version does not treat the differences in test coverage very well and has no way of incorporating either ancient DNA data or paper trail data. The major change in the proposed solution is quite CPU- and disc-intensive, as it will involve repeatedly calculating the intersection of BED files. The priority here will be to get the SNP side of things working, but we'll need to leave scope for an independent age estimate calculated via STRs, from ancillary sources, and the ability to combine them to produce a final age.

-----------
(default readme below - keeping this here for ideas on how to update this file)

# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* Quick summary
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

First, clone this repo.

The main script is "src/redux.py".

This requires the following things to be installed:
* Python 3, with subpackages
    * numpy
    * pyyaml
    * requests

If your default python interpreter is not Python 3, you might want to
use virtualenv; in the homedir, do:
* mkdir venv
* mkvirtualenv -p /usr/bin/python3 venv
* . venv/bin/activate (REMEMBER THE DOT!)
* pip install numpy pyyaml requests

Then you should be able to run the redux.py script:
* cd src
* export REDUX_PATH=$CWD
* Change config.yaml use_web_api: to True, to pull down initial data
* redux.py --testdrive

This takes a while.

* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact
