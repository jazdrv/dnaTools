#!/bin/python

# 1. The records being stored on James's haplogroup-r.org database contain VCF+BED zipped files with associated meta-data. 
#    A check needs to be made against this data and new files downloaded and unzipped. A data structure "Kits" references data 
#    structure "People", as some people have more than one kit at one company. Kits contains information on the kit (number, 
#    company, coverage, etc.), while people needs to contain aggregated and meta-data for each person: MDKA latitude; longitude; 
#    MDKA date; MDKA date uncertainty; MDKA text string; most-recent known haplogroup; shared testing coverage; shared testing 
#    coverage within age BED region; total testing coverage.

# 2. VCF files need parsed for a list of variants.

# 3. This list of variants needs compared against a database (e.g. YBrowse's snps_hg38.csv) to identify ancestral/derived 
#    values (replacing reference values).

# 4. The list of variants needs parsed to identify and standardise identical indels which are reported differently between tests.

# 5. Multiple kits (e.g. BigY+YElite) from the same person need identified and merged, and new VCF and BED structures created for them.

# 6. A data matrix (Calls) is needed, encoding kits and variants. For each entry, an array is needed showing: (a) postive/negative 
#    assignment (b) positive/negative call, (c) called/uncalled, (d) read quality, (e) mapping quality, (f) read depth, (g) within 
#    the age analysis BED regions.

# 7. The VCF files need parsed to recover (b), (d), (e) and (f) in Calls

# 8. The BED files need parsed to recover (c) in Calls.

# 9. Comparison to the age analysis BED file is needed to populate (g) in Calls.

# 10. Non-NGS tests can similarly be added to Calls, including those from YSeq, the literature, and SNP packs. For YSeq, or 
#    where corresponding BED information is unavailable, no-calls must be inserted for non-positive variants.

print "hello"
