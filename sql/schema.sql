-- schema
-- DDL for redux programs and utilities

/* enable write-ahead log (this is a persistent database setting) */
PRAGMA journal_mode=WAL;
/* 20MiB default cache size (persistent db setting) */
PRAGMA default_cache_size=20000;


/* a person who supplied DNA, or if needed, some other person */
drop table if exists person;
create table person(
    ID INTEGER PRIMARY KEY,
    firstName TEXT DEFAULT NULL,
    middleName TEXT DEFAULT NULL,
    surname TEXT DEFAULT 'Unknown',
    maidenName TEXT DEFAULT NULL,
    yHaplogroupId INTEGER DEFAULT NULL,
    mtHaplogroupId INTEGER DEFAULT NULL,
    UNIQUE(surname,firstName,middleName)
    );

/* data sets uploaded to Haplogroup-R */
drop table if exists dataset;
create table dataset(
    ID INTEGER PRIMARY KEY,
    seq INTEGER DEFAULT 0,     -- can be used to order files arbitrarily
    kitName TEXT,              -- can be a simple name for the dataset
    DNAID INTEGER references person(ID), -- whose DNA is this?
    -- the following fields are directly derived from H-R DW API json
    contact TEXT,              -- contact info for person who uploaded
    kitId TEXT,                -- the ID assigned by the lab
    surnameID INTEGER references surname(ID), -- ancestral surname
    countryID INTEGER references country(ID), -- ancestral country
    birthYr INTEGER DEFAULT NULL, -- birth year of MDKA
    labID INTEGER references lab(id), -- testing lab
    testTypeID INTEGER REFERENCES testtype(ID), -- type of DNA test
    buildID INTEGER references build(ID), -- reference build for the test
    importDt TEXT,             -- when originally imported
    fileNm TEXT,               -- normalized file/path name
    origFileNm TEXT,           -- the original file name from user
    otherInfo TEXT,            -- freeform text entered by user
    normalOrigID INTEGER references origin(ID), -- normalized ancestral origin
    lat TEXT,                  -- ancestral location latitude
    lng TEXT,                  -- ancestral location longitude
    policyVer TEXT DEFAULT NULL, -- in H-R schema
    accessToken TEXT DEFAULT NULL, -- in H-R schema
    updated TEXT,              -- date the file/record was updated
    approxHg TEXT,              -- approximate haplogroup assigned by DW
    unique(kitId)
    );

create index fileidx on dataset(kitId);

/* kits to be used for analysis */
/* this provides a way to analyze a static subset of the loaded kits */
drop table if exists analysis_kits;
create table analysis_kits(
    pID INTEGER references dataset(DNAID) -- person of interest in analysis
    );

/* kits to be excluded */
/* this provides a way to ignore some kits */
drop table if exists exclude_kits;
create table exclude_kits(
    pID INTEGER references dataset(DNAID) -- person to exclude from analysis
    );

/* variants to be excluded */
/* this provides a way to ignore some variants */
drop table if exists exclude_variants;
create table exclude_variants(
    vID INTEGER references variants(ID) -- variant to exclude from analysis
    );

/* surname associated with person and/or kit */
drop table if exists surname;
create table surname(
    ID INTEGER PRIMARY KEY,
    surname TEXT,
    UNIQUE(surname)
    );

/* description of where data comes from; e.g. Big-Y, FGC Elite, pseudo */
drop table if exists TestType;
create table testtype(
    ID INTEGER PRIMARY KEY,
    testNm TEXT,               -- the name of the test
    description TEXT,          -- further description if needed
    isNGS BOOLEAN,             -- is it a nextgen sequencing test?
    tagNm TEXT,                -- field in H-R schema
    priority INTEGER,          -- field in H-R schema
    UNIQUE(testNm,isNGS)
    );

/* regions of interest, e.g. chrY, b38, 57000000 */
drop table if exists Contig;
CREATE TABLE Contig(
    id INTEGER PRIMARY KEY,
    buildID INTEGER references build(ID), -- reference build for the region
    description TEXT,          -- e.g. chrY
    length INTEGER
    );


/* ancestral country names */
drop table if exists country;
create table country(
    ID INTEGER PRIMARY KEY,
    country TEXT,              -- full text of country name
    unique(country)
    );

/* ancestral origin normalized names */
drop table if exists origin;
create table origin(
    ID INTEGER PRIMARY KEY,
    origin TEXT,               -- normalized ancestral origin name
    unique(origin)
    );

/* testing lab names */
drop table if exists lab;
create table lab(
    ID INTEGER PRIMARY KEY,
    labNm TEXT,                -- name of the testing lab
    unique(labNm)
    );

/* the ranges (min,max) as reported in BED files */
drop table if exists bedranges;
create table bedranges(
    ID INTEGER PRIMARY KEY,
    minaddr INTEGER,
    maxaddr INTEGER
    );
create unique index rangeuniq on bedranges(minaddr, maxaddr);
--fixme these indexes may be important to performance
--create index rangeidx1 on bedranges(maxaddr);
--create index rangeidx2 on bedranges(minaddr);

/* ranges covered by tests as reported in the individual's BED file */
drop table if exists bed;
create table bed(
    pID INTEGER REFERENCES dataset(ID),
    bID INTEGER REFERENCES bedranges(ID)
    );
create index bedidx1 on bed(pID);
create index bedidx2 on bed(bID);

/* calls reported by the individual's VCF file */
drop table if exists vcfcalls;
create table vcfcalls(
    pID INTEGER REFERENCES dataset(ID),
    vID INTEGER REFERENCES variants(ID),
    callinfo INTEGER   --some info about this call packed into an int
    );
create index vcfidx1 on vcfcalls(vID);
create index vcfpid2 on vcfcalls(pID);

/* per-kit call statistics */
drop table if exists vcfstats;
create table vcfstats(
    pID INTEGER REFERENCES dataset(ID),
    ny INTEGER,
    nv INTEGER,
    ns INTEGER,
    ni INTEGER,
    nr INTEGER
);

/* per-kit BED coverage statistics */
drop table if exists bedstats;
create table bedstats(
    pID INTEGER REFERENCES dataset(ID),
    coverage1 INTEGER,
    coverage2 INTEGER,
    nranges INTEGER
);

/* info about this entire data set or run */
drop table if exists meta;
create table meta(
    descr TEXT,
    val TEXT
    );

/* table for STR definitions */
drop table if exists strs;
create table strs(
    ID INTEGER PRIMARY KEY,
    strname TEXT,
    ordering INTEGER,           -- can be used to put STRs in particular order
    unique(strname)
    );

/* STR values for testers */
drop table if exists strcalls;
create table strcalls(
    pID INTEGER REFERENCES dataset(ID), -- or maybe people?
    sID INTEGER REFERENCES strs(ID), -- STR definition
    val INTEGER -- STR marker value
    );

/* list of all variants known for this computation */
drop table if exists variants;
create table variants(
    ID INTEGER PRIMARY KEY,
    buildID INTEGER references build(ID),
    pos INTEGER,
    anc INTEGER references alleles(ID),
    der INTEGER references alleles(ID),
    UNIQUE(buildID, pos, anc, der)
    );
--indexes that are likely to help query plans in priority order
--create index varidx1 on variants(pos);
--create index varidx2 on variants(der);
--create index varidx3 on variants(anc);
--create index varidx4 on variants(buildID);

/* reference positives */
/* a listing in this table means reverse the meaning of anc and der */
drop table if exists refpos;
create table refpos(
    vID INTEGER references variants(ID),
    UNIQUE(vid)
    );

/* allele values, strings of DNA letters */
drop table if exists alleles;
create table alleles(
    ID INTEGER PRIMARY KEY,
    allele TEXT,                  -- the allele value, e.g. "A" or "TTGT"
    UNIQUE(allele)
    );

/* names and aliases associated with variants */
drop table if exists snpnames;
create table snpnames(
    vID INTEGER REFERENCES variants(ID),
    snpname TEXT,
    unique(snpname,vID)
    );
create index snpidx1 on snpnames(snpname);
create index snpidx2 on snpnames(vID);

/* build (reference genome assembly) associated with data sets */
drop table if exists build;
create table build(
    ID INTEGER PRIMARY KEY,
    buildNm TEXT,
    unique(buildNm)
    );

/* ranges from age.bed that are used in age calculations */
drop table if exists agebed;
create table agebed(
    ID INTEGER PRIMARY KEY,
    bID INTEGER REFERENCES bedranges(ID)
    );

/* tree data structure - this probably still needs work */

/*
 * For sql processing of the tree, use a closure table.
 * E.g. see https://www.slideshare.net/billkarwin/sql-antipatterns-strike-back/68-Naive_Trees_Solution_3_Closure
 * This table lists paths all nodes to all of their descendants.
 */
drop table if exists treepaths;
CREATE TABLE treepaths(
    ancestor INTEGER REFERENCES treeclade(ID),
    descendant INTEGER REFERENCES treeclade(ID),
    treedepth integer,
    UNIQUE(ancestor,descendant,treedepth)
    );
create index treeancidx on treepaths(ancestor);
create index treedesidx on treepaths(descendant);
create index treedepidx on treepaths(treedepth);

/*
 * Tree clade, aka block, part of tree definition
 * A clade represents a set of variants that a set of people (DNAIDs)
 * all have.  Additional information about a clade can go in this
 * table. E.g. we might want to store age information here.
 */
drop table if exists treeclade;
CREATE TABLE treeclade(
    ID INTEGER PRIMARY KEY,
    cladeName CHARACTER --the name we prefer, e.g. dominant SNP name
    );

/*
 * Clade variants, part of tree definition
 * This simple table links variant IDs to a clade
 */
drop table if exists cladevariants;
CREATE TABLE cladevariants(
    cladeID INTEGER REFERENCES treeclade(ID),
    vID INTEGER REFERENCES variants(ID),
    UNIQUE(cladeID,vID)
    );
create index cladevidx on cladevariants(cladeid);

/*
 * Clade kits, part of tree definition
 * This simple table links DNAIDs to a clade
 */
drop table if exists cladekits;
CREATE TABLE cladekits(
    cladeID INTEGER REFERENCES treeclade(ID),
    pID INTEGER REFERENCES person(ID),
    UNIQUE(cladeID,pID)
    );
create index cladekidx on cladekits(cladeid);

