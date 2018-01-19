-- schema
-- DDL for redux programs and utilities


/* data sets, aka kits, aka people */
drop table if exists datasets;
create table datasets(
    ID INTEGER PRIMARY KEY,
    kitname TEXT,              -- text name, can be easy-read alias
    seq INTEGER DEFAULT 0,     -- number can be used to order files arbitrarily
    buildID INTEGER REFERENCES builds(ID),
    -- the following fields are as presented from the DW API
    kitId TEXT,                -- previously dnaID, the ID assigned by the lab
    uploaded TEXT,             -- date the file was uploaded
    dataFile TEXT,             -- normalized file/path name
    surname TEXT,              -- ancestral surname associated with the kit
    country INTEGER references countries(ID), -- ancestral country
    normalOrig TEXT,           -- normalized ancestral origin
    long TEXT,                 -- ancestral location longitude
    lat TEXT,                  -- ancestral location latitude
    otherInfo TEXT,            -- freeform text entered by user
    lab INTEGER references labs(id),            -- name of testing lab
    origFileName TEXT,         -- the original file name from user
    birthYear INTEGER,         -- birth year of MDKA
    approxHg TEXT,             -- approximate haplogroup assigned by DW
    build INTEGER references builds(ID),        -- name of the reference build
    isNGS BOOLEAN,             -- T/F is this a NGS test result
    testType INTEGER REFERENCES testtypes(ID)   -- name of the DNA test
    );

create index fileidx on datasets(kitname);

/* description of where data comes from; e.g. Big-Y, FGC Elite, pseudo */
drop table if exists testtypes;
create table testtypes(
    ID INTEGER PRIMARY KEY,
    originname TEXT,
    description TEXT
    );

/* ancestral origins country names */
drop table if exists countries;
create table countries(
    ID INTEGER PRIMARY KEY,
    countryname TEXT           -- full text of country name
    );

/* the ranges (min,max) as reported in BED files */
drop table if exists bedranges;
create table bedranges(
    ID INTEGER PRIMARY KEY,
    minaddr INTEGER,
    maxaddr INTEGER
    );

create index rangeidx on bedranges(minaddr);
create index rangeidx2 on bedranges(minaddr);

/* ranges covered by tests as reported in the individual's BED file */
drop table if exists bed;
create table bed(
    pID INTEGER REFERENCES datasets(ID),
    bID INTEGER REFERENCES bedranges(ID)
    );

/* calls reported by the individual's VCF file */
drop table if exists vcfcalls;
create table vcfcalls(
    pID INTEGER REFERENCES datasets(ID),
    vID INTEGER REFERENCES variants(ID)
    );

create index vcfidx on vcfcalls(vID);

/* per-kit call statistics */
drop table if exists vcfstats;
create table vcfstats(
    pID INTEGER REFERENCES datasets(ID),
    ny INTEGER,
    nv INTEGER,
    ns INTEGER,
    ni INTEGER,
    nr INTEGER
);

/* per-kit BED coverage statistics */
drop table if exists bedstats;
create table bedstats(
    pID INTEGER REFERENCES datasets(ID),
    coverage1 INTEGER,
    coverage2 INTEGER,
    nranges INTEGER
);

/* per-kit calls that are classified REJECTs */
drop table if exists vcfrej;
create table vcfrej(
    pID INTEGER REFERENCES datasets(ID),
    vID INTEGER REFERENCES variants(ID)
    );

create index vcfridx on vcfrej(vid);

/* info about this entire data set or run */
drop table if exists meta;
create table meta(
    descr TEXT,
    val TEXT
    );

/* list of all variants known for this computation */
drop table if exists variants;
create table variants(
    ID INTEGER PRIMARY KEY,
    buildID INTEGER references builds(ID),
    pos INTEGER,
    ref INTEGER references alleles(ID),
    alt INTEGER references alleles(ID),
    UNIQUE(buildID, pos, ref, alt)
    );

create index varidx on variants(buildID, pos);

/* allele values, strings of DNA letters */
drop table if exists alleles;
create table alleles(
    ID INTEGER PRIMARY KEY,
    allele TEXT                   -- the allele value, e.g. "A" or "TTGT"
    );

/* names and aliases associated with variants */
drop table if exists snpnames;
create table snpnames(
    vID INTEGER REFERENCES variants(ID),
    snpname TEXT
    );

/* builds (reference genome assembly) associated with data sets */
drop table if exists builds;
create table builds(
    ID INTEGER PRIMARY KEY,
    buildname TEXT
    );

