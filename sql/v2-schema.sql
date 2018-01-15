/* data sets, aka kits, aka people */

create table files(
    ID INTEGER PRIMARY KEY,
    kitname TEXT,              -- surname or other name associated with file
    dnaID TEXT,                -- the "kit ID" from FTDNA or other testing lab
    seq INTEGER DEFAULT 0,     -- number can be used to order files arbitrarily
    oID INTEGER REFERENCES origins(ID),  -- where data came from
    buildID INTEGER REFERENCES builds(ID)
    );

create index fileidx on files(kitname);

/* description of where data comes from; e.g. Big-Y, FGC Elite, pseudo */

create table origins(
    ID INTEGER PRIMARY KEY,
    originname TEXT,
    description TEXT
    );

/* the ranges (min,max) as reported in BED files */

create table bedranges(
    ID INTEGER PRIMARY KEY,
    minaddr INTEGER,
    maxaddr INTEGER
    );

create index rangeidx on bedranges(minaddr);
create index rangeidx2 on bedranges(minaddr);

/* ranges covered by tests as reported in the individual's BED file */

create table bed(
    pID INTEGER REFERENCES files(ID),
    bID INTEGER REFERENCES bedranges(ID)
    );

/* calls reported by the individual's VCF file */

create table vcfcalls(
    pID INTEGER REFERENCES files(ID),
    vID INTEGER REFERENCES variants(ID)
    );

create index vcfidx on vcfcalls(vID);

/* per-kit call statistics */

create table vcfstats(
    pID INTEGER REFERENCES files(ID),
    ny INTEGER,
    nv INTEGER,
    ns INTEGER,
    ni INTEGER,
    nr INTEGER
);

/* per-kit BED coverage statistics */

create table bedstats(
    pID INTEGER REFERENCES files(ID),
    coverage1 INTEGER,
    coverage2 INTEGER,
    nranges INTEGER
);

/* per-kit calls that are classified REJECTs */

create table vcfrej(
    pID INTEGER REFERENCES files(ID),
    vID INTEGER REFERENCES variants(ID)
    );

create index vcfridx on vcfrej(vid);

/* info about this entire data set or run */

create table meta(
    descr TEXT,
    val TEXT
    );

/* list of all variants known for this computation */
/* ref,alt: normalize? maybe not, as this table's probably not huge */

create table variants(
    ID INTEGER PRIMARY KEY,
    buildID INTEGER references builds(ID),
    pos INTEGER,
    ref TEXT,
    alt TEXT,
    UNIQUE(buildID, pos, ref, alt)
    );

create index varidx on variants(buildID, pos);

/* names and aliases associated with variants */

create table snpnames(
    vID INTEGER REFERENCES variants(ID),
    snpname TEXT
    );

/* builds (reference genome assembly) associated with data sets */

create table builds(
    ID INTEGER PRIMARY KEY,
    buildname TEXT
    );

