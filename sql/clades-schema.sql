-- clades.py schema

-- DROPS

drop index if exists c_rangeidx;
drop index if exists c_rangeidx2;
drop index if exists c_fileidx;
drop index if exists c_vcfidx;
drop index if exists c_vcfidx2;
drop index if exists c_vcfridx;
drop index if exists c_varidx;

drop table if exists c_files;
drop table if exists c_bed;
drop table if exists c_vcfcalls;
drop table if exists c_vcfstats;
drop table if exists c_bedstats;
drop table if exists c_vcfrej;
drop table if exists c_meta;
drop table if exists c_variants;
drop table if exists c_snpnames;
 
-- CREATES

create table c_files(
    ID INTEGER PRIMARY KEY, 
    name TEXT, 
    seq INTEGER DEFAULT 0, 
    kit INTEGER DEFAULT 0
);
create table c_bed
    ID INTEGER REFERENCES files(ID), 
    minaddr INTEGER, 
    maxaddr INTEGER
);

-- vcfcalls mostly comes from the .vcf files, so origin is "k"
-- there are some tweaks from implications.txt, and origin is set to "^" or "<"

create table c_vcfcalls(
    ID INTEGER REFERENCES files(ID), 
    vid REFERENCES VARIANTS(id), 
    origin TEXT DEFAULT "k"
);
create table c_vcfstats(
    ID INTEGER REFERENCES files(ID), 
    ny INTEGER, 
    nv INTEGER, 
    ns INTEGER, 
    ni INTEGER, 
    nr INTEGER
);

create table c_bedstats(
    ID INTEGER REFERENCES files(ID), 
    coverage1 INTEGER, 
    coverage2 INTEGER, 
    nranges INTEGER
);
create table c_vcfrej(
    ID INTEGER REFERENCES files(ID), 
    vid INTEGER REFERENCES variants(id), 
    origin TEXT DEFAULT "k"
);
create table c_meta(
    key TEXT, 
    val TEXT
);
create table c_variants(
    ID INTEGER PRIMARY KEY, 
    pos INTEGER, 
    ref TEXT, 
    alt TEXT, 
    UNIQUE(pos,ref,alt)
);
create table c_snpnames(
    ID INTEGER REFERENCES variants(ID), 
    name TEXT
);

create index c_rangeidx on bed(id, minaddr, maxaddr);
create index c_rangeidx2 on bed(id, minaddr);
create index c_fileidx on files(name);
create index c_vcfidx on vcfcalls(id, vid);
create index c_vcfidx2 on vcfcalls(id);
create index c_vcfridx on vcfrej(id,vid);
create index c_varidx on variants(pos);

