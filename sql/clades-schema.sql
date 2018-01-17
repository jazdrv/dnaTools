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
 
-- Note: these are new ones in the v2 schema
-- builds
-- bedranges
-- origins

-- CREATES

create table c_files(
    ID INTEGER PRIMARY KEY, 
    name TEXT,  -- v2 kitname
    -- v2 kitname text
    seq INTEGER DEFAULT 0, 
    -- v2 oID int  
    kit INTEGER DEFAULT 0 -- not in v2
    -- v2 buildId int (ref builds(id))
);

-- v2 tbl
-- create table c_origins(
--    id int primary key,
--    originname text,
--    description text    
--};

-- v2 tbl
-- create table bedranges(
--    id int primary key,
--    minaddr int,
--    maxaddr int,
-- );

create table c_bed
    ID INTEGER REFERENCES c_files(ID), 
    minaddr INTEGER,  -- not in v2
    maxaddr INTEGER   -- not in v2
    -- v2 bId int (ref bedranges(ID))
);

-- vcfcalls mostly comes from the .vcf files, so origin is "k"
-- there are some tweaks from implications.txt, and origin is set to "^" or "<"

create table c_vcfcalls(
    ID INTEGER REFERENCES c_files(ID), 
    vid REFERENCES c_VARIANTS(id), 
    origin TEXT DEFAULT "k" -- not in v2
);
create table c_vcfstats(
    ID INTEGER REFERENCES c_files(ID), 
    ny INTEGER, 
    nv INTEGER, 
    ns INTEGER, 
    ni INTEGER, 
    nr INTEGER
);

create table c_bedstats(
    ID INTEGER REFERENCES c_files(ID), 
    coverage1 INTEGER, 
    coverage2 INTEGER, 
    nranges INTEGER
);
create table c_vcfrej(
    ID INTEGER REFERENCES c_files(ID), 
    vid INTEGER REFERENCES c_variants(id), 
    origin TEXT DEFAULT "k" -- not in v2
);
create table c_meta(
    key TEXT,  -- not in v2
    -- v2 descr text    
    val TEXT 
);
create table c_variants(
    ID INTEGER PRIMARY KEY, 
    -- v2 buildId (ref builds(id))
    pos INTEGER, 
    ref TEXT, 
    alt TEXT, 
    UNIQUE(pos,ref,alt) -- v2 also includes buildId
);
create table c_snpnames(
    ID INTEGER REFERENCES c_variants(ID), 
    name TEXT -- v2 snpname text
);

-- v2 tbl
-- create table c_builds (
--    id int primary key,
--    buildname text
-- );

create index c_rangeidx on c_bed(id, minaddr, maxaddr);
create index c_rangeidx2 on c_bed(id, minaddr);
create index c_fileidx on c_files(name); -- v2 kitname
create index c_vcfidx on c_vcfcalls(id, vid);
create index c_vcfidx2 on c_vcfcalls(id);
create index c_vcfridx on c_vcfrej(id,vid); -- v2 only vid
create index c_varidx on c_variants(pos); -- v2 no buildId

-- v2 index
-- create index rangeidx on c_bedranges(minaddr);
-- create index rangeidx2 on c_bedranges(minaddr);


