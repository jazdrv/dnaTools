-- sort prototype schema

-- DROP VIEWS {{{

DROP VIEW IF EXISTS v_imx_kits;
DROP VIEW IF EXISTS v_imx_variants;
DROP VIEW IF EXISTS v_imx_variants_test1;
DROP VIEW IF EXISTS v_imx_variants_test2;
DROP VIEW IF EXISTS v_imx_variants_all;
DROP VIEW IF EXISTS v_imx_variants_with_kits;
DROP VIEW IF EXISTS v_imx_assignments;
DROP VIEW IF EXISTS v_imx_assignments_with_unk;
DROP VIEW IF EXISTS v_unq_pos_variants;
DROP VIEW IF EXISTS v_unq_pos_variants_to_kits;
-- DROP VIEW IF EXISTS v_imx_variants_pos;
DROP VIEW IF EXISTS v_imx_variants_with_kits;
DROP VIEW IF EXISTS v_all_calls_with_kit;
DROP VIEW IF EXISTS v_all_calls_with_kits;

DROP VIEW IF EXISTS v_mx_kits;
DROP VIEW IF EXISTS v_mx_variants;

DROP VIEW IF EXISTS v_imx_variants_base;
DROP VIEW IF EXISTS v_imx_variants_lim;

DROP VIEW IF EXISTS v_pos_call_chk;
DROP VIEW IF EXISTS v_neg_call_chk1;
DROP VIEW IF EXISTS v_neg_call_chk2;
DROP VIEW IF EXISTS v_pos_neg_call_chk;
DROP VIEW IF EXISTS v_max_snpnames;
DROP VIEW IF EXISTS v_pos_call_chk_with_kits;
DROP VIEW IF EXISTS v_unk_call_chk_with_kits;

DROP VIEW IF EXISTS v_only_pos_variants;
DROP VIEW IF EXISTS v_only_neg_variants;

DROP VIEW IF EXISTS v_ref_variants;
DROP VIEW IF EXISTS v_ref_variants_hg19;

-- }}}
-- DROP TABLES {{{

DROP TABLE IF EXISTS mx_kits;
DROP TABLE IF EXISTS mx_variants;
DROP TABLE IF EXISTS mx_idxs;
DROP TABLE IF EXISTS mx_calls;
DROP TABLE IF EXISTS mx_dupe_variants;
DROP TABLE IF EXISTS mx_notused_variants;
DROP TABLE IF EXISTS mx_sort_recommendations;
DROP TABLE IF EXISTS mx_sups_subs;
DROP TABLE IF EXISTS mx_variant_stash;
DROP TABLE IF EXISTS mx_clade_priorities;
DROP TABLE IF EXISTS tmp1;
DROP TABLE IF EXISTS tmp2;
DROP TABLE IF EXISTS mx_call_negs;
DROP TABLE IF EXISTS mx_assignments_with_unk;

-- }}}
-- DROP INDEXES {{{

DROP INDEX IF EXISTS snpidx;
DROP INDEX IF EXISTS vcfcallsidx;
DROP INDEX IF EXISTS tmp1idx1;
DROP INDEX IF EXISTS tmp1idx2;

-- }}}
-- CREATE TABLES {{{

create table tmp1 (
    pid int,
    vid int,
    pos int,
    val bool
);

CREATE TABLE tmp2 (
    vid int,
    pid int,
    pidvid text
);

CREATE TABLE mx_kits(
 ID int,
 kitId text
);

CREATE TABLE mx_variants (
 ID int,  
 ref_variant_id int,
 name text,
 pos int,
 consistencyFlg text
);

CREATE TABLE mx_notused_variants (
 vId int,  
 reasonId int -- (other than unk's), 1 = only pos's, -1 = only neg's 
);

CREATE TABLE mx_dupe_variants (
 vId int,  
 dupe_vID int
);

CREATE TABLE mx_idxs(
 type_id int,   -- 0 = variants, 1 = kits (people)
 axis_id int,   -- either the variant id (vID) or the kit_id (pID) 
 idx int
);

CREATE TABLE mx_calls (
 pID int,
 vID int,
 assigned boolean,
 confidence int, -- 1 is ambiguous, 2 is solid
 changed int,    -- original val 
 removal int    -- removal reason: 1 no neg, 2 no pos, 3 dupe
);

CREATE TABLE mx_sort_recommendations (
 vID int,
 instructions text
);

CREATE TABLE mx_sups_subs (
 sup int,
 sub int
);

CREATE TABLE mx_variant_stash (
 ID int
);

create table mx_clade_priorities(
  ID INTEGER PRIMARY KEY, --  autoincrement,
  snpname text,
  vID int
);

create table mx_call_negs (
  vid1 int,
  vid2 int,
  name2 text, 
  pos int,
  pid int,  
  assigned int,
  genotype text
);

-- create table mx_assignments_with_unk (
--   pID int, 
--   name text, 
--   assigned int,
--   pos int,
--   vID int,
--   kitID text,
--   val bool,
--   genotype text
-- );

-- }}}
-- CREATE VIEWS (import to matrix) {{{

--jct CREATE VIEW v_max_snpnames AS
--jct   SELECT max(snpname) as snpname,vID from snpnames group by 2;

--jct CREATE VIEW v_imx_kits AS 
--jct   SELECT DISTINCT C.pID, max(D.kitId) as kitId from vcfcalls C 
--jct   INNER JOIN dataset D 
--jct   ON C.pID = D.ID 
--jct   GROUP BY 1;

--jct CREATE VIEW v_pos_call_chk as 
--jct   SELECT DISTINCT C.vID, V.pos
--jct   FROM vcfcalls C 
--jct   INNER JOIN variants V
--jct   ON C.vID  = V.ID and C.assigned = 1 AND C.genotype = '1/1';

--jct CREATE VIEW v_pos_call_chk_with_kits as 
--jct   SELECT DISTINCT C.vID, V.pos, C.pID, cast(C.pID as varchar(15))||'|'||cast(C.vID as varchar(15)) as pidvid
--jct   FROM vcfcalls C 
--jct   INNER JOIN variants V 
--jct   ON
--jct   C.vID  = V.ID and C.assigned = 1 AND C.genotype = '1/1';

--jct CREATE VIEW v_unk_call_chk_with_kits as 
--jct   SELECT DISTINCT C.vID, V.pos, C.pID, cast(C.pID as varchar(15))||'|'||cast(C.vID as varchar(15)) as pidvid
--jct   FROM vcfcalls C 
--jct   INNER JOIN variants V
--jct   ON C.vID  = V.ID and C.assigned = -1;

--jct CREATE VIEW v_all_calls_with_kits as 
--jct   SELECT DISTINCT C.vID, V.pos, C.pID, cast(C.pID as varchar(15))||'|'||cast(C.vID as varchar(15)) as pidvid
--jct   FROM vcfcalls C 
--jct   INNER JOIN variants V
--jct   ON C.vID  = V.ID; 

--jct CREATE VIEW v_neg_call_chk1 as 
--jct   SELECT DISTINCT C.vID
--jct   FROM vcfcalls C
--jct   WHERE C.assigned = 1 AND C.genotype = '0/0'
--jct   union 
--jct   select distinct vid2 as vID from mx_call_negs
--jct   ;

--jct CREATE VIEW v_pos_neg_call_chk as 
--jct   SELECT DISTINCT vID FROM v_pos_call_chk UNION SELECT vID FROM v_neg_call_chk1;

--jct CREATE VIEW v_neg_call_chk2 as 
--jct   SELECT DISTINCT vID FROM v_neg_call_chk1 UNION SELECT DISTINCT vID FROM tmp2;

--jct CREATE VIEW v_imx_variants AS
--jct   SELECT DISTINCT ifnull(S.snpname,V.ID) as name, V.pos, V.ID 
--jct   FROM vcfcalls C, v_pos_call_chk P, v_neg_call_chk2 N, variants V
--jct   LEFT JOIN v_max_snpnames S
--jct   ON S.vID = V.ID
--jct   WHERE N.vID = V.ID AND P.vID = V.ID AND V.ID = C.vID;

-- CREATE VIEW v_imx_variants_pos AS
--   SELECT DISTINCT ifnull(S.snpname,V.ID) as name, V.pos, V.ID 
--   FROM vcfcalls C, v_pos_call_chk P, variants V
--   LEFT JOIN v_max_snpnames S
--   ON S.vID = V.ID
--   WHERE P.vID = V.ID AND V.ID = C.vID;

--jct CREATE VIEW v_only_pos_variants AS
--jct   SELECT DISTINCT P.vID from v_pos_call_chk P 
--jct   LEFT JOIN v_neg_call_chk2 N
--jct   ON P.vID = N.vID
--jct   WHERE N.vID is Null;
  
--jct CREATE VIEW v_only_neg_variants AS
--jct   SELECT DISTINCT N.vID 
--jct   FROM v_neg_call_chk2 N
--jct   LEFT JOIN v_pos_call_chk P
--jct   ON P.vID = N.vID
--jct   WHERE P.vID is Null;
  
--jct CREATE VIEW v_imx_variants_with_kits AS
--jct   SELECT DISTINCT K.pID, PV.name, PV.pos, PV.ID as vID, K.kitId
--jct   FROM v_imx_variants PV
--jct   CROSS JOIN v_imx_kits K;

--jct CREATE VIEW v_imx_assignments AS
--jct   SELECT DISTINCT C.pID, PV.name, PV.pos, PV.ID as vID, C.assigned,C.genotype
--jct   FROM vcfcalls C 
--jct   INNER JOIN v_imx_variants PV
--jct   ON C.vID = PV.ID;

--jct CREATE VIEW v_imx_assignments_with_unk AS
--jct   SELECT DISTINCT PVK.pID, PVK.name, ifnull(PVKA.assigned,0) as assigned, PVK.pos, PVK.vID, PVK.kitId, T.val, PVKA.genotype, CN.assigned, CN.genotype 
--jct   FROM v_imx_variants_with_kits PVK, tmp1 T
--jct   LEFT JOIN v_imx_assignments PVKA
--jct   ON PVK.vID = PVKA.vID AND PVK.pID = PVKA.pID 
--jct   LEFT JOIN mx_call_negs CN
--jct   ON CN.pos = PVK.pos and CN.pID = PVK.pID and CN.vid2 = PVK.vID
--jct   WHERE T.pid = PVK.pID and T.vID = PVK.vID;

--jct CREATE VIEW v_ref_variants AS
--jct   SELECT DISTINCT S.snpname, V.ID, V.pos, B.buildNm, AA.allele as
--jct   anc, DA.allele as der, ifnull(IX.idx,9999999) as idx, D1.vID as vID1, D2.vID as vID2, NU.reasonId, MXV.name
--jct   FROM build B, alleles AA, alleles DA, variants V
--jct   LEFT JOIN mx_idxs IX
--jct   ON IX.axis_id = V.ID and IX.type_id = 0
--jct   LEFT JOIN mx_variants MXV
--jct   ON MXV.ID = V.ID
--jct   LEFT JOIN snpnames S
--jct   ON S.vID = V.ID
--jct   LEFT JOIN mx_dupe_variants D1 -- #is a parent to dupe "children"
--jct   ON D1.vID = V.ID
--jct   LEFT JOIN mx_dupe_variants D2 -- #is a dupe "child" of another vix
--jct   ON D2.dupe_vID = V.ID
--jct   LEFT JOIN mx_notused_variants NU
--jct   ON NU.vID = V.ID
--jct   WHERE
--jct   V.anc = AA.ID and V.der = DA.ID
--jct   and B.buildNm = 'hg38'
--jct   and V.buildID = B.ID
--jct   ORDER BY 1;

--jct CREATE VIEW v_ref_variants_hg19 AS
--jct   SELECT DISTINCT S.snpname,V.ID,V.pos, B.buildNm
--jct   FROM build B, variants V
--jct   LEFT JOIN snpnames S
--jct   ON S.vID = v.ID
--jct   WHERE
--jct   B.buildNm = 'hg19'
--jct   and V.buildID = B.ID
--jct   ORDER BY 1;

-- }}}
-- CREATE VIEWS (inside matrix) {{{

--jct CREATE VIEW v_mx_kits AS 
--jct   SELECT DISTINCT ID,kitId FROM mx_kits;

--jct CREATE VIEW v_mx_variants AS 
--jct   SELECT DISTINCT name, pos, ID FROM mx_variants;

-- }}}
-- CREATE INDEXES{{{

CREATE INDEX snpidx on snpnames(snpname);
--jct CREATE INDEX vcfcallsidx on vcfcalls(assigned,genotype);
CREATE INDEX tmp1idx1 on tmp1(pid,vid);
CREATE INDEX tmp1idx2 on tmp1(pid,vid,val);

/*}}}*/
