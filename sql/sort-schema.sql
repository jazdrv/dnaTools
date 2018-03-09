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
DROP VIEW IF EXISTS v_imx_variants_pos;
DROP VIEW IF EXISTS v_imx_variants_with_kits;
DROP VIEW IF EXISTS v_all_calls_with_kit;

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

-- }}}
-- DROP INDEXES {{{

DROP INDEX IF EXISTS snpidx;
DROP INDEX IF EXISTS vcfcallsidx;

-- }}}
-- CREATE TABLES {{{

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

-- }}}
-- CREATE VIEWS (import to matrix) {{{

CREATE VIEW v_max_snpnames AS
  SELECT max(snpname) as snpname,vID from snpnames group by 2;

CREATE VIEW v_imx_kits AS 
  SELECT DISTINCT C.pID, max(D.kitId) as kitId from vcfcalls C 
  INNER JOIN dataset D 
  ON C.pID = D.ID 
  GROUP BY 1;

CREATE VIEW v_pos_call_chk as 
  SELECT DISTINCT C.vID, V.pos
  FROM vcfcalls C 
  INNER JOIN variants V
  ON C.vID  = V.ID and C.assigned = 1 AND C.genotype = '1/1';

CREATE VIEW v_pos_call_chk_with_kits as 
  SELECT DISTINCT C.vID, V.pos, C.pID, cast(C.pID as varchar(15))||'|'||cast(C.vID as varchar(15)) as pidvid
  FROM vcfcalls C 
  INNER JOIN variants V 
  ON
  C.vID  = V.ID and C.assigned = 1 AND C.genotype = '1/1';

CREATE VIEW v_unk_call_chk_with_kits as 
  SELECT DISTINCT C.vID, V.pos, C.pID, cast(C.pID as varchar(15))||'|'||cast(C.vID as varchar(15)) as pidvid
  FROM vcfcalls C 
  INNER JOIN variants V
  ON C.vID  = V.ID and C.assigned = -1;

CREATE VIEW v_all_calls_with_kits as 
  SELECT DISTINCT C.vID, V.pos, C.pID, cast(C.pID as varchar(15))||'|'||cast(C.vID as varchar(15)) as pidvid
  FROM vcfcalls C 
  INNER JOIN variants V
  ON C.vID  = V.ID; 

CREATE VIEW v_neg_call_chk1 as 
  SELECT DISTINCT C.vID
  FROM vcfcalls C
  WHERE C.assigned = 1 AND C.genotype = '0/0';

CREATE VIEW v_pos_neg_call_chk as 
  SELECT DISTINCT vID FROM v_pos_call_chk UNION SELECT vID FROM v_neg_call_chk1;

CREATE VIEW v_neg_call_chk2 as 
  SELECT DISTINCT vID FROM v_neg_call_chk1 UNION SELECT DISTINCT vID FROM tmp2;

CREATE VIEW v_imx_variants AS
  SELECT DISTINCT ifnull(S.snpname,V.ID) as name, V.pos, V.ID 
  FROM vcfcalls C, v_pos_call_chk P, v_neg_call_chk2 N, variants V
  LEFT JOIN v_max_snpnames S
  ON S.vID = V.ID
  WHERE N.vID = V.ID AND P.vID = V.ID AND V.ID = C.vID;

CREATE VIEW v_imx_variants_pos AS
  SELECT DISTINCT ifnull(S.snpname,V.ID) as name, V.pos, V.ID 
  FROM vcfcalls C, v_pos_call_chk P, variants V
  LEFT JOIN v_max_snpnames S
  ON S.vID = V.ID
  WHERE N.vID = V.ID AND P.vID = V.ID AND V.ID = C.vID;

CREATE VIEW v_only_pos_variants AS
  SELECT DISTINCT P.vID from v_pos_call_chk P 
  LEFT JOIN v_neg_call_chk2 N
  ON P.vID = N.vID
  WHERE N.vID is Null;
  
CREATE VIEW v_only_neg_variants AS
  SELECT DISTINCT N.vID 
  FROM v_neg_call_chk2 N
  LEFT JOIN v_pos_call_chk P
  ON P.vID = N.vID
  WHERE P.vID is Null;
  
CREATE VIEW v_imx_variants_with_kits AS
  SELECT DISTINCT K.pID, PV.name, PV.pos, PV.ID as vID, K.kitId
  FROM v_imx_variants PV
  CROSS JOIN v_imx_kits K;

CREATE VIEW v_imx_assignments AS
  SELECT DISTINCT C.pID, PV.name, PV.pos, PV.ID as vID, C.assigned,C.genotype
  FROM vcfcalls C 
  INNER JOIN v_imx_variants PV
  ON C.vID = PV.ID;

CREATE VIEW v_imx_assignments_with_unk AS
  SELECT DISTINCT PVK.pID, PVK.name, ifnull(PVKA.assigned,0) as assigned, PVK.pos, PVK.vID, PVK.kitId, T.val, PVKA.genotype 
  FROM v_imx_variants_with_kits PVK, tmp1 T
  LEFT JOIN v_imx_assignments PVKA
  ON PVK.vID = PVKA.vID AND PVK.pID = PVKA.pID 
  WHERE T.pid = PVK.pID and T.vID = PVK.vID;

CREATE VIEW v_ref_variants AS
  SELECT DISTINCT S.snpname, V.ID, V.pos, B.buildNm, AA.allele as
  anc, DA.allele as der, ifnull(IX.idx,9999999) as idx, D1.vID as vID1, D2.vID as vID2, NU.reasonId, MXV.name
  FROM build B, alleles AA, alleles DA, variants V
  LEFT JOIN mx_idxs IX
  ON IX.axis_id = V.ID and IX.type_id = 0
  LEFT JOIN mx_variants MXV
  ON MXV.ID = V.ID
  LEFT JOIN snpnames S
  ON S.vID = V.ID
  LEFT JOIN mx_dupe_variants D1 -- #is a parent to dupe "children"
  ON D1.vID = V.ID
  LEFT JOIN mx_dupe_variants D2 -- #is a dupe "child" of another vix
  ON D2.dupe_vID = V.ID
  LEFT JOIN mx_notused_variants NU
  ON NU.vID = V.ID
  WHERE
  V.anc = AA.ID and V.der = DA.ID
  and B.buildNm = 'hg38'
  and V.buildID = B.ID
  ORDER BY 1;

CREATE VIEW v_ref_variants_hg19 AS
  SELECT DISTINCT S.snpname,V.ID,V.pos, B.buildNm
  FROM build B, variants V
  LEFT JOIN snpnames S
  ON S.vID = v.ID
  WHERE
  -- V.pos in (%s)
  B.buildNm = 'hg19'
  and V.buildID = B.ID
  ORDER BY 1;

-- }}}
-- CREATE VIEWS (inside matrix) {{{

CREATE VIEW v_mx_kits AS 
  SELECT DISTINCT ID,kitId FROM mx_kits;

CREATE VIEW v_mx_variants AS 
  SELECT DISTINCT name, pos, ID FROM mx_variants;

-- }}}
-- CREATE INDEXES{{{

CREATE INDEX snpidx on snpnames(snpname);
CREATE INDEX vcfcallsidx on vcfcalls(assigned,genotype);
-- CREATE INDEX tmp1idx1 on tmp1(pid,vid);
-- CREATE INDEX tmp1idx2 on tmp1(pid,vid,val);

/*}}}*/
