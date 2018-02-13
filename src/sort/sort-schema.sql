-- sort prototype schema

-- DROP VIEWS

drop view if exists saved_assignments_with_unk;
drop view if exists saved_assignments;
drop view if exists saved_variants_with_kits;
drop view if exists saved_variants;
drop view if exists saved_kits;

drop view if exists perfect_assignments_with_unk;
drop view if exists perfect_assignments;
drop view if exists perfect_variants_with_kits;
drop view if exists perfect_variants;
drop view if exists kits_view;

-- DROP TABLES

drop table if exists s_variants;
drop table if exists s_calls;
drop table if exists s_kits;

drop table if exists s_mx_kits;
drop table if exists s_mx_idxs;
drop table if exists s_mx_variants;
drop table if exists s_mx_calls;

-- CREATE TABLES

create table s_variants (
 variant_id int,  
 variant_loc varchar(10),
 name varchar(20)
);

create table s_calls(
 kit_id varchar(10), 
 variant_loc varchar(10),
 assigned boolean 
);

create table s_kits(
 kit_id  varchar(10)
);

create table s_mx_kits(
 kit_id  varchar(10),
 idx int
);

create table s_mx_variants (
 variant_id int,  
 ref_variant_id int,
 variant_loc varchar(10),
 name varchar(20), 
 idx int
);

create table s_mx_idxs(
 type_id int,           -- 0 = variants, 1 = kits
 axis_id varchar(10),   -- either the variant id or the kit_id 
 idx int
);

create table s_mx_calls (
 kit_id varchar(10),        
 variant_loc varchar(10),
 assigned boolean,
 confidence int,
 changed int
);

-- CREATE VIEWS

create view kits_view AS 
  SELECT DISTINCT kit_id from s_calls;

create view perfect_variants AS
  SELECT DISTINCT V.name, V.variant_loc, V.variant_id
  FROM s_calls C, s_variants V
  WHERE (C.assigned = -1 OR V.name = 'top') AND
  V.variant_loc = C.variant_loc;

create view perfect_variants_with_kits AS
  SELECT K.kit_id, PV.name, PV.variant_loc, PV.variant_id
  FROM perfect_variants PV
  CROSS JOIN kits_view K;

create view perfect_assignments AS
  SELECT C.kit_id, PV.name, PV.variant_loc, PV.variant_id, C.assigned
  FROM s_calls C, perfect_variants PV
  WHERE C.variant_loc = PV.variant_loc;

create view perfect_assignments_with_unk AS
  SELECT PVK.kit_id, PVK.name, ifnull(PVKA.assigned,0) as assigned, PVK.variant_loc, PVK.variant_id
  FROM perfect_variants_with_kits PVK
  LEFT JOIN perfect_assignments PVKA
  ON PVK.variant_loc = PVKA.variant_loc AND
  PVK.kit_id = PVKA.kit_id;

create view saved_kits AS 
  SELECT DISTINCT kit_id from s_mx_kits;

create view saved_variants AS 
  SELECT DISTINCT name, variant_loc, variant_id
  FROM s_mx_variants;

create view saved_variants_with_kits AS
  SELECT SK.kit_id, SV.name, SV.variant_loc, SV.variant_id
  FROM saved_variants SV
  CROSS JOIN saved_kits SK;

create view saved_assignments AS
  SELECT SC.kit_id, SV.name, SV.variant_loc, SV.variant_id, SC.assigned
  FROM s_mx_calls SC, saved_variants SV
  WHERE SC.variant_loc = SV.variant_loc;

create view saved_assignments_with_unk AS
  SELECT SVK.kit_id, SVK.name, ifnull(SVKA.assigned,0) as assigned, SVK.variant_loc, SVK.variant_id, VI.idx, KI.idx
  FROM s_mx_idxs VI, s_mx_idxs KI, saved_variants_with_kits SVK
  LEFT JOIN saved_assignments SVKA
  ON SVK.variant_loc = SVKA.variant_loc AND
  SVK.kit_id = SVKA.kit_id
  WHERE VI.type_id = 0 and VI.axis_id = SVK.name AND 
  KI.type_id = 1 and KI.axis_id = SVK.kit_id
  ORDER BY 6,7;

