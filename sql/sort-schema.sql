-- sort prototype schema

-- DROPS

/* replace with schema.sql
drop table if exists s_variants;
drop table if exists s_calls;
-- drop table if exists s_call_passes;
-- drop table if exists s_call_fails;
-- drop table if exists s_sort_variants;
drop table if exists s_kits;
-- drop table if exists s_sort_kits;
*/

-- CREATES

/* replace with variants
create table s_variants (
 -- variant_id int, -- not needed for prototype
 variant_loc int,  -- PK
 name varchar(20),
 -- old_reference varchar(2), -- commenting out right now cuz not part of ian's doc
 sort_order int
);
*/

/* replace with vcfcalls
create table s_calls(
 -- call_id int, -- PK - commenting out for now
 kit_id int, 
 variant_loc int,
 assigned boolean 
);
*/
-- unique index (kit_id,variant_loc) -- commented out for now
/* (hide this for now){{{
-- note: I could see calls becoming maybe two tables.

create table call_passes (
 call_id int, -- PK
 variant_loc int -- PK
 -- new_reference varchar(2), -- not needed for prototype
 -- override int -- keep out for now 
);

create table call_fails (
 call_id int, -- PK
 variant_loc int --
 -- override int -- keep out for now   
);

}}} */

-- create table s_variants (
--  variant_loc int,
--  sort_order int
-- );

/* replace with datasets
create table s_kits(
 kit_id  int,  -- later this can be person_id
 sort_order int
);
*/

-- hide this for now {{{

-- create table s_sort_problems (
--  problem_id int -- PK
--  kit_id int, -- later can be person tbl
--  variants text,
--  issue_type tinyint,
--  notes text,
--  resolution text
-- );

-- create table s_call_overrides (
--   override_id -- PK  
--   call_id int -- 
--   old -- not sure how this could look
-- );

-- }}}


