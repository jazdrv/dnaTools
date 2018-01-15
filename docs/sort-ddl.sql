-- DDL

create table variants (
 -- variant_id int,
 variant_id int, -- this should be location, putting it as variant_id to keep things simple
 name text --,
 -- reference varchar(2), -- commenting out right now cuz not part of ian's doc
);

create table calls (
 call_id int,
 hg38 int,
 name text,
 variant_id int
)

create table variants (
 variant_id int,
 location int,
 passed boolean 
);
 
create table sort_variants (
 variant_id int
 sort_order int
);

create table sort_kits(
 kit_id  int,  -- later this can be person_id
 sort_order int
);

-- hide this for now {{{

-- create table sort_problems (
--  kit id int
--  variants text
--  issue_type tinyint
--  notes text
--  resolution text
-- );

-- create table call_overrides (
--   call_id int
--   old
-- );

-- }}}

-- DML

-- insert into variants (variant_id,reference) values (3019783,'M343');
-- insert into variants (variant_id,refernece) values (6920349

