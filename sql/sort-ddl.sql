-- DDL

create table variants (
 -- variant_id int, -- not needed for prototype
 variant_loc int,  -- PK
 name varchar(20),
 -- old_reference varchar(2), -- commenting out right now cuz not part of ian's doc
);

create table calls(
 call_id int, -- PK
 kit_id int, 
 variant_loc int,
 passed boolean 
};
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

create table sort_variants (
 variant_loc int
 sort_order int
);

create table sort_kits(
 kit_id  int,  -- later this can be person_id
 sort_order int
);

-- hide this for now {{{

-- create table sort_problems (
--  problem_id int -- PK
--  kit_id int, -- later can be person tbl
--  variants text,
--  issue_type tinyint,
--  notes text,
--  resolution text
-- );

-- create table call_overrides (
--   override_id -- PK  
--   call_id int -- 
--   old -- not sure how this could look
-- );

-- }}}

-- DML

-- insert into variants (variant_id,reference) values (3019783,'M343');
-- insert into variants (variant_id,refernece) values (6920349

