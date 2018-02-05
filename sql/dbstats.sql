.stats
.header on
create temporary table foo as select distinct vid from vcfcalls;
select count(*) as num_called_variants from foo;
select count(distinct v.pos) as num_positions from variants v inner join foo f on f.vid=v.id;
drop table foo;
select count(*) as num_variants from variants;
select count(*) as num_calls from vcfcalls;
select count(*) as num_alleles from alleles;
select count(*) as num_bedrange from bedranges;
select count(*) as num_bed from bed;
select count(distinct pid) as num_kits from vcfcalls;
