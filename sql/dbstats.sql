.stats
.header on
select count(*) as num_called_variants from (select v.id,v.anc,v.der,count(c.pid) from variants v inner join vcfcalls c on c.vid = v.id group by 1,2,3);
select count(distinct pos) as num_positions from (select pos from variants v inner join vcfcalls c on c.vid=v.id);
select count(*) as num_variants from variants;
select count(*) as num_calls from vcfcalls;
select count(*) as num_alleles from alleles;
select count(*) as num_bedrange from bedranges;
select count(*) as num_bed from bed;
select count(distinct pid) as num_kits from vcfcalls;
