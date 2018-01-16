# !/bin/bash
SNP="ZZ10"
NEWLABEL="$SNP"
M1=`grep "$SNP" p312/short-report.csv | head -1 | awk '{print $16+0}' FS=,`
M2=`grep "$SNP" p312/short-report.csv | head -1 | awk '{print $6+$16-1}' FS=,`
C1=`echo "$M1" | awk '{print $1-17}'`
C2=`echo "$M2" | awk '{print $1-17+1}'`
echo "$M1 $M2" 
awk 'NR<16 || $1+0==0 || $2~"P312" || $2~"L21" || $2~"L459" || $2~"DF13" || $16>=m1 && $16<=m2' m1="$M1" m2="$M2" FS=, p312/short-report.csv | cut -d, -f1-17,"$M1"-"$M2" > p312/short-report-l21-zz10.csv

awk '$0~"Row 1" {intable=1} \
$0~"</table>" {intable=0} \
intable==0 {print} \
intable==1 {alldata=$0} \
intable==1 && $0~"<tr>" {print} \
intable==1 && ($0~"<td col" || $0~"<TD col") {split($2,x,"="); split($NF,y,"="); gsub(/"/,"",y[2]); gsub(/>/,"",y[2]); split(y[2],z,","); \
    if ((z[2]<m1 && z[3]>m1) || (z[2]<m2 && z[3]>m2)) {z[2]=z[2]<m1?m1:z[2]; z[3]=z[3]>m2?m2:z[3]; newcols=z[3]-z[2]; old2=new2=$2; gsub(/[0-9]/,"",new2); new2=new2""newcols; gsub(old2,new2,alldata); print "<!",new2,">"} \
    if (z[2]>=m1 && z[3]<=m2 && z[3]-z[2]>0) {print alldata; if ($0!~"null") n=1}} \
intable==1 && $0~"</td>" && n==1 {print; n=0} \
intable==1 && $0~"<td bg" {split($11,x,"="); gsub(/"/,"",x[2]); if (x[2]+0>=m1 && x[2]+0<m2) print}' m1="$C1" m2="$C2" p312/report.html > p312/report-l21-zz10.html



SNP="DF49"
NEWLABEL="$SNP"
M1=`grep "$SNP" p312/short-report.csv | head -1 | awk '{print $16+0}' FS=,`
M2=`grep "$SNP" p312/short-report.csv | head -1 | awk '{print $6+$16-1}' FS=,`
C1=`echo "$M1" | awk '{print $1-17}'`
C2=`echo "$M2" | awk '{print $1-17+1}'`
echo "$M1 $M2"
awk 'NR<16 || $1+0==0 || $2~"P312" || $2~"L21" || $2~"L459" || $2~"DF13" || $16>=m1 && $16<=m2' m1="$M1" m2="$M2" FS=, p312/short-report.csv | cut -d, -f1-17,"$M1"-"$M2" > p312/short-report-l21-df49.csv

awk '$0~"Row 1" {intable=1} \
$0~"</table>" {intable=0} \
intable==0 {print} \
intable==1 {alldata=$0} \
intable==1 && $0~"<tr>" {print} \
intable==1 && ($0~"<td col" || $0~"<TD col") {split($2,x,"="); split($NF,y,"="); gsub(/"/,"",y[2]); gsub(/>/,"",y[2]); split(y[2],z,","); \
    if ((z[2]<m1 && z[3]>m1) || (z[2]<m2 && z[3]>m2)) {z[2]=z[2]<m1?m1:z[2]; z[3]=z[3]>m2?m2:z[3]; newcols=z[3]-z[2]; old2=new2=$2; gsub(/[0-9]/,"",new2); new2=new2""newcols; gsub(old2,new2,alldata); print "<!",new2,">"} \
    if (z[2]>=m1 && z[3]<=m2 && z[3]-z[2]>0) {print alldata; if ($0!~"null") n=1}} \
intable==1 && $0~"</td>" && n==1 {print; n=0} \
intable==1 && $0~"<td bg" {split($11,x,"="); gsub(/"/,"",x[2]); if (x[2]+0>=m1 && x[2]+0<m2) print}' m1="$C1" m2="$C2" p312/report.html > p312/report-l21-df49.html



SNP="DF21"
NEWLABEL="$SNP"
M1=`grep "$SNP" p312/short-report.csv | head -1 | awk '{print $16+0}' FS=,`
M2=`grep "$SNP" p312/short-report.csv | head -1 | awk '{print $6+$16-1}' FS=,`
C1=`echo "$M1" | awk '{print $1-17}'`
C2=`echo "$M2" | awk '{print $1-17+1}'`
echo "$M1 $M2" 
awk 'NR<16 || $1+0==0 || $2~"P312" || $2~"L21" || $2~"L459" || $2~"DF13" || $16>=m1 && $16<=m2' m1="$M1" m2="$M2" FS=, p312/short-report.csv | cut -d, -f1-17,"$M1"-"$M2" > p312/short-report-l21-df21.csv

awk '$0~"Row 1" {intable=1} \
$0~"</table>" {intable=0} \
intable==0 {print} \
intable==1 {alldata=$0} \
intable==1 && $0~"<tr>" {print} \
intable==1 && ($0~"<td col" || $0~"<TD col") {split($2,x,"="); split($NF,y,"="); gsub(/"/,"",y[2]); gsub(/>/,"",y[2]); split(y[2],z,","); \
    if ((z[2]<m1 && z[3]>m1) || (z[2]<m2 && z[3]>m2)) {z[2]=z[2]<m1?m1:z[2]; z[3]=z[3]>m2?m2:z[3]; newcols=z[3]-z[2]; old2=new2=$2; gsub(/[0-9]/,"",new2); new2=new2""newcols; gsub(old2,new2,alldata); print "<!",new2,">"} \
    if (z[2]>=m1 && z[3]<=m2 && z[3]-z[2]>0) {print alldata; if ($0!~"null") n=1}} \
intable==1 && $0~"</td>" && n==1 {print; n=0} \
intable==1 && $0~"<td bg" {split($11,x,"="); gsub(/"/,"",x[2]); if (x[2]+0>=m1 && x[2]+0<m2) print}' m1="$C1" m2="$C2" p312/report.html > p312/report-l21-df21.html



SNP="L513"
NEWLABEL="$SNP"
M1=`grep "$SNP" p312/short-report.csv | head -1 | awk '{print $16+0}' FS=,`
M2=`grep "$SNP" p312/short-report.csv | head -1 | awk '{print $6+$16-1}' FS=,`
C1=`echo "$M1" | awk '{print $1-17}'`
C2=`echo "$M2" | awk '{print $1-17+1}'`
echo "$M1 $M2" 
awk 'NR<16 || $1+0==0 || $2~"P312" || $2~"L21" || $2~"L459" || $2~"DF13" || $16>=m1 && $16<=m2' m1="$M1" m2="$M2" FS=, p312/short-report.csv | cut -d, -f1-17,"$M1"-"$M2" > p312/short-report-l21-l513.csv

awk '$0~"Row 1" {intable=1} \
$0~"</table>" {intable=0} \
intable==0 {print} \
intable==1 {alldata=$0} \
intable==1 && $0~"<tr>" {print} \
intable==1 && ($0~"<td col" || $0~"<TD col") {split($2,x,"="); split($NF,y,"="); gsub(/"/,"",y[2]); gsub(/>/,"",y[2]); split(y[2],z,","); \
    if ((z[2]<m1 && z[3]>m1) || (z[2]<m2 && z[3]>m2)) {z[2]=z[2]<m1?m1:z[2]; z[3]=z[3]>m2?m2:z[3]; newcols=z[3]-z[2]; old2=new2=$2; gsub(/[0-9]/,"",new2); new2=new2""newcols; gsub(old2,new2,alldata); print "<!",new2,">"} \
    if (z[2]>=m1 && z[3]<=m2 && z[3]-z[2]>0) {print alldata; if ($0!~"null") n=1}} \
intable==1 && $0~"</td>" && n==1 {print; n=0} \
intable==1 && $0~"<td bg" {split($11,x,"="); gsub(/"/,"",x[2]); if (x[2]+0>=m1 && x[2]+0<m2) print}' m1="$C1" m2="$C2" p312/report.html > p312/report-l21-l513.html


SNP="FGC11134"
NEWLABEL="$SNP"
M1=`grep "$SNP" p312/short-report.csv | head -1 | awk '{print $16+0}' FS=,`
M2=`grep "$SNP" p312/short-report.csv | head -1 | awk '{print $6+$16-1}' FS=,`
C1=`echo "$M1" | awk '{print $1-17}'`
C2=`echo "$M2" | awk '{print $1-17+1}'`
echo "$M1 $M2" 
awk 'NR<16 || $1+0==0 || $2~"P312" || $2~"L21" || $2~"L459" || $2~"DF13" || $16>=m1 && $16<=m2' m1="$M1" m2="$M2" FS=, p312/short-report.csv | cut -d, -f1-17,"$M1"-"$M2" > p312/short-report-l21-fgc11134.csv

awk '$0~"Row 1" {intable=1} \
$0~"</table>" {intable=0} \
intable==0 {print} \
intable==1 {alldata=$0} \
intable==1 && $0~"<tr>" {print} \
intable==1 && ($0~"<td col" || $0~"<TD col") {split($2,x,"="); split($NF,y,"="); gsub(/"/,"",y[2]); gsub(/>/,"",y[2]); split(y[2],z,","); \
    if ((z[2]<m1 && z[3]>m1) || (z[2]<m2 && z[3]>m2)) {z[2]=z[2]<m1?m1:z[2]; z[3]=z[3]>m2?m2:z[3]; newcols=z[3]-z[2]; old2=new2=$2; gsub(/[0-9]/,"",new2); new2=new2""newcols; gsub(old2,new2,alldata); print "<!",new2,">"} \
    if (z[2]>=m1 && z[3]<=m2 && z[3]-z[2]>0) {print alldata; if ($0!~"null") n=1}} \
intable==1 && $0~"</td>" && n==1 {print; n=0} \
intable==1 && $0~"<td bg" {split($11,x,"="); gsub(/"/,"",x[2]); if (x[2]+0>=m1 && x[2]+0<m2) print}' m1="$C1" m2="$C2" p312/report.html > p312/report-l21-fgc11134.html



SNP="L1335"
NEWLABEL="$SNP"
M1=`grep "$SNP" p312/short-report.csv | head -1 | awk '{print $16+0}' FS=,`
M2=`grep "$SNP" p312/short-report.csv | head -1 | awk '{print $6+$16-1}' FS=,`
C1=`echo "$M1" | awk '{print $1-17}'`
C2=`echo "$M2" | awk '{print $1-17+1}'`
echo "$M1 $M2" 
awk 'NR<16 || $1+0==0 || $2~"P312" || $2~"L21" || $2~"L459" || $2~"DF13" || $16>=m1 && $16<=m2' m1="$M1" m2="$M2" FS=, p312/short-report.csv | cut -d, -f1-17,"$M1"-"$M2" > p312/short-report-l21-l1335.csv

awk '$0~"Row 1" {intable=1} \
$0~"</table>" {intable=0} \
intable==0 {print} \
intable==1 {alldata=$0} \
intable==1 && $0~"<tr>" {print} \
intable==1 && ($0~"<td col" || $0~"<TD col") {split($2,x,"="); split($NF,y,"="); gsub(/"/,"",y[2]); gsub(/>/,"",y[2]); split(y[2],z,","); \
    if ((z[2]<m1 && z[3]>m1) || (z[2]<m2 && z[3]>m2)) {z[2]=z[2]<m1?m1:z[2]; z[3]=z[3]>m2?m2:z[3]; newcols=z[3]-z[2]; old2=new2=$2; gsub(/[0-9]/,"",new2); new2=new2""newcols; gsub(old2,new2,alldata); print "<!",new2,">"} \
    if (z[2]>=m1 && z[3]<=m2 && z[3]-z[2]>0) {print alldata; if ($0!~"null") n=1}} \
intable==1 && $0~"</td>" && n==1 {print; n=0} \
intable==1 && $0~"<td bg" {split($11,x,"="); gsub(/"/,"",x[2]); if (x[2]+0>=m1 && x[2]+0<m2) print}' m1="$C1" m2="$C2" p312/report.html > p312/report-l21-l1335.html



SNP="DF41"
NEWLABEL="$SNP"
M1=`grep "$SNP" p312/short-report.csv | head -1 | awk '{print $16+0}' FS=,`
M2=`grep "$SNP" p312/short-report.csv | head -1 | awk '{print $6+$16-1}' FS=,`
C1=`echo "$M1" | awk '{print $1-17}'`
C2=`echo "$M2" | awk '{print $1-17+1}'`
echo "$M1 $M2" 
awk 'NR<16 || $1+0==0 || $2~"P312" || $2~"L21" || $2~"L459" || $2~"DF13" || $16>=m1 && $16<=m2' m1="$M1" m2="$M2" FS=, p312/short-report.csv | cut -d, -f1-17,"$M1"-"$M2" > p312/short-report-l21-df41.csv

awk '$0~"Row 1" {intable=1} \
$0~"</table>" {intable=0} \
intable==0 {print} \
intable==1 {alldata=$0} \
intable==1 && $0~"<tr>" {print} \
intable==1 && ($0~"<td col" || $0~"<TD col") {split($2,x,"="); split($NF,y,"="); gsub(/"/,"",y[2]); gsub(/>/,"",y[2]); split(y[2],z,","); \
    if ((z[2]<m1 && z[3]>m1) || (z[2]<m2 && z[3]>m2)) {z[2]=z[2]<m1?m1:z[2]; z[3]=z[3]>m2?m2:z[3]; newcols=z[3]-z[2]; old2=new2=$2; gsub(/[0-9]/,"",new2); new2=new2""newcols; gsub(old2,new2,alldata); print "<!",new2,">"} \
    if (z[2]>=m1 && z[3]<=m2 && z[3]-z[2]>0) {print alldata; if ($0!~"null") n=1}} \
intable==1 && $0~"</td>" && n==1 {print; n=0} \
intable==1 && $0~"<td bg" {split($11,x,"="); gsub(/"/,"",x[2]); if (x[2]+0>=m1 && x[2]+0<m2) print}' m1="$C1" m2="$C2" p312/report.html > p312/report-l21-df41.html


SNP="L21"
NEWLABEL="L21 and DF13 minor clades"
M1=`grep "DF41" p312/short-report.csv | head -1 | awk '{print $6+$16}' FS=,`
M2=`grep "$SNP" p312/short-report.csv | head -1 | awk '{print $6+$16-1}' FS=,`
C1=`echo "$M1" | awk '{print $1-17}'`
C2=`echo "$M2" | awk '{print $1-17+1}'`
echo "$M1 $M2" 
echo "$C1 $C2" 
awk 'NR<16 || $1+0==0 || $2~"P312" || $2~"L21" || $16>=m1 && $16<=m2' m1="$M1" m2="$M2" FS=, p312/short-report.csv | cut -d, -f1-17,"$M1"-"$M2" > p312/short-report-l21x.csv

awk '$0~"Row 1" {intable=1} \
$0~"</table>" {intable=0} \
intable==0 {print} \
intable==1 {alldata=$0} \
intable==1 && $0~"<tr>" {print} \
intable==1 && ($0~"<td col" || $0~"<TD col") {split($2,x,"="); split($NF,y,"="); gsub(/"/,"",y[2]); gsub(/>/,"",y[2]); split(y[2],z,","); \
    if ((z[2]<m1 && z[3]>m1) || (z[2]<m2 && z[3]>m2)) {z[2]=z[2]<m1?m1:z[2]; z[3]=z[3]>m2?m2:z[3]; newcols=z[3]-z[2]; old2=new2=$2; gsub(/[0-9]/,"",new2); new2=new2""newcols; gsub(old2,new2,alldata); print "<!",new2,">"} \
    if (z[2]>=m1 && z[3]<=m2 && z[3]-z[2]>0) {print alldata, "<!",z[2], z[3],">"; if ($0!~"null") n=1}} \
intable==1 && $0~"</td>" && n==1 {print; n=0} \
intable==1 && $0~"<td bg" {split($11,x,"="); gsub(/"/,"",x[2]); if (x[2]+0>=m1 && x[2]+0<m2) print}' m1="$C1" m2="$C2" p312/report.html > p312/report-l21x.html

SNP="P312"
NEWLABEL="P312 minor clades, including DF27 and U152"
M1=`grep "L21" p312/short-report.csv | head -1 | awk '{print $6+$16}' FS=,`
M2=`head -1 p312/short-report.csv | awk '{print NF}' FS=,`
C1=`echo "$M1" | awk '{print $1-17}'`
C2=`echo "$M2" | awk '{print $1-17+1}'`
echo "$M1 $M2" 
echo "$C1 $C2" 
awk 'NR<16 || $1+0==0 || $2~"P312" || $16>=m1 && $16<=m2' m1="$M1" m2="$M2" FS=, p312/short-report.csv | cut -d, -f1-17,"$M1"-"$M2" > p312/short-report-p312x.csv

awk '$0~"Row 1" {intable=1} \
$0~"</table>" {intable=0} \
intable==0 {print} \
intable==1 {alldata=$0} \
intable==1 && $0~"<tr>" {print} \
intable==1 && ($0~"<td col" || $0~"<TD col") {split($2,x,"="); split($NF,y,"="); gsub(/"/,"",y[2]); gsub(/>/,"",y[2]); split(y[2],z,","); \
    if ((z[2]<m1 && z[3]>m1) || (z[2]<m2 && z[3]>m2)) {z[2]=z[2]<m1?m1:z[2]; z[3]=z[3]>m2?m2:z[3]; newcols=z[3]-z[2]; old2=new2=$2; gsub(/[0-9]/,"",new2); new2=new2""newcols; gsub(old2,new2,alldata); print "<!",new2,">"} \
    if (z[2]>=m1 && z[3]<=m2 && z[3]-z[2]>0) {print alldata, "<!",z[2], z[3],">"; if ($0!~"null") n=1}} \
intable==1 && $0~"</td>" && n==1 {print; n=0} \
intable==1 && $0~"<td bg" {split($11,x,"="); gsub(/"/,"",x[2]); if (x[2]+0>=m1 && x[2]+0<m2) print}' m1="$C1" m2="$C2" p312/report.html > p312/report-p312x.html
