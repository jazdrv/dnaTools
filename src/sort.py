# license + libs{{{

# Purpose: Y-DNA NGS analytics
# Git repo: https://github.com/jazdrv/dnaTools
# For free distribution under the terms of the GNU General Public License,
# version 3 (29 June 2007) https://www.gnu.org/licenses/gpl.html

import sys,os,yaml,csv,json,time,numpy as np
from beautifultable import BeautifulTable
from collections import OrderedDict
import pandas as pd
from array_api import *

#}}}

try:
    config = yaml.load(open(os.environ['REDUX_CONF_ZAK']))
except:
    print("Missing environment variable REDUX_CONF_ZAK. Aborting.")
    sys.exit()
sys.path.append(config['REDUX_PATH'])

def l2s(lst):
    return ",".join(str(x) for x in lst)

class Variant(object):

    def __init__(self):
        self.dbo = None

    def proc(self,vname):
        #Note: just process one variant 
        self.sort.restore_mx_data()
        self.vix = self.proc_vname(vname)
        if self.vix is None:
            print ("No variants could be found. Exiting.")
            sys.exit()
        self.proc_chk(allowImperfect=False)
        
    def proc_nonsplits(self):
        #Note: calls proc_chk on all imperfect_variants, and processes all
        #variants that don't have any split inconsistencies

        self.sort.restore_mx_data()
        vixs = self.sort.get_imperfect_variants_idx()
        
        if len(vixs):

            #need to push to something to something other than vix here...
            vids = self.sort.get_vid_by_vix(vixs)
            cntErr = 0

            #and then convert it back, because the loop will keep chopping away vixes
            for v in vids:
                self.vix = self.sort.get_vixs_by_vids(v)
                newErr = self.proc_chk(allowImperfect=False,auto_nonsplits=True)
                cntErr = cntErr + newErr
                if newErr == 0:
                    print("No errors - so running: upd_unk, remove_dupes, vandh_sort, save_mx")
                    self.upd_unk(auto_nonsplits=True)
                    self.sort.mx_remove_dupes()
                    self.sort.mx_vandh_sort()
                    self.sort.save_mx()
            if cntErr>0:
                print("\n%s consistency problem(s) seen. Please resolve.\n"%cntErr)
            else:
                print("Consistency check: OK.\n")

            #matrix counts
            print("\nTotal matrix kits: %s" % len(self.sort.KITS))
            print("Total perfect matrix variants: %s" % len(self.sort.get_perfect_variants_idx()))
            print("Total imperfect matrix variants: %s" % len(self.sort.get_imperfect_variants_idx()))
            print("Total matrix variants: %s\n" % len(self.sort.VARIANTS))
        
    def info(self,vname):
        self.sort.restore_mx_data()
        allowImperfect = config['allowImperfectWithVInfo']
        self.vix = self.proc_vname(vname)
        self.set_info(lev=2,allowImperfect=allowImperfect)
        self.stdout_info()
        
    def matrix(self,argL,perfOnly=False,imperfOnly=False):
        #Note: calls stdout_matrix

        self.sort.restore_mx_data()

        #Note: perfect/imperfect only matrix
        if perfOnly:
            vix = self.sort.get_perfect_variants_idx()
        elif imperfOnly:
            vix = self.sort.get_imperfect_variants_idx()
        if perfOnly or imperfOnly:
            if len(argL)>0:
                outK = argL[0]
                outKA = outK.split(",")
                if str(self.proc_kname(outKA[0])) != outKA[0]:
                    kix = self.sort.get_kix_by_name(outKA)
                else:
                    kix = list(map(int, outKA))
                if kix is None:
                    print ("No kits could be found. Exiting.")
                    sys.exit()
            else: #default (show everything)
                kix = None

        #Note: col/row custom matrix (if applicable)
        else:
            if len(argL)>0:
                outV = argL[0]
                outVA = outV.split(",")
                if outVA[0] == '-':
                    vix = None
                elif str(self.proc_vname(outVA[0])) != outVA[0]:
                    vix = self.sort.get_vix_by_name(outVA)
                else:
                    vix = list(map(int, outVA))
                if vix is None:
                    print ("No variants could be found. Exiting.")
                    sys.exit()
            else: #default (show everything)
                vix = None
            if len(argL)>1:
                outK = argL[1]
                outKA = outK.split(",")
                if str(self.proc_kname(outKA[0])) != outKA[0]:
                    kix = self.sort.get_kix_by_name(outKA)
                else:
                    kix = list(map(int, outKA))
                if kix is None:
                    print ("No kits could be found. Exiting.")
                    sys.exit()
            else: #default (show everything)
                kix = None
        self.sort.stdout_matrix(vix=vix,kix=kix)
    def clade_priority(self,argL):
        #argL parsing
        argL = [x.upper() for x in argL]
        for a in argL[:]:
            if a.find(','):
                argL = argL + a.split(",")
                argL.remove(a)
            if a.find('/'):
                argL = argL + a.split("/")
                argL.remove(a)

        if len(argL) > 0:
            self.sort.restore_mx_data()
            for x in argL:

                #get relevant values out of argL
                argA = x.split("^")
                if len(argA) != 2 or argA[1].isdigit() is False or argA[0].isdigit() is True:
                    print("\nImproper format. Exiting.\n")
                    sys.exit()
                new_snpname = argA[0]
                new_vID = int(argA[1])

                #delete any old mx_clade_priorities refs  (this tbl might not be necessary)
                sql = "delete from mx_clade_priorities where snpname = '%s' and vID = %s " % (new_snpname,new_vID)
                self.dbo.sql_exec(sql)
                sql = '''
                    insert into mx_clade_priorities (ID,snpname,vID) select
                    null,snpname,vID from snpnames where snpname = '%s' and vID='%s';
                    ''' % (new_snpname,new_vID)
                self.dbo.sql_exec(sql)

                #this checks for a different snpname/vID combo
                sql = '''
                    select distinct D.vID,V.pos as dupe_pos
                    from snpnames S, mx_dupe_variants D,variants V
                    where D.dupe_vID=S.vID and S.snpname = '%s' and S.vID = %s
                    and V.ID = D.dupe_vID
                    ''' % (new_snpname,new_vID)
                self.dbo.sql_exec(sql)
                F = self.dbo.fetchall()
                if len(F) > 0:
                    new_vid = True
                else:

                    #this checks for same vID, but different snpname
                    vix = self.sort.get_vixs_by_vids(new_vID)
                    vn = self.sort.get_vname_by_vix(vix)
                    if vn == new_snpname:
                        #Nothing to do, this is the existing setting 
                        return 
                    sql = '''
                        select distinct D.vID,V.pos as new_pos
                        from snpnames S, mx_dupe_variants D,variants V
                        where D.vID=S.vID and S.snpname = '%s' and S.vID = %s
                        and V.ID = D.vID
                        ''' % (new_snpname,new_vID)
                    self.dbo.sql_exec(sql)
                    F = self.dbo.fetchall()
                    if len(F) > 0:
                        new_vid = False

                if len(F) > 0:

                    #vars
                    old_vID = F[0][0]
                    new_vpos = F[0][1]
                    vix = self.sort.get_vixs_by_vids(old_vID)
                    old_snpname = self.sort.get_vname_by_vix(vix)

                    #update self.sort.VARIANTS
                    self.sort.VARIANTS.pop(old_snpname,None)
                    self.sort.VARIANTS[new_snpname] = (int(new_vID),vix,new_vpos)

                    #TODO: perhaps this sups/subs code should really be refreshing all
                    #the relations from scratch for this clade
                    #update mx_sups_subs - sup references
                    if new_vid:
                        sql = "update mx_sups_subs set sup = %s where sup = %s" % (new_vID,old_vID)
                        self.dbo.sql_exec(sql)

                    #update mx_sups_subs - sub references
                    if new_vid:
                        sql = "update mx_sups_subs set sub = %s where sub = %s" % (old_vID,new_vID)
                        self.dbo.sql_exec(sql)

                    #update panda version of mx_sups_subs
                    if new_vid:
                        sql = "select distinct * from mx_sups_subs;"
                        self.dbo.sql_exec(sql)
                        sups_ = self.dbo.fetchall()
                        self.SS = pd.DataFrame(sups_,columns=['sup','sub'])

                    #update mx_dupe_variants - vID references
                    if new_vid:
                        sql = "update mx_dupe_variants set vID = %s where vID = %s" % (new_vID,old_vID)
                        self.dbo.sql_exec(sql)

                    #update mx_dupe_variants - dupe_vID references
                    if new_vid:
                        sql = "update mx_dupe_variants set dupe_vID = %s where dupe_vID = %s" % (old_vID,new_vID)
                        self.dbo.sql_exec(sql)

                    #update mx_variants table
                    sql = "delete from mx_variants;"
                    self.dbo.sql_exec(sql)
                    sql = "insert into mx_variants (ID,name,pos) values (?,?,?);"
                    self.dbo.sql_exec_many(sql,[(tuple([vid,nm,pos])) for (n,(nm,(vid,idx,pos))) in enumerate(self.sort.get_axis('variants'))])

                    #update mx_idxs table
                    if new_vid:
                        sql = "update mx_idxs set axis_id = %s where axis_id= %s and type_id=0" % (new_vID,old_vID)
                        self.dbo.sql_exec(sql)

                    #reset mx_sort_recommendations
                    if new_vid:
                        sql = "delete from mx_sort_recommendations;"
                        self.dbo.sql_exec(sql)

                    print("Done making changes. Check the matrix and clade detail views to make sure it worked!")

    def ref(self,argL,name=False,pos=False,id=False,vix=False,strT=False,clade=False,calls=False,hg19=True):
        #argL parsing
        argL = [x.upper() for x in argL]
        for a in argL[:]:
            if a.find(','):
                argL = argL + a.split(",")
                argL.remove(a)
            if a.find('/'):
                argL = argL + a.split("/")
                argL.remove(a)

        #strT/sqlw clause
        if strT:
            sqlw_ = "'"+"','".join(str(x) for x in sorted(list(set(argL))))+"'"
        else:
            sqlw_ = ",".join(str(x) for x in sorted(list(set(argL))))

        #ref_clade (if applicable)
        if clade:
            self.ref_clade(name=name,pos=pos,id=id,vix=vix,sqlw_=sqlw_)
            return

        #sqlw clause
        if name: sqlw = "RV.snpname in (%s)" % sqlw_
        elif pos: sqlw = "RV.pos in (%s)" % sqlw_
        elif id: sqlw = "RV.ID in (%s)" %  sqlw_
        elif vix: sqlw = "RV.idx in (%s)" % sqlw_

        #ref_calls (if applicable)
        if calls:
            self.ref_calls(sqlw)
            return

        #get the data
        sql = '''
            select distinct RV.snpname, RV.ID, RV.pos, RV.buildNm, RV.anc,
            RV.der, RV.idx, RV.vID1, RV.vID2, RV.reasonId
            from v_ref_variants RV where %s order by 7,1,3;
            ''' % sqlw
        self.dbo.sql_exec(sql)
        F = self.dbo.fetchall()

        #stdout hg38 tbl data
        if len(F) > 0:
            print("")
            table = BeautifulTable(max_width=90)
            table.column_headers = ['vix']+['build']+['name']+['id']+['pos']+['anc']+['der']+['dupeP']+['nouse']
            for row in F:
                if row[7] == None and row[8] == None:
                    dupeP = 'N'
                elif row[7] != None:
                    dupeP = 'Y'
                elif row[8] != 'None':
                    dupeP = row[8]
                if row[9] == None:
                    nouse = '-'
                elif row[9] == 1:
                    nouse = 'P'
                elif row[9] == -1:
                    nouse = 'N'
                table.append_row([str(row[6]).replace('9999999','-')]+[str(row[3]).replace('None','-')]+[str(row[0]).replace('None','-')]+[str(row[1])]+[row[2]]+[row[4]]+[row[5]]+[dupeP]+[nouse])
                table.row_seperator_char = ''
                table.column_seperator_char = ''
                table.column_alignments['name'] = BeautifulTable.ALIGN_LEFT
            print(table)
            print("")

        #stdout hg19 tbl
        if hg19:
            sql = '''
                select distinct RV.snpname,RV.ID,RV.pos,RV.buildNm from v_ref_variants_hg19 RV
                where %s order by 1;
                ''' % sqlw
            self.dbo.sql_exec(sql)
            F = self.dbo.fetchall()
            if len(F) > 0:
                table = BeautifulTable()
                table.column_headers = ['build']+['name']+['id']+['pos']
                for row in F:
                    table.append_row([row[3]]+[row[0]]+[row[1]]+[row[2]])
                    table.row_seperator_char = ''
                    table.column_seperator_char = ''
                print(table)
        
    def ref_clade(self,sqlw_,name=False,pos=False,id=False,vix=False):
        self.sort.restore_mx_data()

        #sqlw clause
        if name:
            sqlw = "S.snpname in (%s)" % sqlw_
        elif pos:
            sqlw = "V.pos in (%s)" % sqlw_
        elif id:
            sqlw = "V.ID in (%s)" %  sqlw_
        elif vix:
            sqlw = "IX.idx in (%s)" % sqlw_

        #get vid - base ref
        sql = '''
            select distinct V.id from build B, variants V
            left join mx_idxs IX on IX.axis_id = V.ID and IX.type_id = 0
            left join snpnames S on S.vID = v.ID
            where B.buildNm = 'hg38' and V.buildID = B.ID and %s
            ''' % (sqlw)
        self.dbo.sql_exec(sql)
        vids1 = [i[0] for i in list(self.dbo.fetchall())]

        #get vid - imx_dupe_variants (dupe_vId ref)
        sql = '''
            select distinct vid from mx_dupe_variants D where dupe_vID in (%s)
            '''% l2s(vids1)
        self.dbo.sql_exec(sql)
        vids2 = [i[0] for i in list(self.dbo.fetchall())]

        #get vid - imx_dupe_variants (vId ref)
        sql = '''
            select distinct dupe_vid from mx_dupe_variants D where vID in (%s)
            '''% l2s(vids1)
        self.dbo.sql_exec(sql)
        vids3 = [i[0] for i in list(self.dbo.fetchall())]

        #get all vid supporting details
        vids = list(set(vids1+vids2+vids3))
        sql = '''
            select distinct RV.snpname, RV.ID, RV.pos, RV.buildNm, RV.anc,
            RV.der, RV.idx, RV.vID1, RV.vID2, RV.reasonId, RV.name
            from v_ref_variants RV where RV.ID in (%s) order by 7,1,3;
            ''' % l2s(vids)
        self.dbo.sql_exec(sql)
        F = self.dbo.fetchall()

        #stdout tbl
        if len(F) > 0:

            cnt = 0
            for row in F:
                row = list(F[cnt])
                #vix handling
                if isinstance(row[6],int) and row[6] != 9999999 and row[0] is not None and self.sort.get_vname_by_vix(row[6]) != row[0]:
                    row[6] = 9999999
                F[cnt] = tuple(list(row))
                cnt = cnt + 1        
            #F.sort(key=lambda x: x[0])
            F.sort(key=lambda x: x[6])

            print("")
            table = BeautifulTable(max_width=80)
            cols = ['vix'] + ['build'] + ['name'] + ['id']
            cols = cols + ['pos'] + ['anc'] + ['der'] + ['dupeP'] + ['nouse']
            table.column_headers = cols

            for row in F:
                if row[7] == None and row[8] == None:
                    dupeP = 'N'
                elif row[7] != None:
                    dupeP = 'Y'
                elif row[8] != 'None':
                    dupeP = row[8]
                if row[9] == None:
                    nouse = '-'
                elif row[9] == 1:
                    nouse = 'P'
                elif row[9] == -1:
                    nouse = 'N'
                row_ = [str(row[6]).replace('9999999','-')] + [str(row[3]).replace('None','-')]
                row_ = row_ + [str(row[0]).replace('None','-')] + [str(row[1])]
                row_ = row_ + [row[2]] + [row[4]] + [row[5]] + [dupeP] + [nouse]
                table.append_row(row_)
                table.row_seperator_char = ''
                table.column_seperator_char = ''
                table.column_alignments['name'] = BeautifulTable.ALIGN_LEFT
            print(table)
            print("")
        
    def ref_calls(self,sqlw):
        self.sort.restore_mx_data()

        #get call data
        sql = '''
            select distinct
            RV.snpname, RV.ID, RV.pos, RV.buildNm, RV.anc,
            RV.der, RV.idx, RV.vID1, RV.vID2, RV.reasonId,
            T.pID, C.assigned, C.genotype, T.val, RV.name, D.kitId
            from v_ref_variants RV, tmp1 T, dataset D
            left join vcfcalls C
            on C.vID = RV.ID and C.pID = T.pID
            where D.ID = T.pID and RV.pos = T.pos and %s -- and T.pID is not None
            order by 1
            ''' % (sqlw)
        self.dbo.sql_exec(sql)
        F = self.dbo.fetchall()

        if len(F) > 0:
            cnt = 0

            #some pre-sort hacking to get the sort setup
            for row in F:
                row = list(F[cnt])

                #kix
                if isinstance(row[10],int):
                    kix = self.sort.get_kixs_by_kids(row[10])
                else:
                    kix = 9999999 #just make a big value for sorting purposes
                #coord
                if isinstance(row[6],int) and isinstance(kix,int) and row[6] != 9999999 and kix != 9999999:
                    coord = self.sort.NP[row[6],kix]
                else:
                    coord = '-'
                #vix handling
                if isinstance(row[6],int) and row[6] != 9999999 and row[0] is not None and self.sort.get_vname_by_vix(row[6]) != row[0]:
                    row[6] = 9999999

                if row[0] is None:
                    row[0] = '-'

                F[cnt] = tuple(list(row) + [kix]+[coord]) #16,17
                cnt = cnt + 1        
                
            #sorting
            F.sort(key=lambda x: x[16])
            F.sort(key=lambda x: x[0])
            F.sort(key=lambda x: x[6])

            #stdout tbl
            print("")
            table = BeautifulTable(max_width=110)
            cols = ['vix'] + ['build'] + ['name']
            cols = cols + ['id'] + ['pos'] + ['anc'] + ['der'] + ['dupeP'] + ['nouse']
            cols = cols + ['pID'] + ['kix'] + ['kit'] + ['assign'] + ['geno'] + ['bed'] + ['mxval']
            table.column_headers = cols
            for row in F:
                if row[7] == None and row[8] == None:
                    dupeP = 'N'
                elif row[7] != None:
                    dupeP = 'Y'
                elif row[8] != 'None':
                    dupeP = row[8]
                if row[9] == None:
                    nouse = '-'
                elif row[9] == 1:
                    nouse = 'P'
                elif row[9] == -1:
                    nouse = 'N'
                if row[6] == 9999999:
                    vix = '-'
                else:
                    vix = row[6]
                if row[16] == 9999999:
                    kix = '-'
                else:
                    kix = row[16]
                row_ = [str(vix)] + [str(row[3]).replace('None','-')] + [str(row[0]).replace('None','-')]
                row_ = row_ + [str(row[1])] + [row[2]] + [row[4]] + [row[5]] + [dupeP] + [nouse]
                row_ = row_ + [str(row[10])] + [str(kix)] + [str(row[15])] + [str(row[11])] + [str(row[12])] + [str(row[13])] + [row[17]]
                table.append_row(row_)
                table.row_seperator_char = ''
                table.column_seperator_char = ''
                table.column_alignments['name'] = BeautifulTable.ALIGN_LEFT
            print(table)
            print("")
        
    def ref_name(self,argL,clade=False,calls=False):
        self.ref(argL,name=True,strT=True,clade=clade,calls=calls)
        
    def ref_pos(self,argL,clade=False,calls=False):
        self.ref(argL,pos=True,clade=clade,calls=calls)
                 
    def ref_id(self,argL,clade=False,calls=False):
        self.ref(argL,id=True,clade=clade,calls=calls)
        
    def ref_vix(self,argL,clade=False,calls=False):
        self.ref(argL,vix=True,clade=clade,calls=calls,hg19=False)

    def stash(self,vname):
        #TODO: setup the ability to stash multiple variants at once
        self.sort.restore_mx_data()
        if vname.isdigit() and int(vname)<=len(self.sort.VARIANTS):
            print("\n** Assuming, you're providing a vix ID...")
            vix = int(vname)
        elif vname.isdigit():
            vix = self.sort.get_vix_by_name(vname)
        else:
            vix = self.sort.get_vix_by_name(vname.upper())
        if vix is None:
            print("Invalid variant given. Exiting.")
            sys.exit()

        #vars
        vid = self.sort.get_vid_by_vix(vix)
        vname = self.sort.get_vname_by_vix(vix)

        if vid:

            #move vid to stash table
            sql = "insert into mx_variant_stash (ID) values (%s)" % vid
            self.dbo.sql_exec(sql)

            #mx_remove_vix routine
            self.sort.mx_remove_vix(vix)
            print("Moved to stash: %s [%s]"%(vname,vix))
        
    def unstash(self,vname):
        #TODO: need to set this up
        print("Unstash isn't setup yet. Exiting.")
        sys.exit()

        self.sort.restore_mx_data()
        if vname.isdigit() and int(vname)<=len(self.sort.VARIANTS):
            print("\n** Assuming, you're providing a vix ID...")
            vix = int(vname)
        elif vname.isdigit():
            vix = self.sort.get_vix_by_name(vname)
        else:
            vix = self.sort.get_vix_by_name(vname.upper())
        #move this variant back from stash tbl

    def proc_vname(self,vname):
        if vname.isdigit() and int(vname)<=len(self.sort.VARIANTS):
            print("\n** Assuming, you're providing a vix ID...")
            vix = int(vname)
        elif vname.isdigit():
            vix = self.sort.get_vix_by_name(vname)
        else:
            vix = self.sort.get_vix_by_name(vname.upper())
        return vix
        
    def proc_kname(self,kname):
        if kname.isdigit() and int(kname)<=len(self.sort.KITS):
            print("\n** Assuming, you're providing a kix ID...")
            kix = int(kname)
        elif kname.isdigit():
            kix = self.sort.get_kix_by_name(kname)
        else:
            kix = self.sort.get_kix_by_name(kname.upper())
        return kix
        
    def set_info(self,vix=None,lev=1,allowImperfect=False):
        lev = lev-1
        if vix is not None:
            self.vix = vix
        self.name = self.sort.get_vname_by_vix(self.vix)
        self.vixn = "%s [%s]"%(self.name,self.vix)
        self.kpc = self.sort.get_kixs_by_val(val=1,vix=self.vix)
        self.knc = self.sort.get_kixs_by_val(val=-1,vix=self.vix)
        self.kuc = self.sort.get_kixs_by_val(val=0,vix=self.vix)
        self.kpcn = "kpc: %s [%s]"%(l2s(self.sort.get_kname_by_kix(self.kpc)),l2s(self.kpc))
        self.kncn = "knc: %s [%s]"%(l2s(self.sort.get_kname_by_kix(self.knc)),l2s(self.knc))
        self.kucn = "kuc: %s [%s]"%(l2s(self.sort.get_kname_by_kix(self.kuc)),l2s(self.kuc))
        self.kpnc = sorted(self.kpc + self.knc)
        self.sups = self.get_rel(relType=1,allowImperfect=allowImperfect) #not from cache
        self.subs = self.get_rel(relType=-1,allowImperfect=allowImperfect) #not from cache
        self.eqv = self.get_rel(relType=0,allowImperfect=allowImperfect)
        self.subsn = "subs: %s [%s]" %(l2s(self.sort.get_vname_by_vix(self.subs)),l2s(self.subs))
        self.supsn = "sups: %s [%s]" %(l2s(self.sort.get_vname_by_vix(self.sups)),l2s(self.sups))
        self.eqvn = "eqv: %s [%s]"%(l2s(self.sort.get_vname_by_vix(vix=self.eqv)),l2s(self.eqv))
        
    def stdout_info(self,tp=None):
        print("")
        if tp is None:
            print("---------------------------------------------------------------------")
            print("")
        if tp==1:
            print("[+] %s" %self.vixn)
            sp = "    "
        elif tp==-1:
            print("[-] %s" %self.vixn)
            sp = "    "
        else:
            print("vix: %s" %self.vixn)
            sp = ""
        print("%s%s" %(sp,self.kpcn))
        print("%s%s" %(sp,self.kncn))
        print("%s%s" %(sp,self.kucn))
        print("%s%s" %(sp,self.supsn))
        print("%s%s" %(sp,self.subsn))
        print("%s%s" %(sp,self.eqvn))
        if tp is None:
            print("")
            print("---------------------------------------------------------------------")
            print("")

    def proc_chk(self,allowImperfect,auto_perfVariants=False,auto_nonsplits=False):
        if auto_perfVariants or auto_nonsplits:
            auto_mode = True
        else:
            auto_mode = False

        if auto_mode is not True:
            print("")
            print("---------------------------------------------------------------------")
            print("")

        #if auto_mode is not True or config['SHOW_PROC_CHK_DETAILS']:
        if auto_perfVariants is not True or config['SHOW_PROC_CHK_DETAILS']:
            print("\nvix: %s [%s]"%(self.sort.get_vname_by_vix(self.vix),self.vix))
        if auto_mode is not True or config['SHOW_PROC_CHK_DETAILS']:
            print("")

        #needed vars 
        self.kuc = self.sort.get_kixs_by_val(val=0,vix=self.vix)
        self.eqvs = self.get_rel(relType=0,allowImperfect=allowImperfect)

        #this is used when checking perfect variants (where the subs/sups are cached)
        if 1==2 and auto_perfVariants is True:
            vid = self.sort.get_vid_by_vix(self.vix)
            try:
                self.subs = self.sort.get_vixs_by_vids(list(self.sort.SS.loc[self.sort.SS['sup']==vid,['sub']].as_matrix().T[0]))
            except:
                self.subs = []
            try:
                self.sups = self.sort.get_vixs_by_vids(list(self.sort.SS.loc[self.sort.SS['sub']==vid,['sup']].as_matrix().T[0]))
            except:
                self.sups = []
        #this is used for imperfect variants (ones with unk's)
        else:
            self.sups = self.get_rel(relType=1,allowImperfect=allowImperfect)
            self.subs = self.get_rel(relType=-1,allowImperfect=allowImperfect)

        #debugging
        if config['SHOW_PROC_CHK_DETAILS']:
            self.kucn = "kuc: %s [%s]"%(l2s(self.sort.get_kname_by_kix(self.kuc)),l2s(self.kuc))
            self.supsn = "sups: %s [%s]" %(l2s(self.sort.get_vname_by_vix(self.sups)),l2s(self.sups))
            self.subsn = "subs: %s [%s]" %(l2s(self.sort.get_vname_by_vix(self.subs)),l2s(self.subs))
            self.eqvn = "eqv: %s [%s]"%(l2s(self.sort.get_vname_by_vix(vix=self.eqvs)),l2s(self.eqvs))
            print(self.supsn)
            print(self.subsn)
            print(self.eqvn)
            print(self.kucn)
            print("num of sups: %s"%len(self.sups))
            print("")

        #starting unk list
        unkL = self.kuc
        unkD = {}
        for k in unkL:
            unkD[k] = [0,0,0] #0's are defaults
        othK = {}
        splitReq = False

        #-------------------------------
        #1st arg of set : sup test1 result
        #1 - ambiguous promotion to positive
        #2 - hard promotion to positive
        #3 - hard negative
        
        #2nd arg of set : consistency test2 result
        #2+ - if it gets to 2+ ... it means this unk should be positive to appease other variants

        #3rd arg of set : consistency test2 result (split override)
        #1 - means neg (due to split)
        #-------------------------------

        #look for sups
        if len(self.sups) == 0:
            if auto_perfVariants is not True:
                print("This variant has no matrix supsets; creating a top.")
            self.sups.append(-999) #top handling

        if len(self.sups) > 0:
            if config['SHOW_PROC_CHK_DETAILS']:
                print("unks: %s [%s]" % (l2s(self.sort.get_kname_by_kix(self.kuc)),l2s(self.kuc)))
                print("")
                print("There are %s sups to go through.\n"%(len(self.sups)))

            for sup in reversed(self.sups):
                splitL = []
                if config['SHOW_PROC_CHK_DETAILS']:
                    print("sup: %s [%s]" % (self.sort.get_vname_by_vix(sup),sup)) #can handle top

                #this part is only for variants with unknowns
                if len(self.kuc)>0:

                    #(beg)sup check
                    knc4sup = self.sort.get_kixs_by_val(val=-1,vix=sup)
                    self.knc = self.sort.get_kixs_by_val(val=-1,vix=self.vix)
                    diffKnc = list(set(self.knc)-set(knc4sup))
                    if len(diffKnc) > 0:
                        pN = list(set(self.kuc).intersection(set(knc4sup)))
                        if len(pN) > 0:
                            if config['SHOW_PROC_CHK_DETAILS']:
                                print(" - [1] vix unk %s [%s] is/are neg" % (l2s(self.sort.get_kname_by_kix(pN)),l2s(pN)))
                            for k in pN:
                                if k in unkL:
                                    unkL.remove(k)
                                unkD[k][0] = 3 #hard failure
                        else:
                            if config['SHOW_PROC_CHK_DETAILS']:
                                print(" - [1] sup is truly sup to target variant, all unks open to promotion")
                    else:
                        if config['SHOW_PROC_CHK_DETAILS']:
                            print(" - [1] sup is sub/eq/sup to target variant")

                    #ambiguous promotions
                    for k in unkL:
                        if unkD[k] != 3:
                            unkD[k][0] = 1 
                    if config['SHOW_PROC_CHK_DETAILS']:
                        print(" - [1] remaining unk - ambig pos: %s [%s]" %(l2s(self.sort.get_kname_by_kix(unkL)),l2s(unkL)))
                    #(end)sup check

                #(beg)consistency check 

                if sup == -999: #top handling
                    supsubs_ = self.sort.get_perfect_variants_idx()

                else:
                    vid = self.sort.get_vid_by_vix(sup)
                    try:
                        supsubs_ = self.sort.get_vixs_by_vids(list(self.sort.SS.loc[self.sort.SS['sup']==vid,['sub']].as_matrix().T[0]))
                    except:
                        supsubs_ = []

                #create a temporary kpc for the target variant based on expected unk-pos promotions
                pos_ = []
                for k,v in unkD.items():
                    if v[0] == 3 or v[2] == 1:
                        continue
                    elif v[0] == 2 or v[1] > 1:
                        pos_.append(k) #pos
                    elif v[0] == 1:
                        pos_.append(k) #pos ambiguous
                kpc_ = sorted(self.kpc + pos_)
                subs_ = self.get_rel(relType=-1,kpc=kpc_,allowImperfect=allowImperfect) #not from cache
                eqvs_ = self.get_rel(relType=0,kpc=kpc_,allowImperfect=allowImperfect) #not from cache

                #debugging - kpc and temp kpc
                if config['SHOW_PROC_CHK_DETAILS']:
                    print("\nself.kpc: %s [%s]"%(l2s(self.sort.get_kname_by_kix(self.kpc)),l2s(self.kpc)))
                    print("temp kpc (kpc_): %s [%s]"%(l2s(self.sort.get_kname_by_kix(kpc_)),l2s(kpc_)))

                #kpc of the sup
                supkpc = self.sort.get_kixs_by_val(val=1,vix=sup)

                if kpc_ != supkpc: #no sense in doing consistency check, if target variant is looking to be a dupe of its sup

                    #remove any target subs or target equivalent variants from consideration (use cache)
                    supsubs = list(set(supsubs_)-set(subs_)-set(eqvs_))

                    #do consistency checks on remaining variants
                    if config['SHOW_PROC_CHK_DETAILS']:
                        print("\nThere are %s supsubs for sup: %s"%(len(supsubs),self.sort.get_vname_by_vix(sup)))
                    if len(supsubs):
                        if config['SHOW_PROC_CHK_DETAILS']:
                            print("")

                        #self.sort.get_perfect_variants_idx()
                        perfNP = self.sort.NP[self.sort.get_perfect_variants_idx()]
                        VAR1a = np.argwhere(perfNP[:,kpc_] == 1) #find self.kpc overlaps
                        mask = np.isin(VAR1a[:,0], supsubs)
                        idx = list(list(np.where(mask))[0])
                        VAR2a = np.unique(VAR1a[:,0][idx]) #these are the overlapping supsubs with target v
                        if config['SHOW_PROC_CHK_DETAILS']:
                            print("Overlapping supsubs: [%s]\n"%(l2s(VAR2a)))
                        for v in VAR2a:
                            supsubkpc = self.sort.get_kixs_by_val(val=1,vix=v)
                            #(mtP) common kpc btw supsub and target variant
                            at_common_kpc = sorted(list(set(supsubkpc).intersection(set(kpc_))))
                            #(msP1) common kpc btw vix and common sup
                            as_common_kpc1 = sorted(list(set(kpc_).intersection(set(supkpc))))
                            #(msP2) common kpc btw supsub and common sup
                            as_common_kpc2 = sorted(list(set(supsubkpc).intersection(set(kpc_))))
                            #(xP) diff_common_kpc = [msP(1/2)-mtP]
                            diff_common_kpc1 = sorted(list(set(as_common_kpc1)-set(at_common_kpc)))
                            diff_common_kpc2 = sorted(list(set(as_common_kpc2)-set(at_common_kpc)))
                            #(cP) common_kp c = intersection btw msP(1/2) and mtP
                            common_kpc1 = sorted(list(set(as_common_kpc1).intersection(set(at_common_kpc))))
                            common_kpc2 = sorted(list(set(as_common_kpc2).intersection(set(at_common_kpc))))
                            if config['SHOW_PROC_CHK_DETAILS']:
                                print(" - [2] supsub: %s [%s]"%(self.sort.get_vname_by_vix(v),v))
                                print("       supsubkpc: %s [%s]"%(l2s(self.sort.get_kname_by_kix(supsubkpc)),l2s(supsubkpc)))
                                print("       (msP1) shared btw vix + sup: %s [%s]"%(l2s(self.sort.get_kname_by_kix(as_common_kpc1)),l2s(as_common_kpc1)))
                                print("       (msP2) shared btw supsub + sup: %s [%s]"%(l2s(self.sort.get_kname_by_kix(as_common_kpc2)),l2s(as_common_kpc2)))
                                print("       (mtP) shared btw vix + supsub: %s [%s]"%(l2s(self.sort.get_kname_by_kix(at_common_kpc)),l2s(at_common_kpc)))
                                print("       (xP1) msP1-mtP: %s [%s]"%(l2s(self.sort.get_kname_by_kix(diff_common_kpc1)),l2s(diff_common_kpc1)))
                                print("       (xP2) msP2-mtP: %s [%s]"%(l2s(self.sort.get_kname_by_kix(diff_common_kpc2)),l2s(diff_common_kpc2)))
                                print("       (cP1) common btw msP1+mtP: %s [%s]"%(l2s(self.sort.get_kname_by_kix(common_kpc1)),l2s(common_kpc1)))
                                print("       (cP2) common btw msP2+mtP: %s [%s]"%(l2s(self.sort.get_kname_by_kix(common_kpc2)),l2s(common_kpc2)))

                            #go through common_kpc1/diff_common_kpc1
                            for k in diff_common_kpc1:
                                if k in unkL:
                                    unkD[k][1] = unkD[k][1] + 1
                                    if unkD[k][1] > 1:
                                        if config['SHOW_PROC_CHK_DETAILS']:
                                            print(" - [2] %s [%s] is a req positive" % (self.sort.get_kname_by_kix(k),k))
                            if len(common_kpc1) > 0:
                                splitL.append(common_kpc1)
                            if len(diff_common_kpc1) > 0:
                                splitL.append(diff_common_kpc1)
                            if config['SHOW_PROC_CHK_DETAILS']:
                                print("       splitL after commkpc1 checks: %s"%splitL)

                            #go through common_kpc2/diff_common_kpc2
                            for k in diff_common_kpc2:
                                if k in unkL:
                                    unkD[k][1] = unkD[k][1] + 1
                                    if unkD[k][1] > 1:
                                        if config['SHOW_PROC_CHK_DETAILS']:
                                            print(" - [2] %s [%s] is a req positive" % (self.sort.get_kname_by_kix(k),k))
                            if len(common_kpc2) > 0:
                                splitL.append(common_kpc2)
                            if len(diff_common_kpc2) > 0:
                                splitL.append(diff_common_kpc2)
                            if config['SHOW_PROC_CHK_DETAILS']:
                                print("       splitL after commkpc2 checks: %s"%splitL)
                        
                    #splits
                    uniq_splits = [list(x) for x in set(tuple(x) for x in splitL)]
                    if config['SHOW_PROC_CHK_DETAILS']:
                        print("       uniq_splits: %s"%uniq_splits)
                    if len(uniq_splits) > 0:
                        split_intersects = sorted(list(set.intersection(*map(set, uniq_splits))))
                        if config['SHOW_PROC_CHK_DETAILS']:
                            print("       split_intersects: %s"%l2s(split_intersects))
                        if (len(uniq_splits) > 1 and len(split_intersects) == 0):
                            if config['SHOW_PROC_CHK_DETAILS']:
                                print(" - [2] split required: btw %s" % uniq_splits)
                                print(" - [2] all unk to negative")
                            splitReq = True

                    if config['SHOW_PROC_CHK_DETAILS']:
                        print("")

        #resolution vars - for imperfect variants
        if auto_perfVariants is not True:
            pos = []
            pos_a = []
            neg = []
            rec = dict()

        #resolution vars - for imperfect/perfect variants
        lenS = 1
        s_ = []
        spl = []

        #resolutions - splits
        if splitReq:
            for s in uniq_splits:
                if len(s) > 1 : lenS = len(s)
                else: s_.append(s[0])
            if len(uniq_splits)>1: #has to be over 1 for it to be a true split
                #print(split_intersects)
                #print(uniq_splits)
                #print("splits: %s [%s]" % (l2s(self.sort.get_kname_by_kix(s_)),l2s(s_)))
                spl = s_
            #else:
            #    print("splits: %s" % l2s(uniq_splits))
            for k in self.kuc:
                unkD[k][2] = 1

        #resolutions - other pos/pos_a/neg
        if auto_perfVariants is not True:
            for k,v in unkD.items():
                if v[0] == 3 or v[2] == 1:
                    neg.append(k)
                elif v[0] == 2 or v[1] > 1:
                    pos.append(k)
                elif v[0] == 1:
                    pos_a.append(k)

        #mx_sort_recommendation updates

        if auto_perfVariants is True and len(spl) == 0:
            return 0 #no split errors for perfect variant chk

        #split found
        if len(spl) > 0:
            if auto_perfVariants is True:
                print("[SPLIT ISSUE] vix: %s [%s] - splits: %s [%s]"%(self.sort.get_vname_by_vix(self.vix),self.vix,l2s(self.sort.get_kname_by_kix(spl)),l2s(spl)))
                print(uniq_splits)
                return 1 #split errors for perfect variant chk
            elif auto_perfVariants is False:
                print("[SPLIT ISSUE] splits: %s [%s]" % (l2s(self.sort.get_kname_by_kix(spl)),l2s(spl)))
                rec.update({"spl":spl})

        vid = self.sort.get_vid_by_vix(self.vix)

        if auto_perfVariants is False:

            #deterine new subs + sups (based on conclusions)
            kpc_ = sorted(list(set(self.kpc + pos + pos_a)))
            sups_ = self.get_rel(relType=1,kpc=kpc_,allowImperfect=allowImperfect) #not from cache
            subs_ = self.get_rel(relType=-1,kpc=kpc_,allowImperfect=allowImperfect) #not from cache

            #prep data for mx_sort_recommendations
            if len(pos):
                print("pos: %s [%s]" % (l2s(self.sort.get_kname_by_kix(pos)),l2s(pos)))
                rec.update({"p":pos})
            if len(pos_a):
                print("pos_a: %s [%s]" % (l2s(self.sort.get_kname_by_kix(pos_a)),l2s(pos_a)))
                rec.update({"pa":pos_a})
            if len(neg):
                print("neg: %s [%s]" %  (l2s(self.sort.get_kname_by_kix(neg)),l2s(neg)))
                rec.update({"n":neg})
            if len(sups_):
                print("sups: %s [%s]" % (l2s(self.sort.get_vname_by_vix(sups_)),l2s(sups_)))
                rec.update({"sups":sups_})
            if len(subs_):
                print("subs: %s [%s]" % (l2s(self.sort.get_vname_by_vix(subs_)),l2s(subs_)))
                rec.update({"subs":subs_})

            #reset mx_sort_recommendations with new instructions
            sql = "delete from mx_sort_recommendations where vID=%s"%vid
            self.dbo.sql_exec(sql)
            if len(rec.items()):
                sql = '''insert into mx_sort_recommendations (vID,instructions) values (%s,"%s");'''%(vid,str(rec).replace(" ",""))
                self.dbo.sql_exec(sql)

        if len(spl) > 0:
            return 1 #split errors for imperfect variant chk
        else:
            return 0 #no split errors for imperfect variant chk

        if auto_mode is not True:
            print("")
            print("---------------------------------------------------------------------")
            print("")
        
    def upd_unk(self,vname=None,auto_nonsplits=False):
        if auto_nonsplits is False:
            self.sort.restore_mx_data()

        #Note: arg input handling (not for auto_nonsplits)
        if auto_nonsplits is False:
            self.vix = self.proc_vname(vname)
            if self.vix is None:
                print ("No variants could be found. Exiting.")
                sys.exit()

        #Note: retrieve insructions for this variant (or exit, if nothing)
        vid = self.sort.get_vid_by_vix(self.vix)
        sql = "select instructions from mx_sort_recommendations where vID=%s"%vid
        self.dbo.sql_exec(sql)
        row = self.dbo.fetchone()
        if row is None:
            print("vix (%s [%s]) - there are no saved instructions.  Exiting."%(self.vix,self.sort.get_vname_by_vix(self.vix)))
            if auto_nonsplits is False:
                sys.exit()
            else:
                return 0 #TODO: need to flag these situations
        rec = row[0].replace("'",'"')
        recD = json.loads(rec)
        if auto_nonsplits is False:
            print("")

        #Note: split inconsistency found - shouldn't be getting here, but just in case...
        if auto_nonsplits and 'spl' in recD.keys():
            print("split found! Not going to do anything.")
            return 1 #TODO: need to set up processes to deal with splits

        #Note: positive calls
        if 'p' in recD.keys():
            for kix in recD['p']:
                if auto_nonsplits is False:
                    print("Changing: coord (%s:%s) to POS"%(self.sort.get_vname_by_vix(self.vix),self.sort.get_kname_by_kix(kix)))
                kid = self.sort.get_kid_by_kix(kix)
                sql = "select * from mx_calls where vID=%s and pID=%s"%(vid,kid)
                self.dbo.sql_exec(sql)
                row = self.dbo.fetchone()
                if row is not None:
                    oldval = row[4]
                    sql = "update mx_calls set assigned=1,confidence=2,changed=%s where pID=%s and vID=%s;"%(oldval,kid,vid)
                    self.dbo.sql_exec(sql)
                else:
                    oldval = self.sort.NP[self.vix,kix]
                    sql = "insert into mx_calls (pID,vID,assigned,confidence,changed) values (%s,%s,%s,%s,%s);"%(kid,vid,1,2,oldval)
                    self.dbo.sql_exec(sql)
                self.sort.NP[self.vix,kix] = 1

        #Note: positive ambiguous calls
        if 'pa' in recD.keys():
            for kix in recD['pa']:
                if auto_nonsplits is False:
                    print("Changing: coord (%s:%s) to POS-A"%(self.sort.get_vname_by_vix(self.vix),self.sort.get_kname_by_kix(kix)))
                self.sort.NP[self.vix,kix] = 1
                kid = self.sort.get_kid_by_kix(kix)
                sql = "select * from mx_calls where vID=%s and pID=%s"%(vid,kid)
                self.dbo.sql_exec(sql)
                row = self.dbo.fetchone()
                if row is not None:
                    oldval = row[4]
                    sql = "update mx_calls set assigned=1,confidence=1,changed=%s where pID=%s and vID=%s;"%(oldval,kid,vid)
                    self.dbo.sql_exec(sql)
                else:
                    oldval = self.sort.NP[self.vix,kix]
                    sql = "insert into mx_calls (pID,vID,assigned,confidence,changed) values (%s,%s,%s,%s,%s);"%(kid,vid,1,1,oldval)
                    self.dbo.sql_exec(sql)
                self.sort.NP[self.vix,kix] = 1

        #Note: negative calls
        if 'n' in recD.keys():
            for kix in recD['n']:
                if auto_nonsplits is False:
                    print("Changing: coord (%s:%s) to NEG"%(self.sort.get_vname_by_vix(self.vix),self.sort.get_kname_by_kix(kix)))
                kid = self.sort.get_kid_by_kix(kix)
                sql = "select * from mx_calls where vID=%s and pID=%s"%(vid,kid)
                self.dbo.sql_exec(sql)
                row = self.dbo.fetchone()
                if row is not None:
                    oldval = row[4]
                    sql = "update mx_calls set assigned=-1,confidence=2,changed=%s where pID=%s and vID=%s;"%(oldval,kid,vid)
                    self.dbo.sql_exec(sql)
                else:
                    oldval = self.sort.NP[self.vix,kix]
                    sql = "insert into mx_calls (pID,vID,assigned,confidence,changed) values (%s,%s,%s,%s,%s);"%(kid,vid,-1,2,oldval)
                    self.dbo.sql_exec(sql)
                self.sort.NP[self.vix,kix] = -1

        #Note: vars needed for not used chk
        kuc = self.sort.get_kixs_by_val(val=0,vix=self.vix)
        kpc = self.sort.get_kixs_by_val(val=1,vix=self.vix)
        knc = self.sort.get_kixs_by_val(val=-1,vix=self.vix)
        if auto_nonsplits is False:
            print("")
        if len(kuc) != 0:
            if 1==1 and auto_nonsplits is False:
                print("There was a problem. Exiting.")
                print("vix (%s [%s]) - there was a problem.  Exiting."%(self.vix,self.sort.get_vname_by_vix(self.vix)))
                sys.exit()

        #Note: updated variant reveals it's not usable
        if len(kpc) == 0 or len(knc) == 0:
            if auto_nonsplits is False:
                if len(kpc) == 0 : print("Updated unks for this variant reveals no positives. Achiving.")
                if len(knc) == 0 : print("Updated unks for this variant reveals no negatives. Achiving.")

            #not used variants
            if len(knc) == 0:
                sql1 = "update mx_calls set removal=1 where vid=%s;"%vid
                sql2 = "insert into mx_notused_variants(vId,reasonId) values(%s,1);"%vid
            elif len(kpc) == 0:
                sql1 = "update mx_calls set removal=2 where vid=%s;"%vid
                sql2 = "insert into mx_notused_variants(vId,reasonId) values(%s,2);"%vid
            self.dbo.sql_exec(sql1)
            self.dbo.sql_exec(sql2)

            #remove variant from self.VARIANTS + matrix
            self.sort.mx_remove_vix(self.vix)

        #Note: it is usable
        else:
            #if not a dupe and not a notused variant, then ...
            vid = self.sort.get_vid_by_vix(self.vix)

            #Note: drop any prior sup/sub definitions for this variant
            sql = "delete from mx_sups_subs where sup = %s or sub = %s;"%(vid,vid)
            self.dbo.sql_exec(sql)
            cnt_new_sups_subs = 0

            #new sups
            if 'sups' in recD.keys():
                sups_ = []
                for sup_i in self.sort.get_vid_by_vix(recD['sups']):
                    sups_.append((sup_i,vid))
                #insert new sup info
                sql = "insert into mx_sups_subs (sup,sub) values (?,?);"
                self.dbo.sql_exec_many(sql,sups_)
                cnt_new_sups_subs = cnt_new_sups_subs + len(recD['sups'])

            #new subs
            if 'subs' in recD.keys():
                subs_ = []
                for sub_i in self.sort.get_vid_by_vix(recD['subs']):
                    subs_.append((vid,sub_i))
                #insert new sub info
                sql = "insert into mx_sups_subs (sup,sub) values (?,?);"
                self.dbo.sql_exec_many(sql,subs_)
                cnt_new_sups_subs = cnt_new_sups_subs + len(recD['subs'])

            #Note: recreate panda sups/subs cache
            sql = "select distinct * from mx_sups_subs;"
            self.dbo.sql_exec(sql)
            sups_ = self.dbo.fetchall()
            self.SS = pd.DataFrame(sups_,columns=['sup','sub'])

            #Note: stdout about what's new
            if cnt_new_sups_subs > 0 and auto_nonsplits is False:
                print("Found %s new sup-sub relations" % cnt_new_sups_subs)

        if auto_nonsplits is False:
            self.sort.mx_remove_dupes()
            self.sort.mx_vandh_sort()
            self.sort.save_mx()
            print("\nTotal matrix kits: %s" % len(self.sort.KITS))
            print("Total perfect matrix variants: %s" % len(self.sort.get_perfect_variants_idx()))
            print("Total imperfect matrix variants: %s" % len(self.sort.get_imperfect_variants_idx()))
            print("Total matrix variants: %s\n" % len(self.sort.VARIANTS))

    def get_rel(self,relType,override_val=None,kix=None,kpc=None,allowImperfect=False):
        #-------------------------------------------------------------------------------------------------{{{
        # 3 Ways to use this routine
        #-------------------------------------------------------------------------------------------------
        # 1. can use it with an override_val and a kix -- in this situation, it will override a given kit 
        #    (ie: one with an unkown val, to a positive val) and then determine the subsets of this 
        #    variant -- again with a routine generated kpc
        # 2. can use it with a vix (it gets the default subs that way), the routine creates a kpc
        # 3. can use it with a given kpc to override the kpc generated by this routine if you use methods 
        #    1 or 2. this is useful if you want to check subsets when overriding multiple kit unknown 
        #    values for a given variant
        #-------------------------------------------------------------------------------------------------
        # 3 Types of relations
        #-------------------------------------------------------------------------------------------------
        # 1. relType ==  1 supset
        # 2. relType == -1 subset
        # 3. relType ==  0 equivalent
        #-------------------------------------------------------------------------------------------------}}}

        #method 1
        if override_val is not None and kix is not None:
            if config['DBG_RELS']: print("\n[rels.0] MODE (get_vix_rel) - ONE (override val w/vix)\n")
            overrideData = self.sort.get_row_when_override_coord(override_val,kix=kix,vix=self.vix)
            #pos conditions when override coord with a value
            kpc = self.sort.get_kixs_by_val(1,vix=self.vix,overrideData=overrideData)

        #method 2
        elif kpc is None:
            if config['DBG_RELS']: print("\n[rels.0] MODE (get_vix_rel) - TWO (only use vix)\n")
            #default pos conditions when just use default target vix criteria
            kpc = self.sort.get_kixs_by_val(1,vix=self.vix) #default pos conditions
    
        #method3
        #requires only sending in a kpc
        else:
            if config['DBG_RELS']: print("\n[rels.0] MODE (get_vix_rel) - THREE (sending in a kpc)\n")

        if relType == 1: #supsets
            rels_ = self.get_supsets_or_eqv(kpc,allowImperfect=allowImperfect)
        elif relType == -1: #subsets
            rels_ = self.get_subsets(kpc,allowImperfect=allowImperfect)
        else: # relType == 0: subsets (eq_override flag)
            #Note: this creates data that needs to be compared to supsets, to
            #get equal variants ... this is def a hack!
            rels_ = self.get_supsets_or_eqv(kpc=kpc,eq_override=True,allowImperfect=allowImperfect)
        
        #if there's nothing here, best to just return
        if len(rels_) == 0:
            return []
        #two ways to return the data. with or without counts.
        #rels = self.filter_perfect_variants(vix=rel_.tolist()).tolist()
        if allowImperfect is False:
            rels = self.sort.filter_perfect_variants(vix=rels_[:,0].tolist()) #.tolist()
        else:
            rels = rels_[:,0].tolist()

        #debugging
        if config['DBG_RELS']: print("[rels.1] rels: %s" % rels)
        if config['DBG_RELS']:
            suffix=''
            if override_val is not None and kix is not None:
                suffix = 'P' if override_val == 1 else 'N'
            print("[rels.2] rels%s: %s kpc: %s" % (suffix,",".join([str(i) for i in rels]),kpc))

        return list(rels)
        
    def get_subsets(self,kpc,allowImperfect=False):

        #imperfect toggle
        if allowImperfect:
            NP = self.sort.NP
        else:
            NP = self.sort.NP[self.sort.get_perfect_variants_idx()]

        if config['DBG_SUBS_IN']:
            print("\n---------------------------------------------------------------------\n")
            print("[subin.0] vix: %s"%self.vix)
            print("[subin.0] vix(name): %s"%self.sort.get_vname_by_vix(self.vix))
            print("[subin.0] kpc: %s"%kpc)
            print("[subin.0] kpc(names): %s"%self.sort.get_kname_by_kix(kpc))

        #looking for variants w/pos assignments that overlap the target variant
        VAR1 = np.argwhere(NP[:,kpc]==1)[:,0]
        if config['DBG_SUBS_IN']: print("[subin.1] VAR1: %s"%VAR1)

        #idx make sure we exclude the target variant
        idx = np.argwhere(VAR1==self.vix)
        VAR2 = np.delete(VAR1, idx)
        if config['DBG_SUBS_IN']: print("[subin.2] VAR2: %s"%VAR2)

        #get uniques of VAR2 with counts
        unique_elements, counts_elements = np.unique(VAR2, return_counts=True)
        VAR3 = np.asarray((unique_elements, counts_elements)).T #get uniques
        if config['DBG_SUBS_IN']: print("\n[subin.3] VAR3: %s"%VAR3)

        #for the following "adding technique" -- we need to exclude unk situations for the comparison (special handling)
        #adding technique - got this idea here: https://stackoverflow.com/questions/30041286/sum-rows-where-value-equal-in-column
        unq, unq_inv = np.unique(VAR3[:,0], return_inverse=True)
        out = np.zeros((len(unq), VAR3.shape[1]), dtype=VAR3.dtype) #create empty array to put the added values
        out[:, 0] = unq #fill the first column
        np.add.at(out[:, 1:], unq_inv, VAR3[:, 1:])
        if config['DBG_SUBS_IN']: print("[subin.4] out: %s"%out)

        #sorting without fields - https://stackoverflow.com/questions/2828059/sorting-arrays-in-numpy-by-column
        out1 = out[out[:,1].argsort()[::-1]] #reverse sort (2nd col) -- to get the max ones first
        if config['DBG_SUBS_IN']:
            print("[subin.4a] out1: %s"%out1)
            print("[subin.5] kpc: %s"%kpc)
            print("[subin.6] kpc(names): %s"%self.sort.get_kname_by_kix(kpc))

        #subsets need to have less kpc than the target variant (cuz subsets)
        VAR4 = np.argwhere(out1[:,1]<len(kpc)) #[:,0]
        if config['DBG_SUBS_IN']: print("[subin.7] VAR4: %s"%VAR4)

        #separate the cols - for debugging (this is the subset variants)
        out2a = out1[:,0]
        if config['DBG_SUBS_IN']: print("[subin.8] out2a: %s"%out2a)

        #and this is their kpc count data
        out2b = out1[:,1]
        if config['DBG_SUBS_IN']: print("[subin.9] out2b: %s"%out2b)

        #the subset variants
        VAR5 = out2a[list(VAR4.T[0])]
        if config['DBG_SUBS_IN']:
            print("[subin.10] VAR5: %s"%VAR5)
            print("[subin.11] VAR5(names): %s"%self.sort.get_vname_by_vix(VAR5))

        #the subset variant counts
        VAR6 = out2b[list(VAR4.T[0])]
        if config['DBG_SUBS_IN']: print("[subin.12] VAR6: %s"%VAR6)

        #subsets shouldn't have diff intersecting positives with common supersets
        subsc = list(VAR5)
        invalid_subs = []
        self.sups = self.get_rel(relType=1,kpc=kpc,allowImperfect=False)
        if len(self.sups) == 0:
            self.sups.append(-999) #top handling
        #self.kpc = self.sort.get_kixs_by_val(val=1,vix=self.vix)
        self.kpc = kpc
        for v in subsc:
            vid = self.sort.get_vid_by_vix(v)
            try:
                v1sups = self.sort.get_vixs_by_vids(list(self.sort.SS.loc[self.sort.SS['sub']==vid,['sup']].as_matrix().T[0]))
            except:
                v1sups = []
            valid_sub = 1
            if len(v1sups) == 0:
                v1sups.append(-999) #top handling
            if len(v1sups):
                common_sups = set(v1sups).intersection(set(self.sups))
                if config['DBG_SUBS_IN']:
                    print("[subin.10a] self_sups: %s"%self.sups)
                    print("[subin.10b] v1sups: %s"%v1sups)
                    print("[subin.10c] common_sups: %s"%common_sups)
                v1kpc = self.sort.get_kixs_by_val(val=1,vix=v)
                for sup in common_sups:
                    v2 = Variant()
                    v2.dbo = self.dbo
                    v2.sort = self.sort
                    v2.vix = sup
                    v2kpc = self.sort.get_kixs_by_val(val=1,vix=sup)
                    common_kpc1 = set(self.kpc).intersection(set(v2kpc))
                    common_kpc2 = set(v1kpc).intersection(set(v2kpc))
                    for c in common_kpc2:
                        if c not in common_kpc1:
                            valid_sub = 0
                            invalid_subs.append(v)
                            break 
                    if valid_sub == 0 : break
        if len(invalid_subs):
            for v in invalid_subs:
                subsc.remove(v)
        search = np.array(invalid_subs)
        #(beg) technique found here: https://stackoverflow.com/questions/32191029/getting-the-indices-of-several-elements-in-a-numpy-array-at-once/32191125#32191125
        sorter = np.argsort(VAR5)
        idx = sorter[np.searchsorted(VAR5,invalid_subs,sorter=sorter)]
        #(end)
        VAR5a = np.delete(VAR5, idx)
        VAR6a = np.delete(VAR6, idx)

        #merged for return
        VAR7 = np.asarray((VAR5a,VAR6a)).T

        if config['DBG_SUBS_IN']: print("[subin.13] VAR7: %s"%VAR7)

        return VAR7
        
    def get_supsets_or_eqv(self,kpc,eq_override=False,allowImperfect=False):

        #imperfect toggle
        if allowImperfect:
            NP = self.sort.NP
        else:
            NP = self.sort.NP[self.sort.get_perfect_variants_idx()]

        if len(kpc) == 0: #if there's no positive kpc to compare
            return [] #then return nothing

        if config['DBG_SUPS_IN']:
            print("\n---------------------------------------------------------------------\n")
            print("[supin.0] vix: %s"%self.vix)
            print("[supin.0] vix(name): %s"%self.sort.get_vname_by_vix(self.vix))
            print("[supin.0] kpc: %s"%kpc)
            print("[supin.0] kpc(names): %s"%self.sort.get_kname_by_kix(kpc))

        #look for variants that intersect target variant along its kpc (exclude target
        #variant from results) -- these are the superset candidates
        vixIntersects = np.argwhere(np.all(NP[:,kpc]==[1]*len(kpc),axis=1)==True)[:,0]
        if len(vixIntersects) == 0: return [] #if there are no kpc intersects, end here
        if config['DBG_SUPS_IN']: print("[supin.1] vixIntersects: %s"%vixIntersects)

        #get the kpc for each of these superset candidates
        allPosVix = np.argwhere(NP[vixIntersects]==1)[:,0]
        unqA, cntA = np.unique(allPosVix, return_counts=True)
        IC = np.asarray((vixIntersects,cntA)).T #IC = intersect counts
        if config['DBG_SUPS_IN']: print("[supin.2] intersectCnts: %s"%IC)

        #now check to see if these superset candidate vix have more/equal kpc to target variant
        for itm in reversed([(idx,val) for idx,val in enumerate(vixIntersects.tolist())]):
            #this logic drops equivalents (eq_override arg)
            if eq_override is False and itm[1] == IC[itm[0]][0] and IC[itm[0]][1] == len(kpc):
                IC = np.delete(IC,(itm[0]),axis=0)
            #this logic drops everything but equivalents
            elif eq_override is True and itm[1] == IC[itm[0]][0] and IC[itm[0]][1] != len(kpc):
                IC = np.delete(IC,(itm[0]),axis=0)

        #return filtered results ... includes counts (in case they're needed)
        if config['DBG_SUPS_IN']: print("[supin.3] intersectCnts: %s"%IC)
        return IC

class Sort(object):

    def __init__(self):
        self.dbo = None #db object
        self.KITS = None
        self.VARIANTS = None
        self.NP = None #matrix
        self.SS = None #sups+subs

    # stdout

    def stdout_matrix(self,vix=None,kix=None):
        if vix is not None and kix is not None:
            vix = sorted(vix)
            kix = sorted(kix)
            NP = self.NP[vix]
            NP = NP[:,kix]
            VARIANTS = self.get_vix_by_vix(vix)
            KITS = self.get_kname_by_kix(kix)
        elif kix is not None:
            kix = sorted(kix)
            NP = self.NP[:,kix]
            VARIANTS = self.get_axis('variants')
            KITS = self.get_kname_by_kix(kix)
        elif vix is not None:
            vix = sorted(vix)
            NP = self.NP[vix]
            VARIANTS = self.get_vix_by_vix(vix)
            KITS = self.get_axis('kits',keysOnly=True)
        else:
            NP = self.NP
            VARIANTS = self.get_axis('variants')
            KITS = self.get_axis('kits',keysOnly=True)

        print("beg stdout-matrix: %s" % format(time.clock()))
        print("")

        #matrix size chk
        if config['SHOW_STDOUT_MATRIX_OVER_50_ROWS'] is False and len(NP) > 50:
            print("Matrix too large for presentation.")
            np.set_printoptions(threshold=np.inf)
            print(NP)
            return

        #debugging
        if config['DBG_MATRIX']:
            print("---------------------------------------")
            print("Matrix: bef stdout")
            print("---------------------------------------")
            print(NP)
            print("")

        #col nums
        if kix:
            lenCols = len(kix)
            kixCols = kix
        else:
            lenCols = len(self.KITS)
            kixCols = [str(x) for x in list(range(lenCols))]

        #matrix
        table = BeautifulTable(max_width=155)
        table.column_headers = ['c']+['p']+['v']+[str(x) for x in KITS]
        table.append_row(['']+['']+['']+kixCols)
        table.append_row(['']+['']+['']+['']*lenCols)
        cntV = 0
        for K,V in VARIANTS:
            table.append_row([self.get_vix_by_name(K)]+[str(V[2])]+[K]+[str(x).replace('-1','') for x in NP[cntV,:].tolist()[0]])
            table.row_seperator_char = ''
            table.column_seperator_char = ''
            table.column_alignments['v'] = BeautifulTable.ALIGN_LEFT
            cntV = cntV + 1
        print(table)

        #Note: counts for stdout
        print("\nTotal matrix kits: %s" % len(self.KITS))
        print("Total perfect matrix variants: %s" % len(self.get_perfect_variants_idx()))
        print("Total imperfect matrix variants: %s" % len(self.get_imperfect_variants_idx()))
        print("Total matrix variants: %s\n" % len(self.VARIANTS))

        #done
        print("")
        print("end stdout-matrix: %s" % format(time.clock()))
        
    def stdout_unknowns(self):
        self.dbo.db = self.dbo.db_init()
        self.dbo.dc = self.dbo.cursor()
        self.restore_mx_data()
        print("")
        print("---------------------------------------------------------------------")
        print("")
        cnt = 1
        if len(self.get_imperfect_variants_idx()) == 0:
            print("There are no imperfect variants in this matrix.")
            print("")

        for vix in self.get_imperfect_variants_idx():
            name = self.get_vname_by_vix(vix)
            kuc = self.get_kixs_by_val(val=0,vix=vix)
            kucn = "kuc: %s [%s]"%(l2s(self.get_kname_by_kix(kuc)),l2s(kuc))
            print("[%s]  var: %s" %(vix,name))
            print("     %s" %kucn)
            print("")
            cnt = cnt + 1
        print("---------------------------------------------------------------------")
        print("")

    # matrix 

    def sort_matrix(self):
        #db
        self.dbo.db = self.dbo.db_init()
        self.dbo.dc = self.dbo.cursor()

        #get data
        self.create_mx_data()

        #sort
        self.mx_vandh_sort()

        #cache sups + subs
        self.reset_ss_data()

        #save data
        self.save_mx()

    def get_vix_by_vix(self,vix):
        newVix = []
        for vo in list(vix):
            for itm in list(self.VARIANTS.items()):
                if itm[1][1] == vo:
                    newVix.append(itm)
                    break
        return(newVix)
        
    def get_kix_by_kix(self,kix):
        newKix = []
        for ko in list(kix):
            for itm in list(self.KITS.items()):
                if itm[1][1] == ko:
                    newKix.append(itm)
                    break
        return(newKix)
        
    def get_vid_by_vix(self,vix):
        intFlg = True
        try:
            value = int(vix)
        except:
            intFlg = False
        if intFlg:
            vix = [vix]
        vidList = []
        for vo in vix:
            for itm in list(self.VARIANTS.items()):
                if itm[1][1] == vo:
                    vidList.append(itm[1][0])
                    break
        if intFlg:
            return vidList[0]
        else:
            return vidList
        
    def get_kid_by_kix(self,kix):
        intFlg = True
        try:
            value = int(kix)
        except:
            intFlg = False
        if intFlg:
            kix = [kix]
        kidList = []
        for ko in kix:
            for itm in list(self.KITS.items()):
                if itm[1][1] == ko:
                    kidList.append(itm[1][0])
                    break
        if intFlg:
            return kidList[0]
        else:
            return kidList
        
    def get_vname_by_vix(self,vix):
        intFlg = True
        try:
            value = int(vix)
        except:
            intFlg = False
        if intFlg:
            vix = [vix]
        vnList = []
        for vo in vix:
            if vo == -999: #top handling
                vnList.append('top')
            else:
                for itm in list(self.VARIANTS.items()):
                    if itm[1][1] == vo:
                        vnList.append(itm[0])
                        break
        if intFlg:
            return vnList[0]
        else:
            return vnList
        
    def get_kname_by_kix(self,kix):
        intFlg = True
        try:
            value = int(kix)
        except:
            intFlg = False
        if intFlg:
            kix = [kix]
        knList = []
        for ko in kix:
            for itm in list(self.KITS.items()):
                if itm[1][1] == ko:
                    knList.append(itm[0])
                    break
        if intFlg:
            return knList[0]
        else:
            return knList

    def get_vix_by_name(self,vnames):
        vixs = []
        strFlg = isinstance(vnames, str)
        if strFlg:
            vnames = [vnames]
        for vn in vnames:
            if vn.upper() in self.VARIANTS.keys():
                vixs.append(self.VARIANTS[vn.upper()][1])
        if strFlg:
            if len(vixs):
                return vixs[0]
            else:
                return None
        else:
            return vixs
        
    def get_kix_by_name(self,knames):
        kixs = []
        strFlg = isinstance(knames, str)
        if strFlg:
            knames = [knames]
        for kn in knames:
            kixs.append(self.KITS[kn][1])
        if strFlg:
            return kixs[0]
        else:
           return kixs
        
    def get_kid_by_name(self,kname):
        return self.KITS[kname][0]
        
    def get_vid_by_name(self,vname):
        return self.VARIANTS[vname][0]

    def get_vixs_by_vids(self,vids):
        intFlg = True
        try:
            value = int(vids)
        except:
            intFlg = False
        if intFlg:
            vids = [vids]
        vixs = []
        for v in vids:
            for itm in list(self.VARIANTS.items()):
                if itm[1][0] == v:
                    vixs.append(itm[1][1])
                    break
        if intFlg:
            return vixs[0]
        else:
            return vixs
    def get_kixs_by_kids(self,kids):
        intFlg = True
        try:
            value = int(kids)
        except:
            intFlg = False
        if intFlg:
            kids = [kids]
        kixs = []
        for k in kids:
            for itm in list(self.KITS.items()):
                if itm[1][0] == k:
                    kixs.append(itm[1][1])
                    break
        if intFlg:
            return kixs[0]
        else:
            return kixs

    def get_kixs_by_val(self,val,vix=None,vname=None,overrideData=None):
        if vname is not None:
            vix = self.get_vix_by_name(vname)
        if vix is not None and overrideData is not None: # we're sending in a custom evaluation
            return list(np.argwhere(overrideData[0,] == val).T[1,]) #with override data, there's only one line evaluated - 1d datset
        elif vix is not None: #no override -- use self.NP (all data)
            if vix == -999: #top handling
                if val == 1:
                    return list(range(len(self.KITS)))
                if val in [-1,0]:
                    return []
            return list(np.argwhere(self.NP[vix,] == val).T[1,]) #default data, it's the entire matrix - 2d dataset 
        
    def get_vixs_by_val(self,val,kix=None,kname=None,overrideData=None):
        if kname is not None:
            kix = self.get_kix_by_name(kname)
        if kix is not None and overrideData is not None:
            return np.argwhere(overrideData[:,0] == val).T[0,] #with override data, there's only one line evaluated - 1d dataset
        elif kix is not None: #no override -- use self.NP (all data)
            return np.argwhere(self.NP[:,kix] == val).T[0,] #default data, it's the entire matrix - 2d dataset
        
    def get_knames_by_val(self,val,vix=None,vname=None,overrideData=None):
        return self.get_kname_by_kix(self.get_kix_by_val(val,vix,vname,overrideData))
        
    def get_vnames_by_val(self,val,kix=None,kname=None,overrideData=None):
        return self.get_vname_by_vix(self.get_vixs_by_val(val,kix,kname,overrideData))

    def mx_vandh_sort(self):
        #Note: vertical/horizontal sort of the matrix

        def return_counts(idx,inv):
            count = np.zeros(len(idx), np.int)
            np.add.at(count, inv, 1)
            return count

        #Note: debugging line outs
        if config['DBG_MATRIX']:
            print("\n=================================")
            print("inside vandh sort")
            print("=================================\n")
            print("\n---------------------------------------")
            print("Matrix: beg vandh sort")
            print("---------------------------------------")
            print(self.NP)
            print("")

        #Note: perfect variants
        prfVix = self.get_perfect_variants_idx()

        #Note: needs a certain amount of perfect vix (or it breaks)
        if len(prfVix) < 2:
            print("Only %s perfect vix found. Not enough. Exiting.\n" % len(prfVix))
            sys.exit()

        #Note: panda special handling to coord perfect vix + kix sorting
        C = pd.DataFrame(np.argwhere(self.NP[prfVix]==1),columns=['vix2','kix'])
        VX = pd.DataFrame(prfVix,columns=['vix1'])
        VX['vix2'] = range(len(VX))
        VC = pd.DataFrame(return_counts(list(range(len(prfVix))),np.argwhere(self.NP[prfVix]==1)[:,0]),columns=['cntV'])
        KC = pd.DataFrame(return_counts(list(range(len(self.KITS))),np.argwhere(self.NP[prfVix]==1)[:,1]),columns=['cntK'])
        KC['kix'] = range(len(KC))
        VC2 = pd.concat([VX,VC], axis=1, join_axes=[VC.index])
        C_VX = pd.merge(pd.merge(C,VC2, on='vix2'),KC,on='kix').sort_values(['cntV','cntK'], ascending=[False, False]).reset_index(drop=True)
        prfVixSorted = list(C_VX.groupby('vix1',as_index=False).head(1)['vix1'].as_matrix())
        KixSorted = list(C_VX.groupby('kix',as_index=False).head(1)['kix'].as_matrix())

        #Note: imperfect variants
        impVix = self.get_imperfect_variants_idx()
        impVixPosV = np.argwhere(self.NP[impVix]==1)[:,0]
        unqIV, cntIV = np.unique(impVixPosV, return_counts=True)
        allIV = np.asarray((impVix,cntIV)) #idx w/counts
        allIV = allIV[:,np.argsort(-allIV[1])] #vsort
        impVixSorted = allIV.T[:,0] #sorted vix

        #Note: vsort
        VixSorted = np.concatenate((prfVixSorted,impVixSorted))
        self.NP = self.NP[np.concatenate((prfVixSorted,impVixSorted)),]

        #Note: if there are any kits w/o positives, figure it out here
        cnts = return_counts(list(range(len(self.KITS))),KixSorted)
        if 0 in list(cnts):
            KixSorted = KixSorted + [i for i in list(range(len(self.KITS))) if cnts[i] == 0]

        #re-map variant names with new indexes
        vnamesL = self.get_vname_by_vix(VixSorted)
        #re-map kit names with new indexes
        knamesL = self.get_kname_by_kix(KixSorted)

        #new matrix
        if config['DBG_MATRIX']:
            print("---------------------------------------")
            print("Matrix: vandh - aft vsort")
            print("---------------------------------------")
            print(self.NP)

        #hsort
        self.NP = self.NP[:,KixSorted]
        if config['DBG_MATRIX']:
            print("\n---------------------------------------")
            print("Matrix: vandh - aft hsort")
            print("---------------------------------------")
            print(self.NP)
            print("")

        #variants
        for ix,v in enumerate(vnamesL):
            #self.VARIANTS[v][1] = ix
            self.VARIANTS[v] = (self.VARIANTS[v][0],ix,self.VARIANTS[v][2])

        #kits
        for ix,k in enumerate(knamesL):
            #self.KITS[k][1] = ix
            self.KITS[k] = (self.KITS[k][0],ix)

    def get_axis(self,orderByType=None,keysOnly=False):
        #Note: this is a useful function for being able to enumerate the VARIANTS/KITS dictionaries

        if orderByType in ['variants','kits']:
            if orderByType == 'variants' : SCH = self.VARIANTS
            if orderByType == 'kits' : SCH = self.KITS
            if keysOnly:
                return [i[0] for i in sorted(SCH.items(), key=lambda e: e[1][1])]
            else:
                return sorted(SCH.items(), key=lambda e: e[1][1])

    def get_mx_row_as_list(self,rownum,NP=None):
        #TODO: this NP override is a hack, need to work out why I need this

        if NP is not None:
            NP = self.NP
        else:
            NP = NP
        return NP[rownum,:].tolist()[0]
        
    def get_mx_kit_data(self,vix=None,vname=None):
        if vname is not None:
            vix = self.get_vix_by_name(vname)
        if vix is not None:
            return self.NP[vix,]
        
    def get_mx_variant_data(self,kix=None,kname=None):
        if kname is not None:
            kix = self.get_kix_by_name(kname)
        if kix is not None:
            return self.NP[:,kix].T

    def get_perfect_variants_idx(self):
        idx = list(range(len(self.VARIANTS)))
        prf_idx = idx[:]
        unk_idx = list(np.unique(np.argwhere(self.NP==0)[:,0]))
        for x in idx:
            if x in unk_idx:
                prf_idx.remove(x)
        return prf_idx
        
    def get_imperfect_variants_idx(self):
        return list(np.unique(np.argwhere(self.NP==0)[:,0])) #make it a copy
        
    def filter_perfect_variants(self,vix):
        if any(isinstance(el, list) for el in vix):
            for itm in reversed([[n,v] for (n,(v,c)) in enumerate(vix)]):
                if itm[1] in self.get_imperfect_variants_idx():
                    vix.remove(vix[itm[0]])
            return vix
        else:
            for itm in reversed([[n,v] for n,v in enumerate(vix)]):
                if itm[1] in self.get_imperfect_variants_idx():
                    vix.remove(vix[itm[0]])
            return vix

    def get_row_when_override_kixs(self,override_val,vix,kixs):
        row = self.get_mx_kit_data(vix=vix)
        if len(kixs) == 0:
            return row
        rowO = np.empty_like(row)
        rowO[:] = row #make duplicate copy - important!
        for kx in kixs:
            rowO[0,kx] = override_val
        return rowO
        
    def get_row_when_override_coord(self,override_val,kix=None,vix=None,kname=None,vname=None):
        row = self.get_mx_kit_data(vix=vix)
        rowO = np.empty_like(row)
        rowO[:] = row #make duplicate copy - important!
        rowO[0,kix] = override_val
        return rowO

    def mx_remove_vix(self,vix):
        vid_ = self.get_vid_by_vix(vix)
        vname = self.get_vname_by_vix(vix)

        #TODO: need to check if there are other dupes to represent the missing variant. 

        #remove variant from self.VARIANTS
        idx_ = list(range(len(self.VARIANTS)))
        idx_.pop(vix)
        self.VARIANTS.pop(vname,None)

        #reset the matrix idxs
        for ix,vn in enumerate(self.get_vname_by_vix(idx_)):
            self.VARIANTS[vn] = (self.VARIANTS[vn][0],ix,self.VARIANTS[vn][2])

        #update mx_variants table
        sql = "delete from mx_variants;"
        self.dbo.sql_exec(sql)
        sql = "insert into mx_variants (ID,name,pos) values (?,?,?);"
        self.dbo.sql_exec_many(sql,[(tuple([vid,nm,pos])) for (n,(nm,(vid,idx,pos))) in enumerate(self.get_axis('variants'))])

        #remove vid from mx_sups_subs
        sql = "delete from mx_sups_subs where sup = %s or sub = %s" % (vid_,vid_)
        self.dbo.sql_exec(sql)

        #update panda version of mx_sups_subs
        sql = "select distinct * from mx_sups_subs;"
        self.dbo.sql_exec(sql)
        sups_ = self.dbo.fetchall()
        self.SS = pd.DataFrame(sups_,columns=['sup','sub'])

        #remove vid from mx_dupe_variants
        sql = "delete from mx_dupe_variants where vID = %s or dupe_vID = %s" % (vid_,vid_)
        self.dbo.sql_exec(sql)

        #reset the matrix
        self.NP = self.NP[idx_,]
        self.mx_vandh_sort()
        self.save_mx()

        #reset mx_sort_recommendations
        sql = "delete from mx_sort_recommendations"
        self.dbo.sql_exec(sql)

    # data 

    def sort_schema(self):
        #Note: where the sort DDL schema is first run
        self.dbo.db = self.dbo.db_init()
        self.dbo.dc = self.dbo.cursor()
        self.dbo.sql_exec_file('sort-schema.sql')
        
    def save_mx(self):
        #Note: push kit/variant/numpy data into saved/matrix tbls + h5py file

        #deletes
        sql = "delete from mx_variants;"
        self.dbo.sql_exec(sql)
        sql = "delete from mx_kits;"
        self.dbo.sql_exec(sql)
        sql = "delete from mx_idxs;"
        self.dbo.sql_exec(sql)

        #save matrix variants
        sql = "insert into mx_variants (ID,name,pos) values (?,?,?);"
        self.dbo.sql_exec_many(sql,[(tuple([vid,nm,pos])) for (n,(nm,(vid,idx,pos))) in enumerate(self.get_axis('variants'))])

        #save matrix variants (idx order)
        sql = "insert into mx_idxs (type_id,axis_id,idx) values (0,?,?);"
        self.dbo.sql_exec_many(sql,[(tuple([vid,idx])) for (n,(nm,(vid,idx,pos))) in enumerate(self.get_axis('variants'))])

        #save matrix kits
        sql = "insert into mx_kits (ID,kitId) values (?,?);"
        self.dbo.sql_exec_many(sql,[(tuple([kid,nm])) for (n,(nm,(kid,idx))) in enumerate(self.get_axis('kits'))])

        #save matrix kits (idx order)
        sql = "insert into mx_idxs (type_id,axis_id,idx) values (1,?,?);"
        self.dbo.sql_exec_many(sql,[(tuple([kid,idx])) for (n,(nm,(kid,idx))) in enumerate(self.get_axis('kits'))])

        #save numpy data
        devnull = open(os.devnull, 'w')
        from mock import patch
        with patch('sys.stdout', devnull):
            with patch('sys.stderr', devnull):
                import h5py
        h5f = h5py.File('data.h5', 'w')
        h5f.create_dataset('dataset_1', data=self.NP)
        h5f.close()
        
    def restore_mx_data(self):
        print("beg MatrixData restore: %s" % format(time.clock()))
        self.VARIANTS = {}
        self.KITS = {}

        #variants
        sql = '''
            SELECT DISTINCT V.name, V.ID, IX.idx, V.pos
            FROM mx_idxs IX
            INNER JOIN mx_variants V
            ON IX.axis_id = V.ID
            WHERE IX.type_id = 0
            ORDER BY IX.idx;
            ''';
        self.dbo.sql_exec(sql)
        F = self.dbo.fetchall()
        for row in F:
            self.VARIANTS[row[0]] = (row[1],row[2],row[3]) #self.VARIANTS[name] = [vID,idx]

        #kits
        sql = '''
            SELECT DISTINCT K.kitId, K.ID, IX.idx
            FROM mx_idxs IX
            INNER JOIN mx_kits K
            ON IX.axis_id = K.ID
            WHERE IX.type_id = 1
            ORDER BY IX.idx;
            ''';
        self.dbo.sql_exec(sql)
        F = self.dbo.fetchall()
        for row in F:
            self.KITS[row[0]] = (row[1],row[2]) #self.VARIANTS[name] = [vID,idx]

        #numpy data
        devnull = open(os.devnull, 'w')
        from mock import patch
        with patch('sys.stdout', devnull):
            with patch('sys.stderr', devnull):
                import h5py
        h5f = h5py.File('data.h5','r')
        self.NP = np.asmatrix(h5f['dataset_1'][:])
        h5f.close()

        #sups+subs
        sql = "SELECT DISTINCT * FROM mx_sups_subs;"
        self.dbo.sql_exec(sql)
        sups_ = self.dbo.fetchall()
        self.SS = pd.DataFrame(sups_,columns=['sup','sub'])
        print("end MatrixData restore: %s" % format(time.clock()))
        
    def create_mx_data(self):

        print("beg MatrixData create: %s" % format(time.clock()))

        #Note: bedranges - uses in_range routine
        #table tmp1 - the covered spots of the bedranges (for the known pos
        #call areas) are pushed to this table; we're only really interested in
        #those situations where a VCF positive call has been discovered along 
        #a position for these -- since we need at least one POS and one NEG for 
        #a variant to be placed in the matrix

        sql = "drop index if exists tmp1idx1;"
        self.dbo.sql_exec(sql)
        sql = "drop index if exists tmp1idx2;"
        self.dbo.sql_exec(sql)
        sql = "drop table if exists tmp1;"
        self.dbo.sql_exec(sql)
        sql = "create table tmp1 (pid int,vid int,pos int,val bool);"
        self.dbo.sql_exec(sql)
        sql = "create index tmp1idx1 ON tmp1(pid,vid);" #needed?
        self.dbo.sql_exec(sql)
        sql = "create index tmp1idx2 ON tmp1(pid,vid,val);" #needed?
        self.dbo.sql_exec(sql)
        sql = "select distinct C.PId FROM vcfcalls C"
        pids = self.dbo.sql_exec(sql)
        F = self.dbo.fetchall()
        for row in F:
            sql = "select * from v_pos_call_chk;"
            dc = self.dbo.cursor()
            calls = dc.execute(sql)
            rc = self.dbo.cursor()
            sql = '''select minaddr,maxaddr from bedranges r
                inner join bed b on b.bID=r.id
                where b.pid=%s order by 1'''%row[0]
            ranges = rc.execute(sql)
            xl = list([(v[0],v[1]) for v in calls])
            pv = [v[1] for v in xl]
            rv = [v for v in ranges]
            coverage = in_range(pv, rv, spans=None)
            mergedList = [(x,y,z) for (x,y),z in zip(xl, coverage)]
            sql = "insert into tmp1 (pid,vid,pos,val) values (%s,?,?,?);"%row[0]
            self.dbo.sql_exec_many(sql,mergedList)

        #Note: table tmp2 - this is the data that correlates variants pos's to whether kits were
        #tested for them or not (a neg is when it was covered, but it's not
        #showing up in the vcfcalls)
        #the pidvid idea is a hack ... so we can isolate pid/vid combos for
        #situations in the vcfcalls where a positive was uncovered, and compare
        #those results to the covered bed positions

        sql = "drop table if exists tmp2;"
        self.dbo.sql_exec(sql)
        sql = '''
            CREATE TABLE tmp2 as 
            SELECT DISTINCT vid as vID, pid as pID,
            cast(pid as varchar(15))||'|'||cast(vid as varchar(15)) as pidvid
            FROM tmp1 WHERE val=1;'''
        self.dbo.sql_exec(sql)

        #Note: anything that is already in VCF's can't be a newly discovered bed NEG
        sql = "delete from tmp2 where pidvid in (select pidvid from v_all_calls_with_kits);"
        self.dbo.sql_exec(sql)

        #Note: sql to fetch matrix variant data; this targets variants that
        #have at least one POS, one NEG across the known kits
        sql = "select * from v_imx_assignments_with_unk;"
        self.dbo.sql_exec(sql)
        F = self.dbo.fetchall()

        #Note: vars that will be used the loop coming up
        DATA = OrderedDict()
        self.KITS = {}
        self.VARIANTS = {}
        cntV = 0
        cntK = 0

        #Note: retrieve data from sqlite like so: [V][K] [x,x,x,x,x,x,...]
        for row in F:
            if row[1] not in DATA:
                DATA[row[1]] = []

            #calls
            PF = 0

            #the vcf positives
            if row[2] == 1 and row[7] == '1/1':
                PF = 1

            #the vcf negatives
            elif row[2] == 1 and row[7] == '0/0':
                PF = -1

            #the additional negs discovered by checking the beds
            elif row[6] == 1 and row[7] == None:
                PF = -1
            DATA[row[1]].append(PF)

            #kits
            if row[5] not in self.KITS.keys():
                #what this is doing: self.KITS[name] = (pID,idx)
                self.KITS[row[5]] = [row[0],cntK]
                cntK = cntK + 1

            #variants
            if row[1] not in self.VARIANTS.keys():
                #what this is doing: self.VARIANTS[name] = (vID,idx,pos)
                self.VARIANTS[row[1]] = [row[4],cntV,row[3]]
                cntV = cntV + 1

        #debugging
        self.NP = np.matrix(list(DATA.values()))
        if config['DBG_MATRIX']:
            print("\n---------------------------------------")
            print("Matrix: just created data")
            print("---------------------------------------")
            print(self.NP)

        #Note: table mx_notused_variants - these are variants that don't have at least one
        #NEG, one POS - they won't be put in the matrix
        sql = "delete from mx_notused_variants;"
        self.dbo.sql_exec(sql)
        sql = "insert into mx_notused_variants select vID,1 from v_only_pos_variants";
        self.dbo.sql_exec(sql)
        sql = "insert into mx_notused_variants select vID,-1 from v_only_neg_variants";
        self.dbo.sql_exec(sql)
        sql = "select count(*) from mx_notused_variants where reasonId = 1";
        self.dbo.sql_exec(sql)
        row = self.dbo.fetchone()
        nu_pos = row[0]
        sql = "select count(*) from mx_notused_variants where reasonId = -1";
        self.dbo.sql_exec(sql)
        row = self.dbo.fetchone()
        nu_neg = row[0]
        print("\nThere are %s pos not-used variants" % nu_pos)
        print("There are %s neg not-used variants" % nu_neg)

        #Note: remove dupes - they aren't kept in the matrix either
        self.mx_remove_dupes()

        #Note: counts for stdout
        print("Total matrix kits: %s" % len(self.KITS))
        print("Total perfect matrix variants: %s" % len(self.get_perfect_variants_idx()))
        print("Total imperfect matrix variants: %s" % len(self.get_imperfect_variants_idx()))
        print("Total matrix variants: %s\n" % len(self.VARIANTS))

        #Note: done (stdout)
        print("end MatrixData create: %s" % format(time.clock()))
        if config['DBG_MATRIX']:
            print("\n---------------------------------------")
            print("Matrix: end of create data routine")
            print("---------------------------------------")
            print(self.NP)
        
    def mx_remove_dupes(self,auto_nonsplits=False):
        #Note: this is where dupe variant profiles are removed from the matrix and stored in the DB

        def return_counts(idx, inv):
            count = np.zeros(len(idx), np.int)
            np.add.at(count, inv, 1)
            return count

        #create variables that we'll need for finding duplicates among the variants
        b = np.ascontiguousarray(self.NP).view(np.dtype((np.void, self.NP.dtype.itemsize * self.NP.shape[1])))
        _, idx, inv = np.unique(b, return_index=True, return_inverse=True)

        cnts = return_counts(idx,inv)
        dupes = [(i,j) for i,j in enumerate(idx) if cnts[i] > 1]

        if auto_nonsplits is False:
            print("\nThere are %s variants pre-dupe" % len(self.NP))
            if len(dupes) > 0:
                print("There are dupes: moving them\n")
            else:
                print("There are 0 dupes\n")

        idx_uniq = [i for i, j in enumerate(idx) if cnts[i] == 1]
        idx_dupe = [i for i, j in enumerate(idx) if cnts[i] > 1]

        dupe_cnt = 0
        for itm in dupes:

            #this gets the duplicate variants (based on idx order) 
            itms = list(np.delete(np.argwhere(inv==itm[0]),0))
            new_vid = self.get_vid_by_vix(itm[1])
            old_vids = self.get_vid_by_vix(itms)

            #add new dupe variants to mx_dupe_variants table
            sql = "insert into mx_dupe_variants(vID,dupe_vID) values (%s,?);" % new_vid
            dupe_itms = [tuple([l]) for l in old_vids]
            self.dbo.sql_exec_many(sql,dupe_itms)

            #track what's been changed
            dupe_cnt = dupe_cnt + len(dupe_itms)
            if config['DBG_DUPE_MSGS']:
                print("removing %s dupes (total: %s)" % (len(dupe_itms),dupe_cnt))

            #reset primary vID in the mx_dupe_variants table when there are dupe removals
            sql = "update mx_dupe_variants set vID = %s where vID in (%s);" % (new_vid,l2s(old_vids))
            self.dbo.sql_exec(sql)

            #reset any sup references in the mx_sups_subs table for removed variants
            sql = "update mx_sups_subs set sup = %s where sup in (%s);" % (new_vid,l2s(old_vids))
            self.dbo.sql_exec(sql)

            #reset any sub references in the mx_sups_subs table for removed variants
            sql = "update mx_sups_subs set sub = %s where sub in (%s);" % (new_vid,l2s(old_vids))
            self.dbo.sql_exec(sql)

            #remove dupe variants from the self.VARIANTS var
            for k in self.get_vname_by_vix(itms):
                self.VARIANTS.pop(k,None)

        if config['DBG_DUPE_MSGS']:
            print("")

        #reset the matrix idxs for the remaining non-dupe variants
        for ix,vn in enumerate(self.get_vname_by_vix(idx)):
            self.VARIANTS[vn] = (self.VARIANTS[vn][0],ix,self.VARIANTS[vn][2])
            
        #reset the matrix (now that the dupes are out)
        self.NP = self.NP[idx,]
        if dupe_cnt > 0:
            print("Total dupes removed: %s\n" % dupe_cnt)
        
    def reset_ss_data(self):
        #Note: this is where sups and subs are assembled, cached, and put into the DB

        #reset tbl
        sql = "delete from mx_sups_subs"
        self.dbo.sql_exec(sql)

        #perfect variants only
        NP = self.NP[self.get_perfect_variants_idx()]
        nrows, ncols = NP.shape
        sups_ = []
        vt = Variant()
        for vix in list(range(nrows)):
            vid = self.get_vid_by_vix(vix)
            vt.vix = vix
            vt.sort = self
            sups = vt.get_rel(relType=1,allowImperfect=False)
            if len(sups):
                for sup_i in self.get_vid_by_vix(sups):
                    sups_.append((sup_i,vid))

        #sqlite tbl
        sql = "insert into mx_sups_subs (sup,sub) values (?,?);"
        self.dbo.sql_exec_many(sql,sups_)

        #panda tbl
        self.SS = pd.DataFrame(sups_,columns=['sup','sub'])

        #check consistencies of perfect variants (this needs to be after self.SS created)
        print("\nChecking consistencies for perfect variants...")
        cntErr = 0
        for vix in list(range(nrows)):
            vt.vix = vix
            vt.kpc = self.get_kixs_by_val(val=1,vix=vt.vix)
            cntErr = cntErr + vt.proc_chk(allowImperfect=False,auto_perfVariants=True)
        if cntErr>0:
            print("\n%s consistency problem(s) seen. Please resolve.\n"%cntErr)
        else:
            print("Consistency check: OK.\n")
        print("Reset sup-sub relations count: %s\n" % len(sups_))


        
