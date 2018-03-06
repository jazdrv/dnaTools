# license + libs{{{

# Purpose: Y-DNA NGS analytics
# Git repo: https://github.com/jazdrv/dnaTools
# For free distribution under the terms of the GNU General Public License,
# version 3 (29 June 2007) https://www.gnu.org/licenses/gpl.html

import sys,os,yaml,csv,json,numpy as np
from beautifultable import BeautifulTable
from collections import OrderedDict
import time,json
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
        self.sort.restore_mx_data()
        self.vix = self.proc_vname(vname)
        if self.vix is None:
            print ("No variants could be found. Exiting.")
            sys.exit()
        self.proc_chk(allowImperfect=False)
        
    def proc_nonsplits(self):
        self.sort.restore_mx_data()
        vixs = self.sort.get_imperfect_variants_idx()
        
        if len(vixs):
            #need to push to something to something other than vix here...
            vids = self.sort.get_vid_by_vix(vixs)
            cntErr = 0
            for v in vids:
                #and then convert it back, because the loop will keep chopping away vixes
                self.vix = self.sort.get_vixs_by_vids(v)
                cntErr = cntErr + self.proc_chk(allowImperfect=False,auto_nonsplits=True)
                self.upd_unk(auto_nonsplits=True)
                self.sort.mx_vandh_sort()
                self.sort.save_mx()
            if cntErr>0:
                print("\n%s consistency problem(s) seen. Please resolve\n."%cntErr)
            else:
                print("Consistency check: OK.\n")
            print("\nNum - matrix kits: %s" % len(self.sort.KITS))
            print("Num - perfect matrix variants: %s" % len(self.sort.get_perfect_variants_idx()))
            print("Num - total matrix variants: %s\n" % len(self.sort.VARIANTS))
        
    def info(self,vname):
        self.sort.restore_mx_data()
        allowImperfect = config['allowImperfectWithVInfo']
        self.vix = self.proc_vname(vname)
        self.set_info(lev=2,allowImperfect=allowImperfect)
        self.stdout_info()
        
    def matrix(self,argL,perfOnly=False):
        self.sort.restore_mx_data()
        #perfects only matrix
        if perfOnly:
            vix = self.sort.get_perfect_variants_idx()
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
        #col/row custom matrix (if applicable)
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

    def lib_name(self,argL):
        argL = [x.upper() for x in argL]
        for a in argL[:]:
            if a.find(','):
                argL = argL + a.split(",")
                argL.remove(a)
            if a.find('/'):
                argL = argL + a.split("/")
                argL.remove(a)

        sqlw = "'"+"','".join(str(x) for x in sorted(list(set(argL))))+"'"

        #build hg38
        sql = '''
            SELECT S.snpname, V.ID, V.pos, B.buildNm, AA.allele as anc,
            DA.allele as der, IX.idx, D1.vID, D2.vID, NU.reasonId
            FROM snpnames S, build B,
            alleles AA, alleles DA, variants V
            LEFT JOIN mx_idxs IX
            ON IX.axis_id = V.ID and IX.type_id=0
            LEFT JOIN mx_dupe_variants D1 -- #is a parent to dupe "children"
            ON D1.vID = V.ID
            LEFT JOIN mx_dupe_variants D2 -- #is a dupe "child" of another vix
            ON D2.dupe_vID = V.ID
            LEFT JOIN mx_notused_variants NU
            ON NU.vID = V.ID
            WHERE
            V.anc = AA.ID and V.der = DA.ID
            and S.snpname in (%s) and V.ID = S.vID
            and B.buildNm = 'hg38'
            and V.buildID = B.ID
            ORDER BY 1;
            ''' % sqlw
        self.dbo.sql_exec(sql)
        F = self.dbo.fetchall()

        if len(F) > 0:
            print("")
            table = BeautifulTable(max_width=70)
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

                table.append_row([str(row[6]).replace('None','-')]+[row[3]]+[row[0]]+[row[1]]+[row[2]]+[row[4]]+[row[5]]+[dupeP]+[nouse])
                table.row_seperator_char = ''
                table.column_seperator_char = ''
                table.column_alignments['name'] = BeautifulTable.ALIGN_LEFT
            print(table)
            print("")

        #build hg19
        sql = '''
            SELECT S.snpname,V.ID,V.pos, B.buildNm
            FROM snpnames S, build B, variants V
            WHERE
            S.snpname in (%s)
            and V.ID = S.vID
            and B.buildNm = 'hg19'
            and V.buildID = B.ID
            ORDER BY 1;
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

        sql = "SELECT axis_id from mx_idxs where type_id=1;"
        self.dbo.sql_exec(sql)
        F = self.dbo.fetchall()
        #for row in F:
            
                 
    def lib_pos(self,argL):
        argL = [x.upper() for x in argL]
        for a in argL[:]:
            if a.find(','):
                argL = argL + a.split(",")
                argL.remove(a)
            if a.find('/'):
                argL = argL + a.split("/")
                argL.remove(a)

        sqlw = ",".join(str(x) for x in sorted(list(set(argL))))

        #build hg38
        sql = '''
            SELECT DISTINCT S.snpname, V.ID, V.pos, B.buildNm, AA.allele as
            anc, DA.allele as der, IX.idx, D1.vID, D2.vID, NU.reasonId
            FROM build B, alleles AA, alleles DA, variants V
            LEFT JOIN mx_idxs IX
            ON IX.axis_id = V.ID and IX.type_id = 0
            LEFT JOIN snpnames S
            ON S.vID = v.ID
            LEFT JOIN mx_dupe_variants D1 -- #is a parent to dupe "children"
            ON D1.vID = V.ID
            LEFT JOIN mx_dupe_variants D2 -- #is a dupe "child" of another vix
            ON D2.dupe_vID = V.ID
            LEFT JOIN mx_notused_variants NU
            ON NU.vID = V.ID
            WHERE
            V.anc = AA.ID and V.der = DA.ID
            and V.pos in (%s)
            and B.buildNm = 'hg38'
            and V.buildID = B.ID
            ORDER BY 1;
            ''' % sqlw
        self.dbo.sql_exec(sql)
        F = self.dbo.fetchall()

        if len(F) > 0:
            print("")
            table = BeautifulTable(max_width=70)
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
                table.append_row([str(row[6]).replace('None','-')]+[row[3]]+[str(row[0]).replace('None','-')]+[row[1]]+[row[2]]+[row[4]]+[row[5]]+[dupeP]+[nouse])
                table.row_seperator_char = ''
                table.column_seperator_char = ''
                table.column_alignments['name'] = BeautifulTable.ALIGN_LEFT
            print(table)
            print("")

        #build hg19
        sql = '''
            SELECT S.snpname,V.ID,V.pos, B.buildNm
            FROM build B, variants V
            LEFT JOIN snpnames S
            ON S.vID = v.ID
            WHERE
            V.pos in (%s)
            and B.buildNm = 'hg19'
            and V.buildID = B.ID
            ORDER BY 1;
            ''' % sqlw
        self.dbo.sql_exec(sql)
        F = self.dbo.fetchall()

        if len(F) > 0:
            table = BeautifulTable()
            table.column_headers = ['build']+['name']+['id']+['pos']
            for row in F:
                table.append_row([row[3]]+[row[0]]+[str(row[1]).replace('None','-')]+[row[2]])
                table.row_seperator_char = ''
                table.column_seperator_char = ''
            print(table)
                 
    def lib_id(self,argL):
        #dupe kix: 803470|805350
        argL = [x.upper() for x in argL]
        for a in argL[:]:
            if a.find(','):
                argL = argL + a.split(",")
                argL.remove(a)
            if a.find('/'):
                argL = argL + a.split("/")
                argL.remove(a)

        sqlw = ",".join(str(x) for x in sorted(list(set(argL))))

        #build hg38
        sql = '''
            SELECT DISTINCT S.snpname, V.ID, V.pos, B.buildNm,
            AA.allele as anc, DA.allele as der , IX.idx, D1.vID, D2.vID, NU.reasonId
            FROM build B, alleles AA, alleles DA, variants V
            LEFT JOIN mx_idxs IX
            ON IX.axis_id = V.ID and IX.type_id = 0
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
            and V.ID in (%s)
            and B.buildNm = 'hg38'
            and V.buildID = B.ID
            ORDER BY 1;
            ''' % sqlw
        self.dbo.sql_exec(sql)
        F = self.dbo.fetchall()

        if len(F) > 0:
            print("")
            table = BeautifulTable(max_width=70)
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
                table.append_row([str(row[6]).replace('None','-')]+[row[3]]+[str(row[0]).replace('None','-')]+[row[1]]+[row[2]]+[row[4]]+[row[5]]+[dupeP]+[nouse])
                table.row_seperator_char = ''
                table.column_seperator_char = ''
                table.column_alignments['name'] = BeautifulTable.ALIGN_LEFT
            print(table)
            print("")

        #build hg19
        sql = '''
            SELECT S.snpname,V.ID,V.pos, B.buildNm
            FROM build B, variants V
            LEFT JOIN snpnames S
            ON S.vID = V.ID
            WHERE
            V.id in (%s)
            and B.buildNm = 'hg19'
            and V.buildID = B.ID
            ORDER BY 1;
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
        
    def lib_vix(self,argL):
        argL = [x.upper() for x in argL]
        for a in argL[:]:
            if a.find(','):
                argL = argL + a.split(",")
                argL.remove(a)
            if a.find('/'):
                argL = argL + a.split("/")
                argL.remove(a)

        sqlw = ",".join(str(x) for x in sorted(list(set(argL))))

        #build hg38
        sql = '''
            SELECT S.snpname,V.ID,V.pos, B.buildNm, AA.allele as anc, DA.allele
            as der, IX.idx, D1.vID, D2.vID, NU.reasonId
            FROM build B, alleles AA, alleles DA, mx_idxs IX, variants V
            LEFT JOIN snpnames S
            ON V.ID = S.vID
            LEFT JOIN mx_dupe_variants D1 -- #is a parent to dupe "children"
            ON D1.vID = V.ID
            LEFT JOIN mx_dupe_variants D2 -- #is a dupe "child" of another vix
            ON D2.dupe_vID = V.ID
            LEFT JOIN mx_notused_variants NU
            ON NU.vID = V.ID
            WHERE
            IX.axis_id = V.ID and IX.type_id=0
            and V.anc = AA.ID and V.der = DA.ID
            and IX.idx in (%s)
            and B.buildNm = 'hg38'
            and V.buildID = B.ID
            ORDER BY 1;
            ''' % sqlw
        self.dbo.sql_exec(sql)
        F = self.dbo.fetchall()

        if len(F) > 0:
            print("")
            table = BeautifulTable(max_width=70)
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
                table.append_row([row[6]]+[row[3]]+[str(row[0]).replace('None','-')]+[row[1]]+[row[2]]+[row[4]]+[row[5]]+[dupeP]+[nouse])
                table.row_seperator_char = ''
                table.column_seperator_char = ''
                table.column_alignments['name'] = BeautifulTable.ALIGN_LEFT
            print(table)
            print("")

            for row in F:
                print("vix: %s"%row[6])
                print("-----------------")
                np.argwhere(self.NP[vix,] == row[6])
            
        #build hg19
        #sql = '''
        #    SELECT S.snpname,V.ID,V.pos, B.buildNm
        #    FROM snpnames S, build B, variants V
        #    WHERE
        #    V.id in (%s)
        #    and V.ID = S.vID
        #    and B.buildNm = 'hg19'
        #    and V.buildID = B.ID
        #    ORDER BY 1;
        #    ''' % sqlw
        #self.dbo.sql_exec(sql)
        #F = self.dbo.fetchall()

        #table = BeautifulTable()
        #table.column_headers = ['build']+['name']+['id']+['pos']
        #for row in F:
        #    table.append_row([row[3]]+[row[0]]+[row[1]]+[row[2]])
        #    table.row_seperator_char = ''
        #    table.column_seperator_char = ''
        #print(table)

    def stash(self,vname):
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
        vid = self.sort.get_vid_by_vix(vix)
        vname = self.sort.get_vname_by_vix(vix)
        sql = "insert into mx_variant_stash (ID) values (%s)" % vid
        self.dbo.sql_exec(sql)
        #move this variant to a stash tbl 
        self.sort.mx_remove_vix(vix)
        print("Moved to stash: %s [%s]"%(vname,vix))
        #rerun sups+subs and consistency chk 
        self.sort.reset_ss_data()
        #sort
        self.sort.mx_vandh_sort()
        #save data
        self.sort.save_mx()
        
    def unstash(self,vname):
        #TODO
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
        #self.sort.restore_mx_data()
        if vname.isdigit() and int(vname)<=len(self.sort.VARIANTS):
            print("\n** Assuming, you're providing a vix ID...")
            vix = int(vname)
        elif vname.isdigit():
            vix = self.sort.get_vix_by_name(vname)
        else:
            vix = self.sort.get_vix_by_name(vname.upper())
        return vix
        
    def proc_kname(self,kname):
        #self.sort.restore_mx_data()
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
        self.vixn = "%s - [%s]"%(self.name,self.vix)
        self.kpc = self.sort.get_kixs_by_val(val=1,vix=self.vix)
        self.knc = self.sort.get_kixs_by_val(val=-1,vix=self.vix)
        self.kuc = self.sort.get_kixs_by_val(val=0,vix=self.vix)
        self.kpcn = "kpc: %s [%s]"%(l2s(self.sort.get_kname_by_kix(self.kpc)),l2s(self.kpc))
        self.kncn = "knc: %s [%s]"%(l2s(self.sort.get_kname_by_kix(self.knc)),l2s(self.knc))
        self.kucn = "kuc: %s [%s]"%(l2s(self.sort.get_kname_by_kix(self.kuc)),l2s(self.kuc))
        self.kpnc = sorted(self.kpc + self.knc)

        chkKits = ['637069','510668','499807','163973','124134','293507','521793','122898']
        chkKitsP = list(set(chkKits).intersection(set(self.sort.get_kname_by_kix(self.kpc))))
        chkKitsN = list(set(chkKits).intersection(set(self.sort.get_kname_by_kix(self.knc))))
        if len(chkKitsN):
            self.chkKitsN = "chkKitsN: "+l2s(chkKitsN)
            self.chkKitsP = "chkKitsP: "+l2s(chkKitsP)
        else:
            self.chkKitsN = ''
            self.chkKitsP = ''

        self.sups = self.get_rel(relType=1,allowImperfect=allowImperfect)
        self.subs = self.get_rel(relType=-1,allowImperfect=allowImperfect)
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
        if self.chkKitsN != '':
            print("%s%s" %(sp,self.chkKitsN))
            print("%s%s" %(sp,self.chkKitsP))

        print("%s%s" %(sp,self.supsn))
        print("%s%s" %(sp,self.subsn))
        print("%s%s" %(sp,self.eqvn))

        noSupOs = 0
        noSubOs = 0
        try:
            self.supOs
        except:
            noSupOs = 1
        try:
            self.subOs
        except:
            noSubOs = 1

        if noSupOs == 0:
            for sup in self.supOs:
                if sup.name != 'top':
                    sup.stdout_info(tp=1)
        if noSubOs == 0:
            for sub in self.subOs:
                sub.stdout_info(tp=-1)

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
        if auto_perfVariants is True:
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
                print("This variant has no supsets, can't do anything with it!")
        else:

            if config['SHOW_PROC_CHK_DETAILS']:
                print("unks: %s [%s]" % (l2s(self.sort.get_kname_by_kix(self.kuc)),l2s(self.kuc)))
                print("")
                print("There are %s sups to go through.\n"%(len(self.sups)))
            for sup in reversed(self.sups):

                splitL = []
                if config['SHOW_PROC_CHK_DETAILS']:
                    print("sup: %s [%s]" % (self.sort.get_vname_by_vix(sup),sup))

                #This part is only for variants with unknowns
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

                vid = self.sort.get_vid_by_vix(sup)
                try:
                    v1subs_ = self.sort.get_vixs_by_vids(list(self.sort.SS.loc[self.sort.SS['sup']==vid,['sub']].as_matrix().T[0]))
                except:
                    v1subs_ = []

                #kpc of the sup
                v1kpc = self.sort.get_kixs_by_val(val=1,vix=sup)

                #remove any target subs or target equivalent variants from consideration (use cache)
                v1subs = list(set(v1subs_)-set(self.subs)-set(self.eqvs))

                #do consistency checks on remaining variants
                if config['SHOW_PROC_CHK_DETAILS']:
                    print("\nThere are %s v1subs for sup %s"%(len(v1subs),self.sort.get_vname_by_vix(sup)))
                if len(v1subs):

                    if config['SHOW_PROC_CHK_DETAILS']:
                        print("")
                    VAR1a = np.argwhere(self.sort.NP[:,self.kpc] == 1)
                    mask = np.isin(VAR1a[:,0], v1subs)
                    idx = list(list(np.where(mask))[0])
                    VAR2a = np.unique(VAR1a[:,0][idx]) #these are the overlapping v1subs with target v
                    for v in VAR2a:
                        akpc = self.sort.get_kixs_by_val(val=1,vix=v)
                        #(mtP) common kpc btw v1sub and target variant
                        at_common_kpc = list(set(akpc).intersection(set(self.kpc)))
                        #(msP) common kpc btw v1sub and common sup
                        as_common_kpc = list(set(akpc).intersection(set(v1kpc)))
                        #(xP) diff_common_kpc = [msP-mtP]
                        diff_common_kpc = list(set(as_common_kpc)-set(at_common_kpc))
                        #(cP) common_kp c = intersection btw msP and mtP
                        common_kpc = list(set(as_common_kpc).intersection(set(at_common_kpc)))
                        if config['SHOW_PROC_CHK_DETAILS']:
                            print(" - [2] v1sub: %s [%s]"%(self.sort.get_vname_by_vix(v),v))
                            print("       (mtP) shared btw vix + v1sub: %s [%s]"%(l2s(self.sort.get_kname_by_kix(at_common_kpc)),l2s(at_common_kpc)))
                            print("       (msP) shared btw sup + v1sub: %s [%s]"%(l2s(self.sort.get_kname_by_kix(as_common_kpc)),l2s(as_common_kpc)))
                            print("       (xP) msP-mtP: %s [%s]"%(l2s(self.sort.get_kname_by_kix(diff_common_kpc)),l2s(diff_common_kpc)))
                            print("       (cP) common btw msP+mtP: %s [%s]"%(l2s(self.sort.get_kname_by_kix(common_kpc)),l2s(common_kpc)))
                        for k in diff_common_kpc:
                            if k in unkL:
                                unkD[k][1] = unkD[k][1] + 1
                                if unkD[k][1] > 1:
                                    if config['SHOW_PROC_CHK_DETAILS']:
                                        print(" - [2] %s [%s] is a req positive" % (self.sort.get_kname_by_kix(k),k))
                        if len(common_kpc) > 0:
                            splitL.append(common_kpc)
                        
                #splits
                uniq_splits = [list(x) for x in set(tuple(x) for x in splitL)]
                if len(uniq_splits) > 0:
                    split_intersects = list(set.intersection(*map(set, uniq_splits)))
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
        if len(spl) > 0:
            if auto_perfVariants is True:
                print("[SPLIT ISSUE] vix: %s [%s] - splits: %s [%s]"%(self.sort.get_vname_by_vix(self.vix),self.vix,l2s(self.sort.get_vname_by_vix(spl)),l2s(spl)))
                print(uniq_splits)
                return 1 #split errors for perfect variant chk
            if auto_perfVariants is False and auto_nonsplits is False:
                print("splits: %s [%s]" % (l2s(self.sort.get_vname_by_vix(spl)),l2s(spl)))
                rec.update({"spl":spl})

        if auto_perfVariants is not True:
            if len(pos):
                print("pos: %s [%s]" % (l2s(self.sort.get_kname_by_kix(pos)),l2s(pos)))
                rec.update({"p":pos})
            if len(pos_a):
                print("pos_a: %s [%s]" % (l2s(self.sort.get_kname_by_kix(pos_a)),l2s(pos_a)))
                rec.update({"pa":pos_a})
            if len(neg):
                print("neg: %s [%s]" %  (l2s(self.sort.get_kname_by_kix(neg)),l2s(neg)))
                rec.update({"n":neg})
            if len(self.sups):
                print("sups: %s [%s]" % (l2s(self.sort.get_vname_by_vix(self.sups)),l2s(self.sups)))
                rec.update({"sups":self.sups})
            if len(self.subs):
                print("subs: %s [%s]" % (l2s(self.sort.get_vname_by_vix(self.subs)),l2s(self.subs)))
                rec.update({"subs":self.subs})
            sql = "delete from mx_sort_recommendations where vID=%s"%self.vix
            self.dbo.sql_exec(sql)
            if len(rec.items()):
                sql = '''insert into mx_sort_recommendations (vID,instructions) values (%s,"%s");'''%(self.vix,str(rec).replace(" ",""))
                self.dbo.sql_exec(sql)

        if len(spl) > 0 and auto_nonsplits is False:
            return 1 #split errors for imperfect variant chk
        else:
            return 0 #no split errors for imperfect variant chk

        if auto_mode is not True:
            print("")
            print("---------------------------------------------------------------------")
            print("")
        
    def upd_unk(self,argL=None,auto_nonsplits=True):
        if auto_nonsplits is False:
            self.sort.restore_mx_data()
        #arg input handling
        if auto_nonsplits is False:
            if len(argL) == 1:
                vname = argL[0]
                if vname.isdigit() and int(vname)<=len(self.sort.VARIANTS):
                    print("\n** Assuming, you're providing a vix ID...")
                    vix = int(vname)
                elif vname.isdigit():
                    vix = self.sort.get_vix_by_name(vname)
                else:
                    vix = self.sort.get_vix_by_name(vname.upper())
                self.vix = vix
            else:
                print("Missing vix criteria. Exiting.")
                sys.exit()
        #retrieve insructions for this variant (or exit, if nothing)
        sql = "select instructions from mx_sort_recommendations where vID=%s"%self.vix
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
        #positive calls
        if 'p' in recD.keys():
            for kix in recD['p']:
                if auto_nonsplits is False:
                    print("Changing: coord (%s:%s) to POS"%(self.sort.get_vname_by_vix(self.vix),self.sort.get_kname_by_kix(kix)))
                kid = self.sort.get_kid_by_kix(kix)
                sql = "select * from mx_calls where vID=%s and pID=%s"%(self.vix,kid)
                self.dbo.sql_exec(sql)
                row = self.dbo.fetchone()
                if row is not None:
                    oldval = row[4]
                    sql = "update mx_calls set assigned=1,confidence=2,changed=%s where pID=%s and vID=%s;"%(oldval,kid,self.vix)
                    self.dbo.sql_exec(sql)
                else:
                    oldval = self.sort.NP[self.vix,kix]
                    sql = "insert into mx_calls (pID,vID,assigned,confidence,changed) values (%s,%s,%s,%s,%s);"%(kid,self.vix,1,2,oldval)
                    self.dbo.sql_exec(sql)
                self.sort.NP[self.vix,kix] = 1
        #positive ambiguous calls
        if 'pa' in recD.keys():
            for kix in recD['pa']:
                if auto_nonsplits is False:
                    print("Changing: coord (%s:%s) to POS-A"%(self.sort.get_vname_by_vix(self.vix),self.sort.get_kname_by_kix(kix)))
                self.sort.NP[self.vix,kix] = 1
                kid = self.sort.get_kid_by_kix(kix)
                sql = "select * from mx_calls where vID=%s and pID=%s"%(self.vix,kid)
                self.dbo.sql_exec(sql)
                row = self.dbo.fetchone()
                if row is not None:
                    oldval = row[4]
                    sql = "update mx_calls set assigned=1,confidence=1,changed=%s where pID=%s and vID=%s;"%(oldval,kid,self.vix)
                    self.dbo.sql_exec(sql)
                else:
                    oldval = self.sort.NP[self.vix,kix]
                    sql = "insert into mx_calls (pID,vID,assigned,confidence,changed) values (%s,%s,%s,%s,%s);"%(kid,self.vix,1,1,oldval)
                    self.dbo.sql_exec(sql)
                self.sort.NP[self.vix,kix] = 1
        #negative calls
        if 'n' in recD.keys():
            for kix in recD['n']:
                if auto_nonsplits is False:
                    print("Changing: coord (%s:%s) to NEG"%(self.sort.get_vname_by_vix(self.vix),self.sort.get_kname_by_kix(kix)))
                kid = self.sort.get_kid_by_kix(kix)
                sql = "select * from mx_calls where vID=%s and pID=%s"%(self.vix,kid)
                self.dbo.sql_exec(sql)
                row = self.dbo.fetchone()
                if row is not None:
                    oldval = row[4]
                    sql = "update mx_calls set assigned=-1,confidence=2,changed=%s where pID=%s and vID=%s;"%(oldval,kid,self.vix)
                    self.dbo.sql_exec(sql)
                else:
                    oldval = self.sort.NP[self.vix,kix]
                    sql = "insert into mx_calls (pID,vID,assigned,confidence,changed) values (%s,%s,%s,%s,%s);"%(kid,self.vix,-1,2,oldval)
                    self.dbo.sql_exec(sql)
                self.sort.NP[self.vix,kix] = -1
        #vars needed for not used chk
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
        #updated variant reveals it's not usable
        if len(kpc) == 0 or len(knc) == 0:
            if auto_nonsplits is False:
                if len(kpc) == 0 : print("Updated unks for this variant reveals no positives. Achiving.")
                if len(knc) == 0 : print("Updated unks for this variant reveals no negatives. Achiving.")
            #remove variant from self.VARIANTS + matrix
            self.sort.mx_remove_vix(self.vix)
            #not used variants
            if len(knc) == 0:
                sql1 = "update mx_calls set removal=1 where vid=%s;"%self.vix
                sql2 = "insert into mx_notused_variants(vId,reasonId) values(%s,1);"%self.vix
            elif len(kpc) == 0:
                sql1 = "update mx_calls set removal=2 where vid=%s;"%self.vix
                sql2 = "insert into mx_notused_variants(vId,reasonId) values(%s,2);"%self.vix
            self.dbo.sql_exec(sql1)
            self.dbo.sql_exec(sql2)
            #reset mx_sort_recommendations
            sql = "delete from mx_sort_recommendations;"
            self.dbo.sql_exec(sql)
        else:
            #if not a dupe and not a notused variant, then ...
            vid = self.sort.get_vid_by_vix(self.vix)
            #drop any prior sup/sub definitions for this variant
            sql = "delete from mx_sups_subs where sup = %s or sub = %s;"%(vid,vid)
            self.dbo.sql_exec(sql)
            cnt_new_sups_subs = 0
            #new sups
            if 'sups' in recD.keys():
                sups_ = []
                for sup_i in self.sort.get_vid_by_vix(recD['sups']):
                    sups_.append((sup_i,vid))
                #sqlite tbl
                sql = "insert into mx_sups_subs (sup,sub) values (?,?);"
                self.dbo.sql_exec_many(sql,sups_)
                cnt_new_sups_subs = cnt_new_sups_subs + len(recD['sups'])
            #new subs
            if 'subs' in recD.keys():
                subs_ = []
                for sub_i in self.sort.get_vid_by_vix(recD['subs']):
                    subs_.append((vid,sub_i))
                #sqlite tbl
                sql = "insert into mx_sups_subs (sup,sub) values (?,?);"
                self.dbo.sql_exec_many(sql,subs_)
                cnt_new_sups_subs = cnt_new_sups_subs + len(recD['subs'])
            #recreate panda sups/subs cache
            sql = "select distinct * from mx_sups_subs;"
            self.dbo.sql_exec(sql)
            sups_ = self.dbo.fetchall()
            self.SS = pd.DataFrame(sups_,columns=['sup','sub'])
            #stdout about what's new
            if cnt_new_sups_subs > 0 and auto_nonsplits is False:
                print("Found %s new sup-sub relations" % cnt_new_sups_subs)

        #if dupe:
        #    sql = "insert into mx_dupe_variants(vID,dupe_vID) values (%s,?);" % self.get_vid_by_vix(itm[1])
        #print("kuc:%s"%kuc)
        #if len(kuc) == 0:
        #    sql = "delete from mx_sort_recommendations where vID=%s;"%self.vix)
        #    self.dbo.sql_exec(sql)
        #    self.proc_chk(allowImperfect=False)

        if auto_nonsplits is False:
            self.sort.mx_vandh_sort()
            self.sort.save_mx()

        if auto_nonsplits is False:
            print("\nNum - matrix kits: %s" % len(self.sort.KITS))
            print("Num - perfect matrix variants: %s" % len(self.sort.get_perfect_variants_idx()))
            print("Num - total matrix variants: %s\n" % len(self.sort.VARIANTS))

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

        #imperfect togglge
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
        if config['DBG_SUBS_IN']: print("[subin.12] VAR6: %s"%VAR6) #slow

        #subsets shouldn't have diff intersecting positives with common supersets
        subsc = list(VAR5)
        invalid_subs = []
        self.sups = self.get_rel(relType=1,allowImperfect=False)
        self.kpc = self.sort.get_kixs_by_val(val=1,vix=self.vix)
        for v in subsc:
            #v1 = Variant()
            #v1.dbo = self.dbo
            #v1.sort = self.sort
            #v1.vix = v
            vid = self.sort.get_vid_by_vix(v)
            #print("vid:%s"%vid)
            try:
                v1sups = self.sort.get_vixs_by_vids(list(self.sort.SS.loc[self.sort.SS['sub']==vid,['sup']].as_matrix().T[0]))
            except:
                v1sups = []
            #print("(from cache) sups of sub (%s): %s"%(v,v1sups))
            #print("12.10") #slow
            #v1sups = v1.get_rel(relType=1)
            #print("12.11")
            valid_sub = 1
            if len(v1sups):
                common_sups = set(v1sups).intersection(set(self.sups))
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
        
    def get_vid_by_vix(self,vix): #ok
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
        
    def get_kid_by_kix(self,kix): #ok
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
        
    def get_vname_by_vix(self,vix): #ok
        intFlg = True
        try:
            value = int(vix)
        except:
            intFlg = False
        if intFlg:
            vix = [vix]
        vnList = []
        for vo in vix:
            for itm in list(self.VARIANTS.items()):
                if itm[1][1] == vo:
                    vnList.append(itm[0])
                    break
        if intFlg:
            return vnList[0]
        else:
            return vnList
        
    def get_kname_by_kix(self,kix): #ok
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
            if vn in self.VARIANTS.keys():
                vixs.append(self.VARIANTS[vn][1])
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

    def get_kixs_by_val(self,val,vix=None,vname=None,overrideData=None):
        if vname is not None:
            vix = self.get_vix_by_name(vname)
        if vix is not None and overrideData is not None: # we're sending in a custom evaluation
            return list(np.argwhere(overrideData[0,] == val).T[1,]) #with override data, there's only one line evaluated - 1d datset
        if vix is not None: #no override -- use self.NP (all data)
            return list(np.argwhere(self.NP[vix,] == val).T[1,]) #default data, it's the entire matrix - 2d dataset 
        
    def get_vixs_by_val(self,val,kix=None,kname=None,overrideData=None):
        if kname is not None:
            kix = self.get_kix_by_name(kname)
        if kix is not None and overrideData is not None:
            return np.argwhere(overrideData[:,0] == val).T[0,] #with override data, there's only one line evaluated - 1d dataset
        if kix is not None: #no override -- use self.NP (all data)
            return np.argwhere(self.NP[:,kix] == val).T[0,] #default data, it's the entire matrix - 2d dataset
        
    def get_knames_by_val(self,val,vix=None,vname=None,overrideData=None,toStr=True):
        if toStr:
            return l2s(self.get_kname_by_kix(self.get_kix_by_val(val,vix,vname,overrideData)))
        else:
            return self.get_kname_by_kix(self.get_kix_by_val(val,vix,vname,overrideData))
        
    def get_vnames_by_val(self,val,kix=None,kname=None,overrideData=None,toStr=True):
        if toStr:
            return l2s(self.get_vname_by_vix(self.get_vixs_by_val(val,kix,kname,overrideData)))
        else:
            return self.get_vname_by_vix(self.get_vixs_by_val(val,kix,kname,overrideData))

    def mx_vandh_sort(self):
        if config['DBG_MATRIX']:
            print("\n=================================")
            print("inside vandh sort")
            print("=================================\n")
        def return_counts(idx,inv):
            count = np.zeros(len(idx), np.int)
            np.add.at(count, inv, 1)
            return count
        #debugging line outs
        if config['DBG_MATRIX']:
            print("\n---------------------------------------")
            print("Matrix: beg vandh sort")
            print("---------------------------------------")
            print(self.NP)
            print("")
        #perfect variants
        prfVix = self.get_perfect_variants_idx()
        #Note: needs a certain amount of perfect vix (or it breaks)
        if len(prfVix) < 2:
            print("Only %s perfect vix found. Not enough. Exiting.\n" % len(prfVix))
            sys.exit()

        #panda special handling to coord perfect vix + kix sorting
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

        #imperfect variants
        impVix = self.get_imperfect_variants_idx()
        impVixPosV = np.argwhere(self.NP[impVix]==1)[:,0]
        unqIV, cntIV = np.unique(impVixPosV, return_counts=True)
        allIV = np.asarray((impVix,cntIV)) #idx w/counts
        allIV = allIV[:,np.argsort(-allIV[1])] #vsort
        impVixSorted = allIV.T[:,0] #sorted vix

        #vsort
        VixSorted = np.concatenate((prfVixSorted,impVixSorted))
        self.NP = self.NP[np.concatenate((prfVixSorted,impVixSorted)),]

        #if there are any kits w/o positives, figure it out here
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

    def get_mx_row_as_list(self,rownum,noneToStr=True,NP=None):
        #TODO: this NP override is a hack, need to work out why I need this
        if NP is not None:
            NP = self.NP
        else:
            NP = NP
        if noneToStr:
            return ['None' if v is None else v for v in NP[rownum,:].tolist()[0]]
        else:
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
        #remove variant from self.VARIANTS
        idx = list(range(len(self.VARIANTS)))
        idx.pop(vix)
        vname = self.get_vname_by_vix(vix)
        self.VARIANTS.pop(vname,None)
        #reset the matrix idxs
        for ix,vn in enumerate(self.get_vname_by_vix(idx)):
            self.VARIANTS[vn] = (self.VARIANTS[vn][0],ix,self.VARIANTS[vn][2])
        #reset the matrix
        self.NP = self.NP[idx,]
        self.mx_vandh_sort()
        self.save_mx()
        #reset mx_sort_recommendations
        sql = "delete from mx_sort_recommendations"
        self.dbo.sql_exec(sql)

    # data 

    def sort_schema(self):
        self.dbo.db = self.dbo.db_init()
        self.dbo.dc = self.dbo.cursor()
        self.dbo.sql_exec_file('sort-schema.sql')
        
    def save_mx(self):
        #push kit/variant/numpy data into saved/matrix tbls + h5py file
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
            SELECT V.name, V.ID, IX.idx, V.pos
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
            SELECT K.kitId, K.ID, IX.idx
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

        def return_counts(idx, inv):
            count = np.zeros(len(idx), np.int)
            np.add.at(count, inv, 1)
            return count

        print("beg MatrixData create: %s" % format(time.clock()))

        #bedranges - uses in_range routine
        sql = "drop index if exists tmp1idx1;"
        self.dbo.sql_exec(sql)
        sql = "drop index if exists tmp1idx2;"
        self.dbo.sql_exec(sql)
        sql = "drop table if exists tmp1;"
        self.dbo.sql_exec(sql)
        sql = "create table tmp1 (pid int,vid int,pos int,val bool);"
        self.dbo.sql_exec(sql)
        sql = "CREATE INDEX tmp1idx1 ON tmp1(pid,vid);"
        self.dbo.sql_exec(sql)
        sql = "CREATE INDEX tmp1idx2 ON tmp1(pid,vid,val);"
        self.dbo.sql_exec(sql)
        sql = "SELECT DISTINCT C.pid FROM vcfcalls C"
        pids = self.dbo.sql_exec(sql)
        F = self.dbo.fetchall()
        for row in F:
            sql = "SELECT * from v_pos_call_chk;"
            calls = self.dbo.dc.execute(sql)
            self.dbo.rc = self.dbo.cursor()
            sql = '''select minaddr,maxaddr from bedranges r
                inner join bed b on b.bID=r.id
                where b.pid=%s order by 1'''%row[0]
            ranges = self.dbo.rc.execute(sql)
            xl = list([(v[0],v[1]) for v in calls])
            pv = [v[1] for v in xl]
            rv = [v for v in ranges]
            coverage = in_range(pv, rv, spans=None)
            mergedList = [(x,y,z) for (x,y),z in zip(xl, coverage)]
            sql = "insert into tmp1 (pid,vid,pos,val) values (%s,?,?,?);"%row[0]
            self.dbo.sql_exec_many(sql,mergedList)

        sql = "drop table if exists tmp2;"
        self.dbo.sql_exec(sql)
        #this is the data that correlates variants pos's to whether kits were
        #tested for them or not (a neg is when it was covered, but it's not
        #showing up in the vcfcalls)
        sql = '''
            CREATE TABLE tmp2 as 
            SELECT DISTINCT vid as vID, pid as pID,
            cast(pid as varchar(15))||'|'||cast(vid as varchar(15)) as pidvid
            FROM tmp1 WHERE val=1;'''
        self.dbo.sql_exec(sql)
        #we don't need to chk for variants we already know have negs (that we
        #can see in the vcfs)
        sql = "delete from tmp2 where vid in (select distinct vID from v_neg_call_chk1);"
        self.dbo.sql_exec(sql)
        #Note: this needs to be done because some pos variant/kit combos in the vcf's are missing in the beds
        sql = "delete from tmp2 where pidvid in (select pidvid from v_pos_call_chk_with_kits);"
        self.dbo.sql_exec(sql)
        #Note: I believe this also needs to be done, since there are some ambiguous calls
        sql = "delete from tmp2 where pidvid in (select pidvid from v_unk_call_chk_with_kits);"
        self.dbo.sql_exec(sql)

        #sql to fetch matrix variant data
        sql = "select * from v_imx_assignments_with_unk;"
        self.dbo.sql_exec(sql)
        F = self.dbo.fetchall()

        #retrieve data from sqlite like so: [V][K] [x,x,x,x,x,x,...]
        DATA = OrderedDict()
        self.KITS = {}
        self.VARIANTS = {}
        cntV = 0
        cntK = 0

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

        self.NP = np.matrix(list(DATA.values()))
        if config['DBG_MATRIX']:
            print("\n---------------------------------------")
            print("Matrix: just created data")
            print("---------------------------------------")
            print(self.NP)

        #mx_notused_variants - these are variants that don't have at least one
        #neg and one positive
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

        #create variables that we'll need for finding duplicates among the variants
        b = np.ascontiguousarray(self.NP).view(np.dtype((np.void, self.NP.dtype.itemsize * self.NP.shape[1])))
        _, idx, inv = np.unique(b, return_index=True, return_inverse=True)

        cnts = return_counts(idx,inv)
        dupes = [(i,j) for i,j in enumerate(idx) if cnts[i] > 1]

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
            #this puts these duplicate variants into the dupes tbl next to the one in their series we're keeping
            sql = "insert into mx_dupe_variants(vID,dupe_vID) values (%s,?);" % self.get_vid_by_vix(itm[1])
            dupe_itms = [tuple([l]) for l in self.get_vid_by_vix(itms)]
            dupe_cnt = dupe_cnt+len(dupe_itms)
            if config['DBG_DUPE_MSGS']:
                print("removing %s dupes (total: %s)" % (len(dupe_itms),dupe_cnt))
            self.dbo.sql_exec_many(sql,dupe_itms)
            #remove dupe variants from the self.VARIANTS var
            for k in self.get_vname_by_vix(itms):
                self.VARIANTS.pop(k,None)

        if config['DBG_DUPE_MSGS']:
            print("")

        #reset the matrix idxs for the remaining non-dupe variants
        for ix,vn in enumerate(self.get_vname_by_vix(idx)):
            self.VARIANTS[vn][1] = ix
            
        #reset the matrix (now that the dupes are out)
        self.NP = self.NP[idx,]
        print("Total dupes removed: %s\n" % dupe_cnt)
        print("Num - matrix kits: %s" % len(self.KITS))
        print("Num - perfect matrix variants: %s" % len(self.get_perfect_variants_idx()))
        print("Num - total matrix variants: %s\n" % len(self.VARIANTS))

        print("end MatrixData create: %s" % format(time.clock()))
        #self.stdout_matrix()
        if config['DBG_MATRIX']:
            print("\n---------------------------------------")
            print("Matrix: end of create data routine")
            print("---------------------------------------")
            print(self.NP)

        #save matrix data
        #self.save_mx()
        
    def reset_ss_data(self):
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
            print("\n%s consistency problem(s) seen. Please resolve\n."%cntErr)
        else:
            print("Consistency check: OK.\n")
        print("Found %s new sup-sub relations\n" % len(sups_))


        
