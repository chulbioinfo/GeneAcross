## Title: maf2omicCircosInput
# Purpose: from maf file, to generate input files omicCircos (R) including seg.f and link polygon (link.pg)
# Writer: Chul Lee clee03@rockefeller.edu
# Date: 20221025

# Global variable

fNAME_tre = 'shortNAME.tre'
UpDownStreamSize = 10000
iLen_block_min = 0
# Library
import sys
import os
import glob
from itertools import combinations

# Functions
# Function which returns subset or r length from n
def tre2SpeciesList(fNAME_tre):
    fpin = open(fNAME_tre,'r')
    line = fpin.readline()
    part = line.strip('\n').split(",")
    sID_list = []
    for sInfo in part:
        sInfo = sInfo.replace("(","")
        sID = sInfo.split(":")[0]
        sID_list.append(sID)
    return(sID_list)

def tre2seg(fNAME_tre,UpDownStreamSize):
    fpin = open(fNAME_tre,'r')
    line = fpin.readline()
    part = line.strip('\n').split(",")
    sID_list = []
    for sInfo in part:
        sInfo = sInfo.replace("(","")
        sID = sInfo.split(":")[0]
        sID_list.append(sID)
    fpout = open("omicCircos/input/seg.f.tsv",'w')
    headline = "seg.name\tseg.Start\tseg.End\tthe.v\tNO\n"
    fpout.write(headline)
    #s1	0	10000	NA	NA
    #s2	0	10000	NA	NA
    #s3	0	10000	NA	NA
    #s4	0	10000	NA	NA
    for sID in sID_list:
        tmpline = sID+'\t'+"0\t"+str(UpDownStreamSize)+"\tNA\tNA\n"
        fpout.write(tmpline)
    fpout.close()
    

def maf2linkmatrix(fNAME_tre,iLen_block_min):
    sID_list = tre2SpeciesList(fNAME_tre)
    for i in range(len(sID_list)):
        sID = sID_list[i]
        tmp_list = sID_list
        fNAME_Maf = "maf/RBFOX1_"+sID+".maf"
        fNAME_LinkMatrix = "omicCircos/work/"+sID+"link.matrix.tsv"
        fpin = open(fNAME_Maf,'r')
        fpout = open(fNAME_LinkMatrix,'w')
        # col: species / row: block
        for line in fpin:
            if len(line)>0:
                if line[0]=="a":
                    fpout.write('\t'.join(tmp_list)+'\n')
                    tmp_list = []
                    for i in range(len(sID_list)):
                        tmp_list.append('')
                elif line[0]=="s":
                    if not ".Anc" in line:
                        if ".N" in line:
                            part = line.strip('\n').split('\t')
                            sID = part[1].split(".N")[0]
                            iStart = int(part[2])+1
                            iLen_block = int(part[3])-1
                            nStrand = part[4]
                            nSeq = part[-1]
                            iCNT_gap = nSeq.count("-")
                            if iLen_block >= iLen_block_min:
                                if iCNT_gap >= iLen_block/2.0:
                                    if nStrand == "+":
                                        iEnd = iStart + iLen_block
                                    elif nStrand == "-":
                                        iStart = iStart -iLen_block
                                        iEnd = iStart
                                    if tmp_list[sID_list.index(sID)]=="":
                                        tmp_list[sID_list.index(sID)] = str(iStart)+"_"+str(iEnd)+"_"+nStrand
                                    else:
                                        tmp_list[sID_list.index(sID)] += ','+str(iStart)+"_"+str(iEnd)+"_"+nStrand
                                        #print("Duplication_in_a_blcok:"+sID+"_"+str(iStart)+"_"+str(iEnd)+"_"+nStrand)
        fpout.write('\t'.join(tmp_list)+'\n')
        fpin.close()
        fpout.close()

def linkmatrix2linkpg(fNAME_tre,UpDownStreamSize):
    sID_list = tre2SpeciesList(fNAME_tre)
    fNAME_linkpg = "omicCircos/input/link.pg.tsv"
    fpout = open(fNAME_linkpg,'w')
    headline = "seg1\tstart1\tend1\tseg2\tstart2\tend2\n"
    fpout.write(headline)

    for i in range(len(sID_list)):
        sID = sID_list[i]
        fNAME_LinkMatrix = "omicCircos/work/"+sID+"link.matrix.tsv"
        fpin = open(fNAME_LinkMatrix,'r')
        fpin.readline()
        for line in fpin:
            part = line.strip('\n').split("\t")
            BlockElement_list = []
            for i in range(len(sID_list)):
                sInfo = part[i]
                if not sInfo == '':
                    if not ',' in sInfo:
                        start1 = int(sInfo.split("_")[0])
                        end1 = int(sInfo.split("_")[1])
                        strand1 = sInfo.split("_")[2]
                        if start1 <= UpDownStreamSize:
                            if end1 <= UpDownStreamSize:
                                BlockElement_list.append(str(i)+"_"+sInfo)
                            else:
                                adj_sInfo = str(start1)+"_"+str(UpDownStreamSize)+"_"+strand1
                                BlockElement_list.append(str(i)+"_"+adj_sInfo)
                    else:
                        sInfos = sInfo.split(",")
                        for sInfo in sInfos:
                            start1 = int(sInfo.split("_")[0])
                            end1 = int(sInfo.split("_")[1])
                            strand1 = sInfo.split("_")[2]
                            if start1 <= UpDownStreamSize:
                                if end1 <= UpDownStreamSize:
                                    BlockElement_list.append(str(i)+"_"+sInfo)
                                else:
                                    adj_sInfo = str(start1)+"_"+str(UpDownStreamSize)+"_"+strand1
                                    BlockElement_list.append(str(i)+"_"+adj_sInfo)                         

            if len(BlockElement_list)>=2:
                for link_pair in list(combinations(BlockElement_list, 2)):
                    sInfo1_list = link_pair[0].split("_")
                    s1_seg   = sID_list[int(sInfo1_list[0])]
                    s1_start = sInfo1_list[1]
                    s1_end   = sInfo1_list[2]
                    
                    sInfo2_list = link_pair[1].split("_")
                    s2_seg   = sID_list[int(sInfo2_list[0])]
                    s2_start = sInfo2_list[1]
                    s2_end   = sInfo2_list[2]

                    tmpline = s1_seg +"\t"+ s1_start +"\t"+ s1_end +"\t"+ s2_seg +"\t"+ s2_start +"\t"+ s2_end +"\n"
                    fpout.write(tmpline)
        fpin.close()
    fpout.close()
            
def main():
    tre2seg(fNAME_tre,UpDownStreamSize)
    maf2linkmatrix(fNAME_tre,iLen_block_min)
    linkmatrix2linkpg(fNAME_tre,UpDownStreamSize)

if __name__=="__main__":
    main()

# End of code
