# Title: maf2comp
# Purpose: from maf file, to generate comparision files for genoPlotR
# Writer: Chul Lee clee03@rockefeller.edu
# Date: 20221020

# Global variable
nGene = "RBFOX1"
fNAME_tre = 'shortNAME.tre'

# Library
import sys
import os
import glob

# Functions
def make_sIDlist(fNAME_tre):
    fpin = open(fNAME_tre,'r')
    line = fpin.readline()
    part = line.strip('\n').split(",")
    sID_list = []
    for sInfo in part:
        sInfo = sInfo.replace("(","")
        sID = sInfo.split(":")[0]
        sID_list.append(sID)
    return(sID_list)

def maf2matrix(nGene):
    sID_list = make_sIDlist(fNAME_tre)
    for i in range(len(sID_list)-1):
        sID = sID_list[i]
        sID2 = sID_list[i+1]
        tmp_list = sID_list
        fNAME_Maf = "paired_maf/"+sID+"_"+sID2+".maf"
        fNAME_MafMatrix = "paired_matrix/"+sID+"_"+sID2+".maf.matrix"
        fpin = open(fNAME_Maf,'r')
        fpout = open(fNAME_MafMatrix,'w')
        # col: species / row: block
        for line in fpin:
            if len(line)>0:
                if line[0]=="a":
                    fpout.write('\t'.join(tmp_list)+'\n')
                    tmp_list = []
                    for i in range(len(sID_list)):
                        tmp_list.append('') 
                elif line[0]=="s":
                    if ".N" in line:
                        part = line.strip('\n').split('\t')
                        sID = part[1].split(".N")[0]
                        iStart = int(part[2])+1
                        iLen_block = int(part[3])-1
                        nStrand = part[4]
                        if nStrand == "+":
                            iEnd = iStart + iLen_block
                        elif nStrand == "-":
                            iStart = iStart -iLen_block
                            iEnd = iStart
                        if tmp_list[sID_list.index(sID)]=="":
                            tmp_list[sID_list.index(sID)] = str(iStart)+"_"+str(iEnd)+"_"+nStrand
                        else:
                            tmp_list[sID_list.index(sID)] += ','+str(iStart)+"_"+str(iEnd)+"_"+nStrand
        fpout.write('\t'.join(tmp_list)+'\n')
        fpin.close()
        fpout.close()

def maf2comp(nGene):
    maf2matrix(nGene)
    sID_list = make_sIDlist(fNAME_tre)
    for i in range(len(sID_list)-1):
        s1 = sID_list[i]
        s2 = sID_list[i+1]
        fNAME_comp = "comp/comparison"+str(i+1)+".comp"
        fpout = open(fNAME_comp,'w')
        headline = "start1\tend1\tstart2\tend2\tdirection\n"
        fpout.write(headline)
        
        fNAME_MafMatrix =  "paired_matrix/"+s1+"_"+s2+".maf.matrix"
        fpin = open(fNAME_MafMatrix,'r')
        fpin.readline()
        for line in fpin:
            part = line.strip('\n').split("\t")
            s1_info = part[i]
            s2_info = part[i+1]
            if not s1_info == '':
                if not s2_info =='':
                    start1 = s1_info.split("_")[0]
                    end1 = s1_info.split("_")[1]
                    strand1 = s1_info.split("_")[2]

                    start2 = s2_info.split("_")[0]
                    end2 = s2_info.split("_")[1]
                    strand2 = s2_info.split("_")[2]

                    if strand1 == strand2:
                        direction = "1"
                    else:
                        direction = "-1"
                    tmpline = start1+'\t'+end1+'\t'+start2+'\t'+end2+'\t'+direction+'\n'
                    fpout.write(tmpline)
        fpin.close()
        fpout.close()
            
def main():
    maf2comp(nGene)

if __name__=="__main__":
    main()

# End of code
