# Title: gff2seg.py
# Writer: Chul Lee clee03@rockefeller.edu
# Date: 20221020

import sys
import os
import glob

nGene = "RBFOX1"

def gff2seg(nGene):
    flist = glob.glob("gff/*.adj.gff")
    for fNAME in flist:
        tmplines = ''
        fpin = open(fNAME,'r')
        iStart_list = []
        iEnd_list = []
        for line in fpin:
            part = line.strip('\n').split('\t')
            nType = part[2]

            # Set up Genic region with up/down streams
            if nType == "mRNA":
                iStart = int(part[3])
                iEnd = int(part[4])
                iStrand = part[6]
                if iStrand == "+":
                    iStrand = 1
                elif iStrand == "-":
                    iStrand = -1
                else:
                    print("Error: +-strand",line)
                    sys.exit()

                fLen = (abs(iStart-iEnd) -2 )/3.0
                tmplines += nGene +'\t'+ str(iStart-10000) +'\t'+ str(iEnd+10000) +'\t'+ str(iStrand) +'\t'+ str(fLen) +'\t'+ 'NA' +'\t'+ nGene +'\t'+ "NA" +'\t'+ "NA" +'\t'+ "NA" +'\t'+ "blocks" +'\t'+ "blocks" +'\t'+ "lightgrey" +'\t'+ "1" +'\t'+ "1" +'\t'+ "8" +'\t'+ "1" +'\t'+ 'lightgrey'+'\n'

            # Read exons
            elif nType == "exon":
                iStart = int(part[3])
                iEnd = int(part[4])
                iStart_list.append(iStart)
                iEnd_list.append(iEnd)
                
                iStrand = part[6]
                if iStrand == "+":
                    iStrand = 1
                elif iStrand == "-":
                    iStrand = -1
                else:
                    print("Error: +-strand",line)
                    sys.exit()

                fLen = (abs(iStart-iEnd) -2 )/3.0
                tmplines += nGene +'\t'+ str(iStart) +'\t'+ str(iEnd) +'\t'+ str(iStrand) +'\t'+ str(fLen) +'\t'+ 'NA' +'\t'+ nGene +'\t'+ "NA" +'\t'+ "NA" +'\t'+ "NA" +'\t'+ "CDS" +'\t'+ "exons" +'\t'+ "black" +'\t'+ "1" +'\t'+ "1" +'\t'+ "8" +'\t'+ "1" +'\t'+ 'black'+'\n'
        fpin.close()

        # Read introns
        for i in range(len(iStart_list)-1):
            if iStrand == 1:
                intron_Start = iEnd_list[i] +1
                intron_End = iStart_list[i+1] -1
            else: #iStrand == -1
                intron_End   = iStart_list[i] -1
                intron_Start = iEnd_list[i+1] +1
            fLen = (abs(intron_Start-intron_End) -2 )/3.0
            tmplines += nGene +'\t'+ str(intron_Start) +'\t'+ str(intron_End) +'\t'+ str(iStrand) +'\t'+ str(fLen) +'\t'+ 'NA' +'\t'+ nGene +'\t'+ "NA" +'\t'+ "NA" +'\t'+ "NA" +'\t'+ "CDS_intron" +'\t'+ "introns" +'\t'+ "black" +'\t'+ "1" +'\t'+ "1" +'\t'+ "8" +'\t'+ "1" +'\t'+ 'black'+'\n'

        if not tmplines == '':
            sID = os.path.basename(fNAME).split(".adj.gff")[0]
            oNAME = "seg/"+sID+".seg"
            fpout = open(oNAME,'w')
            headline = "\t".join(["name","start","end","strand","length","pid","gene","synonym","product","proteinid","feature","gene_type","col","lty","lwd","pch","cex","fill"])+'\n'
            fpout.write(headline)
            fpout.write(tmplines)
            fpout.close()

def main():
    gff2seg(nGene)


if __name__ == "__main__":
    main()

# End of code
