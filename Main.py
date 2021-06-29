#!/usr/bin/python3
######################################
#     MEMDIS (MEMbrane DISorder)     #
#  Prediction of disordered regions  #
#     in transmembrane proteins      #
#            Laszlo Dobson           #
#       Institute of Enzymology      #
#                2021                #
######################################
import sys
import os

installDir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(installDir)
import utils
##########################
#MODIFY THE FOLLOWING PATHS
AAINDEX=installDir+"/indices"
#USAGE
#Main.py [FASTA] [TOPOLOGY] [SURFACE] [MODE]
#e.g.
#sensitive settings
#Main.py test.fas [test.xml/test.top] test.rsa test.psi test.json sens
#specific settings
#Main.py test.fas [test.xml/test.top] test.rsa test.psi test.json sens
##########################
FASTA=sys.argv[1]#single sequence fasta file 
TOPOLOGY=sys.argv[2]#CCTOP topology file (XML format, file extension must be xml) or single sequence topology file (fasta format). 
NETSURFP=sys.argv[3]#NETSURFP result file
PSIBLAST=sys.argv[4]#PSIBLAST file
SEGRES=sys.argv[5]#Platoloco res file (json format)
MODE=sys.argv[6]#sensitive/specific
##########################
#STEP 1
[Inside,inside,outside,Outside,lstm]=utils.readmod(installDir)
#STEP 2
[ID,sequence]=utils.readfasta(FASTA)
#STEP 3
#RUN CCTOP or user can define test.top 
#STEP 4
#RUN NETSURFP-1.0 "netsurfp -i test.fas -a -o test.rsa"
#STEP 5
#RUN PSIBLAST "blastpgp -j 3 -i test.fas -d uniprot_sprot.fas -e 1e-3 -Q test.psi"
#STEP 6
#RUN Platoloco  "python3 wrapper -i test.fas -m SEG >test.seg"
#STEP 7
topres=utils.readtopology(TOPOLOGY)
#STEP 8
surfres=utils.readnetsurfp(NETSURFP)
#STEP 9
[nrhits,conservation,pssm]=utils.readpsiblast(PSIBLAST)
#STEP 10
segres=utils.readseg(SEGRES)
#STEP 11
[tmregions,tailregion,segmentregion,distance,membranethickness]=utils.calctop(topres)
#STEP 12
slices=utils.cut(sequence)
#STEP 13 
[molweight,isoelectric,instability]=utils.protparam(slices)
#STEP 14 
nraaindex=utils.readAAIndex(AAINDEX)
#STEP 15
calculatedaaindex=utils.indexing(slices,nraaindex)
#STEP 16
if utils.checktop(topres)==False:
	print("not a TM protein")
else:
#STEP 17
	matrix=utils.castmx(nrhits,conservation,pssm,topres,tmregions,tailregion,segmentregion,distance,membranethickness,surfres,molweight,isoelectric,instability,calculatedaaindex,segres)	
#STEP 18
	[idr,numtop]=utils.cnn(Inside,inside,outside,Outside,sequence,matrix)
	prediction=[]
	if MODE=="spec":
#STEP 19
		prediction=utils.rescale(idr)

	
	if MODE=="sens":
#STEP 20
		idr2=utils.bidirlstm(lstm,idr,numtop)
#STEP 21
		prediction=utils.rescale2(idr2)
#STEP 22
	utils.fasprint(sequence,topres,prediction)


"""
makemx
predict
postprocess
#rmdir
"""
