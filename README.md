######################################
#     MEMDIS (MEMbrane DISorder)     #
#  Prediction of disordered regions  #
#     in transmembrane proteins      #
#            Laszlo Dobson           #
#       Institute of Enzymology      #
#                2021                #
######################################

USAGE
Main.py [FASTA] [TOPOLOGY] [SURFACE] [PSI-BLAST] [LOW-COMPLEXITY] [MODE]
e.g.
Main.py test.fas [test.xml/test.top] test.rsa test.psi test.json sens

FASTA
Input file must be in fasta format

TOPOLOGY
Topology file in CCTOP XML format
see: http://cctop.enzim.ttk.mta.hu/?_=/documents/direct_interface.html

SURFACE
Surface accessibility file in netsurfp-1.0 format
see: http://www.cbs.dtu.dk/services/NetSurfP-1.0/

PSI-BLAST
Psi-blast result

MODE
Sensitive/Specific prediction [sens/spec]

OUTPUT
Number Residue Topology DisorderedScore

DEPENDENCIES
1. Install Keras and TensorFlow: https://phoenixnap.com/kb/how-to-install-keras-on-linux
