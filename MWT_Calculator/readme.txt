=========================================
  NMWT Calculator for MadGraph/MadEvent

     authors: Eugenio Del Nobile
              Tuomas Hapola
              Claudio Pica


     website: http://cp3-origins.dk
=========================================


To compile the calculator type "make".

Modify the model's parameter input file (an example is given in the file "parameter_input"). 
For the parameters not specified in the input file, a default value is used.

To run the calculator and generate the param_card.dat file, type:

./MWT_Calculator parameter_input > param_card.dat
