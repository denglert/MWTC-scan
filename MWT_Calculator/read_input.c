/*
 *  read_input.c
 *  MWT_Calculator
 *
 *  Modified by Claudio Pica on 06/09/2010.
 *  Copyright 2008,2009,2010 All rights reserved.
 *
 */

#include <memory.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "read_input.h"

void read_input(input_record_t irec[], const char *filename) {
  FILE *inputfile;
  fpos_t pos;
  char buf[256];
  int npar=0, *found=NULL;
  input_record_t *crec=NULL;

  if (irec==NULL || irec->name==NULL) return; /* no input parameter specified */

    if((inputfile=fopen(filename,"r"))==NULL) {
        printf("Failed to open input file [%s].\n",filename);
        exit(1);
    }

  /* count the input parameters */
  for (crec=irec, npar=0; crec->name!=NULL; ++crec) { ++npar; }
  found=calloc(npar,sizeof(*found)); /* set found to zero */

  do {
    int count=0;
    int nowarn=0;
    /*skip white space*/
    nowarn=fscanf(inputfile," ");

    fgetpos(inputfile,&pos);

    /* skip comments */
    nowarn=fscanf(inputfile,"//%n",&count);
    if(count==2) { nowarn=fscanf(inputfile,"%*[^\n]"); goto NEXTLOOP; } if(feof(inputfile)) break; fsetpos(inputfile,&pos);

    /* read variables as described in irec */
    for (count=0; count<npar; ++count) {
      if(irec[count].descr==NULL) continue;
      if(fscanf(inputfile, irec[count].descr, irec[count].ptr)==1) {
        found[count]=1;
        goto NEXTLOOP;
      } 
      if(feof(inputfile)) goto ENDLOOP; 
      fsetpos(inputfile,&pos);
    }

    nowarn=fscanf(inputfile,"%s",buf);
    fprintf(stderr,"Ignoring unknown token: [%s]\n",buf);

NEXTLOOP:

    if(ferror(inputfile)) {
      printf("Cannot read input file.\n");
      perror(0);
    }

  } while (!feof(inputfile));

ENDLOOP:

  fclose(inputfile);

  while(npar>0) {
    --npar;
    if (found[npar]==0 && irec[npar].descr!=NULL) 
      fprintf(stderr,
          "Warning: using default value for parameter [%s] not specified in [%s]!\n",
          irec[npar].name, filename );
  }
  free(found);

}





