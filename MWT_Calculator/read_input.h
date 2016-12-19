/*
 *  read_input.h
 *  MWT_Calculator
 *
 *  Created by Claudio Pica on 06/09/2010.
 *  Copyright 2010 All rights reserved.
 *
 */


typedef enum _datatype_t {
    INT_T,
    UNSIGNED_T,
    DOUBLE_T,
    STRING_T
} datatype_t;

typedef struct _input_record_t {
    char *name;
    char *descr;
    datatype_t type;
    void *ptr;
} input_record_t;

void read_input(input_record_t irec[], const char *filename);
