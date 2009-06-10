#!/bin/bash
icc *.c -Nmpi -lVT -I$VT_ROOT/include -L$VT_LIB_DIR $VT_ADD_LIBS
