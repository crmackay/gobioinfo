#include <unity.h>
#include <../align.c>



char *query = "GTGTCAGTCAC"
char *subject = "GGTCTGTCC"

// should get
//
// GTGTCAGTCAC
// | |||x||| |
// G-GTCTGTC-C
