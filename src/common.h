#ifndef COMMON_H
#define COMMON_H










#include <stdio.h>

#define MAX_STRING_LENGTH 1000
#define TSIZE 1048576
#define SEED 1159241
#define HASHFN bitwisehash

typedef double real;
typedef struct cooccur_rec {
    int word1;
    int word2;
    real val;
} CREC;
typedef struct hashrec {
    char *word;
    long long num; 
    struct hashrec *next;
} HASHREC;


int scmp( char *s1, char *s2 );
unsigned int bitwisehash(char *word, int tsize, unsigned int seed);
HASHREC **inithashtable();
int get_word(char *word, FILE *fin);
void free_table(HASHREC **ht);
int find_arg(char *str, int argc, char **argv);
void free_fid(FILE **fid, const int num);


int log_file_loading_error(char *file_description, char *file_name);

#endif 

