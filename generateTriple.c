#define _CRT_SECURE_NO_WARNINGS

#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// windows pthread.h is buggy, but this #define fixes it
#define HAVE_STRUCT_TIMESPEC
#include <pthread.h>


#include <fcntl.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#define _FILE_OFFSET_BITS 64
#define MAX_TRIPLE 1000000
#define MAX_TRIPLE1 10000


typedef struct co_Beaver_triple{
    __uint128_t a1;
    __uint128_t a2;
    __uint128_t b1[200];
    __uint128_t b2[200];
    __uint128_t c1[200];
    __uint128_t c2[200];
}COTRIPLE;

typedef struct Beaver_triple{
    __uint128_t a1;
    __uint128_t a2;
    __uint128_t b1;
    __uint128_t b2;
    __uint128_t c1;
    __uint128_t c2;
}TRIPLE;

__uint128_t llrand() {
    __uint128_t r ;
  
    for (int i = 0; i < 9; ++i) {
        r = (r << 15) | (rand() & 0x7FFF);
    }

    return r & (__uint128_t)(0-1);
    
}


int main(){

   TRIPLE* Triples=(TRIPLE*)malloc(MAX_TRIPLE * sizeof(TRIPLE));
   COTRIPLE* Cotriples=(COTRIPLE*)malloc(MAX_TRIPLE1 * sizeof(COTRIPLE));
    FILE *fin,*fp;
   for (int i=0;i<MAX_TRIPLE;i++){
    __uint128_t a=llrand();
    __uint128_t b=llrand();
     __uint128_t c=a*b;
     Triples[i].a1 =llrand();
    Triples[i].a2=a-Triples[i].a1;

     Triples[i].b1=llrand();
     Triples[i].b2=b-Triples[i].b1;

    Triples[i].c1=llrand();
    Triples[i].c2=c-Triples[i].c1;     
    }
     char triplefile[50]="GenTriples";
    fp = fopen(triplefile, "wb");
    fseek(fp, 0, SEEK_END);
    fwrite(Triples, sizeof(TRIPLE), MAX_TRIPLE, fp);
    fclose(fp);
    printf("%sfinish\n",triplefile);

 for (int i=0;i<MAX_TRIPLE1;i++){

    __uint128_t a=llrand();
    Cotriples[i].a1 =llrand();
    Cotriples[i].a2=a-Cotriples[i].a1;
    for (int j=0;j<200;j++){
    __uint128_t b=llrand();
     __uint128_t c=a*b;
     Cotriples[i].b1[j] =llrand();
     Cotriples[i].b2[j]=b-Cotriples[i].b1[j];

     Cotriples[i].c1[j] =llrand();
     Cotriples[i].c2[j]=c-Cotriples[i].c1[j];
    }


}
    char cotriplefile[50]="GenCoTriples";
    fp = fopen(cotriplefile, "wb");
    fseek(fp, 0, SEEK_END);
    fwrite(Cotriples, sizeof(COTRIPLE), MAX_TRIPLE1, fp);
    fclose(fp);
    printf("%sfinish\n",cotriplefile);

    // for(int i=0;i<MAX_TRIPLE1;i++){
    //     if((Triples[i].a1+Triples[i].a2)*(Triples[i].b1+Triples[i].b2)!=(Triples[i].c1+Triples[i].c2))
    //         printf("false\n");
    // }

    // for(int i=0;i<MAX_TRIPLE;i++){
    //     for(int j=0;j<200;j++)
    //     if((Cotriples[i].a1+Cotriples[i].a2)*(Cotriples[i].b1[j]+Cotriples[i].b2[j])!=(Cotriples[i].c1[j]+Cotriples[i].c2[j]))
    //         printf("false\n");
    // }



}