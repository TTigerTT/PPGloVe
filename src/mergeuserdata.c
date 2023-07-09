#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "common.h"

typedef struct cooccur_rec_id {
    int word1;
    int word2;
    real val;
    int id;
} CRECID;

int verbose = 2; 

char *file_head;





int write_chunk(CREC *cr, long long length, FILE *fout) {
    if (length == 0) return 0;

    long long a = 0;
    CREC old = cr[a];
    
    for (a = 1; a < length; a++) {
        if (cr[a].word1 == old.word1 && cr[a].word2 == old.word2) {
            old.val += cr[a].val;
            continue;
        }
        fwrite(&old, sizeof(CREC), 1, fout);
        old = cr[a];
    }
    fwrite(&old, sizeof(CREC), 1, fout);
    return 0;
}


int compare_crec(const void *a, const void *b) {
    int c;
    if ( (c = ((CREC *) a)->word1 - ((CREC *) b)->word1) != 0) return c;
    else return (((CREC *) a)->word2 - ((CREC *) b)->word2);
    
}


int compare_crecid(CRECID a, CRECID b) {
    int c;
    if ( (c = a.word1 - b.word1) != 0) return c;
    else return a.word2 - b.word2;
}


void swap_entry(CRECID *pq, int i, int j) {
    CRECID temp = pq[i];
    pq[i] = pq[j];
    pq[j] = temp;
}


void insert(CRECID *pq, CRECID new, int size) {
    int j = size - 1, p;
    pq[j] = new;
    while ( (p=(j-1)/2) >= 0 ) {
        if (compare_crecid(pq[p],pq[j]) > 0) {swap_entry(pq,p,j); j = p;}
        else break;
    }
}


void delete(CRECID *pq, int size) {
    int j, p = 0;
    pq[p] = pq[size - 1];
    while ( (j = 2*p+1) < size - 1 ) {
        if (j == size - 2) {
            if (compare_crecid(pq[p],pq[j]) > 0) swap_entry(pq,p,j);
            return;
        }
        else {
            if (compare_crecid(pq[j], pq[j+1]) < 0) {
                if (compare_crecid(pq[p],pq[j]) > 0) {swap_entry(pq,p,j); p = j;}
                else return;
            }
            else {
                if (compare_crecid(pq[p],pq[j+1]) > 0) {swap_entry(pq,p,j+1); p = j + 1;}
                else return;
            }
        }
    }
}


int merge_write(CRECID new, CRECID *old, FILE *fout) {
    if (new.word1 == old->word1 && new.word2 == old->word2) {
        old->val += new.val;
        return 0; 
    }
    fwrite(old, sizeof(CREC), 1, fout);
    *old = new;
    return 1; 
}


int merge_files(int num) {
    int i, size;
    long long counter = 0;
    CRECID *pq, new, old;
    char filename[200];
    FILE **fid, *fout;
    fid = calloc(num, sizeof(FILE));
    pq = malloc(sizeof(CRECID) * num);
    fout = stdout;
    if (verbose > 1) fprintf(stderr, "Merging cooccurrence files: processed 0 lines.");
    
    
    for (i = 0; i < num; i++) {
        sprintf(filename,"%s%d_20mb.bin",file_head,i);
        fid[i] = fopen(filename,"rb");
        if (fid[i] == NULL) {log_file_loading_error("file", filename); free_fid(fid, num); free(pq); return 1;}
        fread(&new, sizeof(CREC), 1, fid[i]);
        new.id = i;
        insert(pq,new,i+1);
    }
    
    
    size = num;
    old = pq[0];
    i = pq[0].id;
    delete(pq, size);
    fread(&new, sizeof(CREC), 1, fid[i]);
    if (feof(fid[i])) size--;
    else {
        new.id = i;
        insert(pq, new, size);
    }
    
    
    while (size > 0) {
        counter += merge_write(pq[0], &old, fout); 
        if ((counter%100000) == 0) if (verbose > 1) fprintf(stderr,"\033[39G%lld lines.",counter);
        i = pq[0].id;
        delete(pq, size);
        fread(&new, sizeof(CREC), 1, fid[i]);
        if (feof(fid[i])) size--;
        else {
            new.id = i;
            insert(pq, new, size);
        }
    }
    fwrite(&old, sizeof(CREC), 1, fout);
    fprintf(stderr,"\033[0GMerging cooccurrence files: processed %lld lines.\n",++counter);
    for (i=0;i<num;i++) {
        sprintf(filename,"%s_%04d_20mb.bin",file_head,i);
        
    }
    fprintf(stderr,"\n");
    free_fid(fid, num);
    free(pq);
    return 0;
}




int main(int argc, char **argv) {
    int i;
    real rlimit, n = 1e5;
    
    file_head = malloc(sizeof(char) * MAX_STRING_LENGTH);

    if (argc == 1) {
  
    
        printf("\t-verbose <int>\n");
        printf("\t\tSet verbosity: 0, 1, 2 (default), or 3\n");
        printf("\t-overflow-file <file>\n");
        printf("\t\tFilename, excluding extension, for temporary files; default overflow\n");
        printf("\t-distance-weighting <int>\n");
        printf("\t\tIf <int> = 0, do not weight cooccurrence count by distance between words; if <int> = 1 (default), weight the cooccurrence count by inverse of distance between words\n");

        printf("\nExample usage:\n");
        printf("./cooccur -verbose 2 -symmetric 0 -window-size 10 -vocab-file vocab.txt -memory 8.0 -overflow-file tempoverflow < corpus.txt > cooccurrences.bin\n\n");
 
        free(file_head);
 
        return 0;
    }

    if ((i = find_arg((char *)"-verbose", argc, argv)) > 0) verbose = atoi(argv[i + 1]);
    if ((i = find_arg((char *)"-overflow-file", argc, argv)) > 0) strcpy(file_head, argv[i + 1]);
    else strcpy(file_head, (char *)"pusercoocur");

    

    
    merge_files(10);
    free(file_head);

    return 1;
}

