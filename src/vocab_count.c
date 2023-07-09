

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"

typedef struct vocabulary {
    char *word;
    long long count;
} VOCAB;

int verbose = 2; 
long long min_count = 1; 
long long max_vocab = 0; 



int CompareVocabTie(const void *a, const void *b) {
    long long c;
    if ( (c = ((VOCAB *) b)->count - ((VOCAB *) a)->count) != 0) return ( c > 0 ? 1 : -1 );
    else return (scmp(((VOCAB *) a)->word,((VOCAB *) b)->word));
    
}


int CompareVocab(const void *a, const void *b) {
    long long c;
    if ( (c = ((VOCAB *) b)->count - ((VOCAB *) a)->count) != 0) return ( c > 0 ? 1 : -1 );
    else return 0;
}


void hashinsert(HASHREC **ht, char *w) {
    HASHREC     *htmp, *hprv;
    unsigned int hval = HASHFN(w, TSIZE, SEED);
    
    for (hprv = NULL, htmp = ht[hval]; htmp != NULL && scmp(htmp->word, w) != 0; hprv = htmp, htmp = htmp->next);
    if (htmp == NULL) {
        htmp = (HASHREC *) malloc( sizeof(HASHREC) );
        htmp->word = (char *) malloc( strlen(w) + 1 );
        strcpy(htmp->word, w);
        htmp->num = 1;
        htmp->next = NULL;
        if ( hprv==NULL )
            ht[hval] = htmp;
        else
            hprv->next = htmp;
    }
    else {
        
        htmp->num++;
        if (hprv != NULL) {
            
            hprv->next = htmp->next;
            htmp->next = ht[hval];
            ht[hval] = htmp;
        }
    }
    return;
}

int get_counts() {
    long long i = 0, j = 0, vocab_size = 12500;
    
    char str[MAX_STRING_LENGTH + 1];
    HASHREC **vocab_hash = inithashtable();
    HASHREC *htmp;
    VOCAB *vocab;
    FILE *fid = stdin;
    
    fprintf(stderr, "BUILDING VOCABULARY\n");
    if (verbose > 1) fprintf(stderr, "Processed %lld tokens.", i);
    
    while ( ! feof(fid)) {
        
        int nl = get_word(str, fid);
        if (nl) continue; 
        if (strcmp(str, "<unk>") == 0) {
            fprintf(stderr, "\nError, <unk> vector found in corpus.\nPlease remove <unk>s from your corpus (e.g. cat text8 | sed -e 's/<unk>/<raw_unk>/g' > text8.new)");
            free_table(vocab_hash);
            return 1;
        }
        hashinsert(vocab_hash, str);
        if (((++i)%100000) == 0) if (verbose > 1) fprintf(stderr,"\033[11G%lld tokens.", i);
    }
    if (verbose > 1) fprintf(stderr, "\033[0GProcessed %lld tokens.\n", i);
    vocab = malloc(sizeof(VOCAB) * vocab_size);
    for (i = 0; i < TSIZE; i++) { 
        htmp = vocab_hash[i];
        while (htmp != NULL) {
            vocab[j].word = htmp->word;
            vocab[j].count = htmp->num;
            j++;
            if (j>=vocab_size) {
                vocab_size += 2500;
                vocab = (VOCAB *)realloc(vocab, sizeof(VOCAB) * vocab_size);
            }
            htmp = htmp->next;
        }
    }
    if (verbose > 1) fprintf(stderr, "Counted %lld unique words.\n", j);
    if (max_vocab > 0 && max_vocab < j)
        
        
        qsort(vocab, j, sizeof(VOCAB), CompareVocab);
    else max_vocab = j;
    qsort(vocab, max_vocab, sizeof(VOCAB), CompareVocabTie); 
    
    for (i = 0; i < max_vocab; i++) {
        if (vocab[i].count < min_count) { 
            if (verbose > 0) fprintf(stderr, "Truncating vocabulary at min count %lld.\n",min_count);
            break;
        }
        printf("%s %lld\n",vocab[i].word,vocab[i].count);
    }
    
    if (i == max_vocab && max_vocab < j) if (verbose > 0) fprintf(stderr, "Truncating vocabulary at size %lld.\n", max_vocab);
    fprintf(stderr, "Using vocabulary of size %lld.\n\n", i);
    free_table(vocab_hash);
    free(vocab);
    return 0;
}

int main(int argc, char **argv) {
    if (argc == 2 &&
        (!scmp(argv[1], "-h") || !scmp(argv[1], "-help") || !scmp(argv[1], "--help"))) {
        printf("Simple tool to extract unigram counts\n");
        printf("Author: Jeffrey Pennington (jpennin@stanford.edu)\n\n");
        printf("Usage options:\n");
        printf("\t-verbose <int>\n");
        printf("\t\tSet verbosity: 0, 1, or 2 (default)\n");
        printf("\t-max-vocab <int>\n");
        printf("\t\tUpper bound on vocabulary size, i.e. keep the <int> most frequent words. The minimum frequency words are randomly sampled so as to obtain an even distribution over the alphabet.\n");
        printf("\t-min-count <int>\n");
        printf("\t\tLower limit such that words which occur fewer than <int> times are discarded.\n");
        printf("\nExample usage:\n");
        printf("./vocab_count -verbose 2 -max-vocab 100000 -min-count 10 < corpus.txt > vocab.txt\n");
        return 0;
    }

    int i;
    if ((i = find_arg((char *)"-verbose", argc, argv)) > 0) verbose = atoi(argv[i + 1]);
    if ((i = find_arg((char *)"-max-vocab", argc, argv)) > 0) max_vocab = atoll(argv[i + 1]);
    if ((i = find_arg((char *)"-min-count", argc, argv)) > 0) min_count = atoll(argv[i + 1]);
    return get_counts();
}

