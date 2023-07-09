

#define _CRT_SECURE_NO_WARNINGS

#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


#define HAVE_STRUCT_TIMESPEC
#include <pthread.h>

#include "common.h"

#define _FILE_OFFSET_BITS 64
#define MAX_TRIPLE1 10000

typedef struct securecooccur_rec {
    int word1;
    int word2;
    __uint128_t pt1;
    __uint128_t pt2;
    __uint128_t lg1;
    __uint128_t lg2;
    __uint128_t f1;
    __uint128_t f2;
} SecureCREC;

typedef struct Beaver_triple{
    __uint128_t a1;
    __uint128_t a2;
    __uint128_t b1;
    __uint128_t b2;
    __uint128_t c1;
    __uint128_t c2;
}TRIPLE;


typedef struct co_Beaver_triple{
    __uint128_t a1;
    __uint128_t a2;
    __uint128_t b1[200];
    __uint128_t b2[200];
    __uint128_t c1[200];
    __uint128_t c2[200];
}COTRIPLE;



int write_header=0; 
int verbose = 2; 
int seed = 0;
int use_unk_vec = 1; 
int num_threads = 8; 
int num_iter = 25; 
int vector_size = 50; 
int save_gradsq = 0; 
int use_binary = 0; 
int model = 2; 
int checkpoint_every = 0; 
int load_init_param = 0; 
int save_init_param = 0; 
int load_init_gradsq = 0; 
real eta = 0.07; 
real real_eta=0.07;
real alpha = 1.0, x_max = 128.0; 
real grad_clip_value = 100.0; 
real *W, *gradsq, *cost;
long long num_lines, *lines_per_thread, vocab_size,lines_actual;
char vocab_file[MAX_STRING_LENGTH];
char input_file[MAX_STRING_LENGTH];
char save_W_file[MAX_STRING_LENGTH];
char save_gradsq_file[MAX_STRING_LENGTH];
char init_param_file[MAX_STRING_LENGTH];
char init_gradsq_file[MAX_STRING_LENGTH];
char maxword1[MAX_STRING_LENGTH];
char maxword2[MAX_STRING_LENGTH];
double max_lines=0.0;
long long err_count=0;
__uint128_t* WW1,*WW2,*lastw1,*lastw2,*gradsqu1,*gradsqu2;
long long communication=0;
const __uint128_t T=((__uint128_t)1)<<22; 
const __uint128_t T_=(((__uint128_t)1)<<106);
 
__uint128_t r_eta;
   __uint128_t t1;
pthread_mutex_t* mutexs;
pthread_mutex_t akkk;

     
         __uint128_t a;
    __uint128_t b;
    __uint128_t c;
    __uint128_t a1;
    __uint128_t a2;

    __uint128_t b1;
    __uint128_t b2;

    __uint128_t c1;
    __uint128_t c2;     

/**
 * Loads a save file for use as the initial values for the parameters or gradsq
 * Return value: 0 if success, -1 if fail
 */

__uint128_t llrand() {
    __uint128_t r ;
    
  
    for (int i = 0; i < 9; ++i) {
        r = (r << 15) | (rand() & 0x7FFF);
    }
    

    return r & (__uint128_t)(0-1);
    
}
__uint128_t truncate(__uint128_t num)
{
    __uint128_t res=0;
    res =(__uint128_t )(num/T);
    if((__int128_t) num<0)
    {
        res =res - T_;
    }
    return res;
}

void mpcMulti(__uint128_t x1, __uint128_t y1 ,__uint128_t x2, __uint128_t y2,__uint128_t* res0,__uint128_t* res1,TRIPLE* tri, int i){
 
   
      
 

    __uint128_t e1;
    __uint128_t e2;

    __uint128_t f1;
    __uint128_t f2;

    __uint128_t f;
    __uint128_t e;
 

    e1=x1-tri[i].a1;
    e2=x2-tri[i].a2;

    f1=y1-tri[i].b1;
    f2=y2-tri[i].b2;

    f =f1 +f2;
    e =e1 +e2;

    *(res0)=truncate(f*x1+e*y1+tri[i].c1);
    *(res1)=truncate(f*x2+e*y2+tri[i].c2-e*f);
}
void mpccoMulti(__uint128_t x1, __uint128_t y1 ,__uint128_t x2, __uint128_t y2,__uint128_t* res0,__uint128_t* res1,COTRIPLE* tri, int i,int j){
 
   
      
 

    __uint128_t e1;
    __uint128_t e2;

    __uint128_t f1;
    __uint128_t f2;

    __uint128_t f;
    __uint128_t e;
 

    e1=x1-tri[i].a1;
    e2=x2-tri[i].a2;

    f1=y1-tri[i].b1[j];
    f2=y2-tri[i].b2[j];

    f =f1 +f2;
    e =e1 +e2;

    *(res0)=truncate(f*x1+e*y1+tri[i].c1[j]);
    *(res1)=truncate(f*x2+e*y2+tri[i].c2[j]-e*f);
}

int load_init_file(char *file_name, real *array, long long array_size) {
    FILE *fin;
    long long a;
    fin = fopen(file_name, "rb");
    if (fin == NULL) {
        log_file_loading_error("init file", file_name);
        return -1;
    }
    for (a = 0; a < array_size; a++) {
        if (feof(fin)) {
            fprintf(stderr, "EOF reached before data fully loaded in %s.\n", file_name);
            fclose(fin);
            return -1;
        }
        fread(&array[a], sizeof(real), 1, fin);
    }
    fclose(fin);
    return 0;
}

void initialize_parameters() {

        a=llrand();
     b=llrand();
     c=a*b;
     a1=llrand();
    a2=a-a1;

     b1=llrand();
     b2=b-b1;

    c1=llrand();
    c2=c-c1;     
    t1=llrand();
    
    if (seed == 0) {
        seed = time(0);
    }
    fprintf(stderr, "Using random seed %d\n", seed);
    srand(seed);
    long long a;
    long long W_size = 2 * vocab_size * (vector_size + 1); 

    
    a = posix_memalign((void **)&W, 128, W_size * sizeof(real)); 
    
    a = posix_memalign((void **)&mutexs, 128, 2*vocab_size * sizeof(pthread_mutex_t));
    for(a=0;a<2*vocab_size;a++){
        int succ=pthread_mutex_init(&mutexs[a], NULL);
        if(succ!=0){
            printf("failed\n");
        }
    }
    pthread_mutex_init(&akkk,NULL);
    if (W == NULL) {
        fprintf(stderr, "Error allocating memory for W\n");
        exit(1);
    }
    a = posix_memalign((void **)&WW1, 128, W_size * sizeof(__uint128_t));
    

    


    a = posix_memalign((void **)&gradsq, 128, W_size * sizeof(real));
     
    a = posix_memalign((void **)&gradsqu1, 128, W_size * sizeof(__uint128_t));
    a = posix_memalign((void **)&gradsqu2, 128, W_size * sizeof(__uint128_t));
    a = posix_memalign((void **)&WW2, 128, W_size * sizeof(__uint128_t));
    if (gradsq == NULL) {
        fprintf(stderr, "Error allocating memory for gradsq\n");
        free(W);
        exit(1);
    }
    if (load_init_param) {
        
        fprintf(stderr, "\nLoading initial parameters from %s \n", init_param_file);
        if (load_init_file(init_param_file, W, W_size)) {
            free(W);
            free(gradsq);
            exit(1);
        }
    } else {
        
        for (a = 0; a < W_size; ++a) {
            
            WW1[a] = (__uint128_t)(long long)(((rand() / (real)RAND_MAX - 0.5) / (2*vector_size))*T);
            WW2[a]=(__uint128_t)(long long)(((rand() / (real)RAND_MAX - 0.5) / (2*vector_size))*T);
            
        }
    }

    if (load_init_gradsq) {
        
        fprintf(stderr, "\nLoading initial squared gradients from %s \n", init_gradsq_file);
        if (load_init_file(init_gradsq_file, gradsq, W_size)) {
            free(W);
            free(gradsq);
            exit(1);
        }
    } else {
        
        for (a = 0; a < W_size; ++a) {
            gradsqu1[a] = (__uint128_t)0; 
            gradsqu2[a] = (__uint128_t)0; 

        }
    }
    r_eta=(__uint128_t)(long long)(eta*T);
}

inline real check_nan(real update) {
    if (isnan(update) || isinf(update)) {
        fprintf(stderr,"\ncaught NaN in update");
        return 0.;
    } else {
        return update;
    }
}
double resotre(__uint128_t a,__uint128_t b){
    return (double)(__int128_t)(a+b)/T;
}
void Bit_X(__uint128_t x, bool*res){

    __uint128_t tmp;
    __uint128_t a=1;
    for(int i=0;i<128;i++){

        tmp=x&1;
        if(tmp==1){
            *(res+i)=true;
        }
        else{
            *(res+i)=false;
        }
        x=x>>1;

    }
}

void *glove_thread(void *vid) {
    FILE*fp;
    int MAX_TRIPLE=1000000;
    TRIPLE* Triples=(TRIPLE*)malloc(MAX_TRIPLE * sizeof(TRIPLE));
     COTRIPLE* CoTriples=(COTRIPLE*)malloc(MAX_TRIPLE * sizeof(COTRIPLE));
  char triplefile[50]="GenTriples";
   char cotriplefile[50]="GenCoTriples";
   fp = fopen(triplefile, "rb");
     fread(Triples, sizeof(TRIPLE), MAX_TRIPLE, fp);

        fp = fopen(cotriplefile, "rb");
     fread(CoTriples, sizeof(COTRIPLE), MAX_TRIPLE1, fp);


    

    
    
    
    
    

    
    

    
    

    
    int triple_count=0;
    int cobeaver_count=0;

    long long a, b ,l1, l2,mutex1,mutex2;
    long long id = *(long long*)vid;
    long long line_count=0;
    long long last_count=0;

    real diff,fdiff;
    SecureCREC cr;
    real fx;

    __uint128_t diff0,diff1, fdiff0,fdiff1, temp10,temp11,temp20,temp21, res0,res1,fx1,fx2,local_eta;
    FILE *fin;


    fin = fopen(input_file, "rb");


    __uint128_t t,crl1,crl2;
    if (fin == NULL) {
        
        log_file_loading_error("input file", input_file);
        pthread_exit(NULL);
    }
    fseeko(fin, (num_lines / num_threads * id) * (sizeof(SecureCREC)), SEEK_SET); 
    cost[id] = 0;
    
    for (a = 0; a < lines_per_thread[id]; a++) {
            

        fread(&cr, sizeof(SecureCREC), 1, fin);
        
        line_count+=1;
        if (cr.word1 < 1 || cr.word2 < 1) { continue; }
        mutex1=cr.word1-1;
        mutex2=(cr.word2-1)+vocab_size;
          
        pthread_mutex_lock(&mutexs[mutex1]);
        pthread_mutex_lock(&mutexs[mutex2]);
        
        
        crl1=cr.lg1;
        crl2=cr.lg2;
     if(line_count-last_count>10000)
    	{
         real_eta = eta * (1 - lines_actual / (real)(num_iter * num_lines + 1));
      if (real_eta < eta * 0.0001) real_eta = eta * 0.0001;

      lines_actual+=line_count-last_count;
      last_count=line_count;
      r_eta=(__uint128_t)(long long)(real_eta*T);
      printf("%cAlpha: %f", 13, real_eta);
        fflush(stdout);
      }
      local_eta=r_eta;
        if (feof(fin)) break;

     
 
            
        l1 = (cr.word1 - 1) * (vector_size + 1); 
        l2 = ((cr.word2 - 1) + vocab_size) * (vector_size + 1); 
  

        
        diff1 = (__uint128_t)0;
        diff0 = (__uint128_t)0;
        fdiff1 = (__uint128_t)0;
        fdiff0 = (__uint128_t)0;
        
        for (b = 0; b < vector_size; b++){
            
            mpcMulti(WW1[b+l1],WW1[b+l2],WW2[b+l1],WW2[b+l2],&res0,&res1,Triples,triple_count);
            triple_count+=1;
            if(triple_count==MAX_TRIPLE)
                triple_count=0;
            
                

            diff0+=res0;
            diff1+=res1;
            
        } 
        diff0 += WW1[vector_size + l1] + WW1[vector_size + l2] - crl1;
        diff1 += WW2[vector_size + l1] + WW2[vector_size + l2] - crl2;
        
        
        
        
        
        
        
        
        fx1=cr.f1;
        fx2=cr.f2;
        mpcMulti(diff0,fx1,diff1,fx2,&res0,&res1,Triples,triple_count);
         triple_count+=1;
            if(triple_count==MAX_TRIPLE)
                triple_count=0;

        fdiff0=res0;
        fdiff1=res1;
        
        

            
        
        fdiff=((double)(__int128_t  )(fdiff0+fdiff1))/T;
        diff=((double)(__int128_t )(diff0+diff1))/T;
        cost[id] += 0.5 * fdiff * diff; 
   

        for (b = 0; b < vector_size; b++) {
            
            mpccoMulti(fdiff0,WW1[b+l2],fdiff1,WW2[b+l2],&res0,&res1,CoTriples,cobeaver_count,b);
            
            
            temp10=truncate(res0*local_eta);
            temp11=truncate(res1*local_eta);
            
            mpccoMulti(fdiff0,WW1[b+l1],fdiff1,WW2[b+l1],&res0,&res1,CoTriples,cobeaver_count,b+100);
            
            
            temp20=truncate(res0*local_eta);
            temp21=truncate(res1*local_eta);


              WW1[b + l1] -= temp10;
              WW2[b+l1]   -= temp11;
              WW1[b + l2] -= temp20;
              WW2[b + l2] -= temp21;
           
        }
             cobeaver_count+=1;
            if(cobeaver_count==MAX_TRIPLE1)
                cobeaver_count=0;
        
        WW1[vector_size + l1] -= truncate(fdiff0*local_eta);
        WW2[vector_size + l1] -= truncate(fdiff1*local_eta);

        WW1[vector_size + l2] -= truncate(fdiff0*local_eta);
        WW2[vector_size + l2] -= truncate(fdiff1*local_eta);
        pthread_mutex_unlock(&mutexs[mutex1]);
        pthread_mutex_unlock(&mutexs[mutex2]);
    
        

        

        
    }

    
    fclose(fin);
    pthread_exit(NULL);
}


int save_params(int nb_iter) {
    /*
     * nb_iter is the number of iteration (= a full pass through the cooccurrence matrix).
     *   nb_iter  > 0 => checkpointing the intermediate parameters, so nb_iter is in the filename of output file.
     *   nb_iter == 0 => checkpointing the initial parameters
     *   else         => saving the final paramters, so nb_iter is ignored.
     */

    long long a, b;
    char format[20];
    char output_file[MAX_STRING_LENGTH+20], output_file_gsq[MAX_STRING_LENGTH+20];
    char *word = malloc(sizeof(char) * MAX_STRING_LENGTH + 1);
    if (NULL == word) {
        return 1;
    }
    FILE *fid, *fout;
    FILE *fgs = NULL;
    for (a = 0; a < 2 * vocab_size * (vector_size + 1); a++){

        W[a]=((double)(long long  )(WW1[a]+WW2[a]))/T; 

    }
    if (use_binary > 0 || nb_iter == 0) {
        
        
        if (nb_iter < 0)
            sprintf(output_file,"%s.bin",save_W_file);
        else
            sprintf(output_file,"%s.%03d.bin",save_W_file,nb_iter);

        fout = fopen(output_file,"wb");
        if (fout == NULL) {log_file_loading_error("weights file", save_W_file); free(word); return 1;}
        for (a = 0; a < 2 * vocab_size * (vector_size + 1); a++) fwrite(&W[a], sizeof(real), 1,fout);
        fclose(fout);
        if (save_gradsq > 0) {
            if (nb_iter < 0)
                sprintf(output_file_gsq,"%s.bin",save_gradsq_file);
            else
                sprintf(output_file_gsq,"%s.%03d.bin",save_gradsq_file,nb_iter);

            fgs = fopen(output_file_gsq,"wb");
            if (fgs == NULL) {log_file_loading_error("gradsq file", save_gradsq_file); free(word); return 1;}
            for (a = 0; a < 2 * vocab_size * (vector_size + 1); a++) fwrite(&gradsq[a], sizeof(real), 1,fgs);
            fclose(fgs);
        }
    }
    if (use_binary != 1) { 
        if (nb_iter < 0)
            sprintf(output_file,"%s.txt",save_W_file);
        else
            sprintf(output_file,"%s.%03d.txt",save_W_file,nb_iter);
        if (save_gradsq > 0) {
            if (nb_iter < 0)
                sprintf(output_file_gsq,"%s.txt",save_gradsq_file);
            else
                sprintf(output_file_gsq,"%s.%03d.txt",save_gradsq_file,nb_iter);

            fgs = fopen(output_file_gsq,"wb");
            if (fgs == NULL) {log_file_loading_error("gradsq file", save_gradsq_file); free(word); return 1;}
        }
        fout = fopen(output_file,"wb");
        if (fout == NULL) {log_file_loading_error("weights file", save_W_file); free(word); return 1;}
        fid = fopen(vocab_file, "r");
        sprintf(format,"%%%ds",MAX_STRING_LENGTH);
        if (fid == NULL) {log_file_loading_error("vocab file", vocab_file); free(word); fclose(fout); return 1;}
        if (write_header) fprintf(fout, "%lld %d\n", vocab_size, vector_size);
        for (a = 0; a < vocab_size; a++) {
            if (fscanf(fid,format,word) == 0) {free(word); fclose(fid); fclose(fout); return 1;}
            
            if (strcmp(word, "<unk>") == 0) {free(word); fclose(fid); fclose(fout);  return 1;}
            fprintf(fout, "%s",word);
            if (model == 0) { 
                for (b = 0; b < (vector_size + 1); b++) fprintf(fout," %lf", W[a * (vector_size + 1) + b]);
                for (b = 0; b < (vector_size + 1); b++) fprintf(fout," %lf", W[(vocab_size + a) * (vector_size + 1) + b]);
            }
            if (model == 1) 
                for (b = 0; b < vector_size; b++) fprintf(fout," %lf", W[a * (vector_size + 1) + b]);
            if (model == 2) 
                for (b = 0; b < vector_size; b++) fprintf(fout," %lf", W[a * (vector_size + 1) + b] + W[(vocab_size + a) * (vector_size + 1) + b]);
            fprintf(fout,"\n");
            if (save_gradsq > 0) { 
                fprintf(fgs, "%s",word);
                for (b = 0; b < (vector_size + 1); b++) fprintf(fgs," %lf", gradsq[a * (vector_size + 1) + b]);
                for (b = 0; b < (vector_size + 1); b++) fprintf(fgs," %lf", gradsq[(vocab_size + a) * (vector_size + 1) + b]);
                fprintf(fgs,"\n");
            }
            if (fscanf(fid,format,word) == 0) {
                
                fclose(fout);
                fclose(fid);
                free(word); 
                return 1;
                } 
        }

        if (use_unk_vec) {
            real* unk_vec = (real*)calloc((vector_size + 1), sizeof(real));
            real* unk_context = (real*)calloc((vector_size + 1), sizeof(real));
            strcpy(word, "<unk>");

            long long num_rare_words = vocab_size < 100 ? vocab_size : 100;

            for (a = vocab_size - num_rare_words; a < vocab_size; a++) {
                for (b = 0; b < (vector_size + 1); b++) {
                    unk_vec[b] += W[a * (vector_size + 1) + b] / num_rare_words;
                    unk_context[b] += W[(vocab_size + a) * (vector_size + 1) + b] / num_rare_words;
                }
            }

            fprintf(fout, "%s",word);
            if (model == 0) { 
                for (b = 0; b < (vector_size + 1); b++) fprintf(fout," %lf", unk_vec[b]);
                for (b = 0; b < (vector_size + 1); b++) fprintf(fout," %lf", unk_context[b]);
            }
            if (model == 1) 
                for (b = 0; b < vector_size; b++) fprintf(fout," %lf", unk_vec[b]);
            if (model == 2) 
                for (b = 0; b < vector_size; b++) fprintf(fout," %lf", unk_vec[b] + unk_context[b]);
            fprintf(fout,"\n");

            free(unk_vec);
            free(unk_context);
        }

        fclose(fid);
        fclose(fout);
        if (save_gradsq > 0) fclose(fgs);
    }
    free(word);
    return 0;
}


int train_glove() {
    long long a, file_size;
    int save_params_return_code;
    int b;
    FILE *fin;
    real total_cost = 0;
    int start,end;
    fprintf(stderr, "TRAINING MODEL\n");
    __uint128_t res0,res1,u08,t08,ua,ub;

    fprintf(stderr,"%f\n",(double)(long long )(res0+res1)/T);
    fin = fopen(input_file, "rb");
    if (fin == NULL) {log_file_loading_error("cooccurrence file", input_file); return 1;}
    fseeko(fin, 0, SEEK_END);
    file_size = ftello(fin);
    num_lines = file_size/(sizeof(SecureCREC)); 
    fclose(fin);
    fprintf(stderr,"Read %lld lines.\n", num_lines);
    lines_actual=0;
    if (verbose > 1) fprintf(stderr,"Initializing parameters...");
    initialize_parameters();
    if (verbose > 1) fprintf(stderr,"done.\n");
    if (save_init_param) {
        if (verbose > 1) fprintf(stderr,"Saving initial parameters... ");
        save_params_return_code = save_params(0);
        if (save_params_return_code != 0)
            return save_params_return_code;
        if (verbose > 1) fprintf(stderr,"done.\n");
    }
    if (verbose > 0) fprintf(stderr,"vector size: %d\n", vector_size);
    if (verbose > 0) fprintf(stderr,"vocab size: %lld\n", vocab_size);
    if (verbose > 0) fprintf(stderr,"x_max: %lf\n", x_max);
    if (verbose > 0) fprintf(stderr,"alpha: %lf\n", alpha);
    pthread_t *pt = (pthread_t *)malloc(num_threads * sizeof(pthread_t));
    lines_per_thread = (long long *) malloc(num_threads * sizeof(long long));
    
    time_t rawtime;
    struct tm *info;
    char time_buffer[80];
        start=time(NULL);

    
    for (b = 0; b < num_iter; b++) {
        total_cost = 0;
        for (a = 0; a < num_threads - 1; a++) lines_per_thread[a] = num_lines / num_threads;
        lines_per_thread[a] = num_lines / num_threads + num_lines % num_threads;
        long long *thread_ids = (long long*)malloc(sizeof(long long) * num_threads);
        for (a = 0; a < num_threads; a++) thread_ids[a] = a;
        for (a = 0; a < num_threads; a++) pthread_create(&pt[a], NULL, glove_thread, (void *)&thread_ids[a]);
        for (a = 0; a < num_threads; a++) pthread_join(pt[a], NULL);
        for (a = 0; a < num_threads; a++) total_cost += cost[a];
        free(thread_ids);

        time(&rawtime);
        info = localtime(&rawtime);
        strftime(time_buffer,80,"%x - %I:%M.%S%p", info);
        fprintf(stderr, "%s, iter: %03d, cost: %lf\n", time_buffer,  b+1, total_cost/num_lines);

        if (checkpoint_every > 0 && (b + 1) % checkpoint_every == 0) {
            fprintf(stderr,"    saving intermediate parameters for iter %03d...", b+1);
            save_params_return_code = save_params(b+1);
            if (save_params_return_code != 0) {
                free(pt);
                free(lines_per_thread);
                return save_params_return_code;
            }
            fprintf(stderr,"done.\n");
        }

    }
     end=time(NULL);
            printf("calculate use %f mins\n",(double)(end-start)/(60));
           
             
    free(pt);
    free(lines_per_thread);
    return save_params(-1);
}

int main(int argc, char **argv) {
    int i;
    FILE *fid;
    int result = 0;
    
    if (argc == 1) {
        printf("GloVe: Global Vectors for Word Representation, v0.2\n");
        printf("Author: Jeffrey Pennington (jpennin@stanford.edu)\n\n");
        printf("Usage options:\n");
        printf("\t-verbose <int>\n");
        printf("\t\tSet verbosity: 0, 1, or 2 (default)\n");
        printf("\t-write-header <int>\n");
        printf("\t\tIf 1, write vocab_size/vector_size as first line. Do nothing if 0 (default).\n");
        printf("\t-vector-size <int>\n");
        printf("\t\tDimension of word vector representations (excluding bias term); default 50\n");
        printf("\t-threads <int>\n");
        printf("\t\tNumber of threads; default 8\n");
        printf("\t-iter <int>\n");
        printf("\t\tNumber of training iterations; default 25\n");
        printf("\t-eta <float>\n");
        printf("\t\tInitial learning rate; default 0.05\n");
        printf("\t-alpha <float>\n");
        printf("\t\tParameter in exponent of weighting function; default 0.75\n");
        printf("\t-x-max <float>\n");
        printf("\t\tParameter specifying cutoff in weighting function; default 100.0\n");
        printf("\t-grad-clip\n");
        printf("\t\tGradient components clipping parameter. Values will be clipped to [-grad-clip, grad-clip] interval\n");
        printf("\t-binary <int>\n");
        printf("\t\tSave output in binary format (0: text, 1: binary, 2: both); default 0\n");
        printf("\t-model <int>\n");
        printf("\t\tModel for word vector output (for text output only); default 2\n");
        printf("\t\t   0: output all data, for both word and context word vectors, including bias terms\n");
        printf("\t\t   1: output word vectors, excluding bias terms\n");
        printf("\t\t   2: output word vectors + context word vectors, excluding bias terms\n");
        printf("\t-input-file <file>\n");
        printf("\t\tBinary input file of shuffled cooccurrence data (produced by 'cooccur' and 'shuffle'); default cooccurrence.shuf.bin\n");
        printf("\t-vocab-file <file>\n");
        printf("\t\tFile containing vocabulary (truncated unigram counts, produced by 'vocab_count'); default vocab.txt\n");
        printf("\t-save-file <file>\n");
        printf("\t\tFilename, excluding extension, for word vector output; default vectors\n");
        printf("\t-gradsq-file <file>\n");
        printf("\t\tFilename, excluding extension, for squared gradient output; default gradsq\n");
        printf("\t-save-gradsq <int>\n");
        printf("\t\tSave accumulated squared gradients; default 0 (off); ignored if gradsq-file is specified\n");
        printf("\t-checkpoint-every <int>\n");
        printf("\t\tCheckpoint a  model every <int> iterations; default 0 (off)\n");
        printf("\t-load-init-param <int>\n");
        printf("\t\tLoad initial parameters from -init-param-file; default 0 (false)\n");
        printf("\t-save-init-param <int>\n");
        printf("\t\tSave initial parameters (i.e., checkpoint the model before any training); default 0 (false)\n");
        printf("\t-init-param-file <file>\n");
        printf("\t\tBinary initial parameters file to be loaded if -load-init-params is 1; (default is to look for vectors.000.bin)\n");
        printf("\t-load-init-gradsq <int>\n");
        printf("\t\tLoad initial squared gradients from -init-gradsq-file; default 0 (false)\n");
        printf("\t-init-gradsq-file <file>\n");
        printf("\t\tBinary initial squared gradients file to be loaded if -load-init-gradsq is 1; (default is to look for gradsq.000.bin)\n");
        printf("\t-seed <int>\n");
        printf("\t\tRandom seed to use.  If not set, will be randomized using current time.");
        printf("\nExample usage:\n");
        printf("./glove -input-file cooccurrence.shuf.bin -vocab-file vocab.txt -save-file vectors -gradsq-file gradsq -verbose 2 -vector-size 100 -threads 16 -alpha 0.75 -x-max 100.0 -eta 0.05 -binary 2 -model 2\n\n");
        result = 0;
    } else {
        if ((i = find_arg((char *)"-write-header", argc, argv)) > 0) write_header = atoi(argv[i + 1]);
        if ((i = find_arg((char *)"-verbose", argc, argv)) > 0) verbose = atoi(argv[i + 1]);
        if ((i = find_arg((char *)"-vector-size", argc, argv)) > 0) vector_size = atoi(argv[i + 1]);
        if ((i = find_arg((char *)"-iter", argc, argv)) > 0) num_iter = atoi(argv[i + 1]);
        if ((i = find_arg((char *)"-threads", argc, argv)) > 0) num_threads = atoi(argv[i + 1]);
        cost = malloc(sizeof(real) * num_threads);
        if ((i = find_arg((char *)"-alpha", argc, argv)) > 0) alpha = atof(argv[i + 1]);
        if ((i = find_arg((char *)"-x-max", argc, argv)) > 0) x_max = atof(argv[i + 1]);
        if ((i = find_arg((char *)"-eta", argc, argv)) > 0) eta = atof(argv[i + 1]);
        if ((i = find_arg((char *)"-grad-clip", argc, argv)) > 0) grad_clip_value = atof(argv[i + 1]);
        if ((i = find_arg((char *)"-binary", argc, argv)) > 0) use_binary = atoi(argv[i + 1]);
        if ((i = find_arg((char *)"-model", argc, argv)) > 0) model = atoi(argv[i + 1]);
        if (model != 0 && model != 1) model = 2;
        if ((i = find_arg((char *)"-save-gradsq", argc, argv)) > 0) save_gradsq = atoi(argv[i + 1]);
        if ((i = find_arg((char *)"-vocab-file", argc, argv)) > 0) strcpy(vocab_file, argv[i + 1]);
        else strcpy(vocab_file, (char *)"vocab.txt");
        if ((i = find_arg((char *)"-save-file", argc, argv)) > 0) strcpy(save_W_file, argv[i + 1]);
        else strcpy(save_W_file, (char *)"vectors");
        if ((i = find_arg((char *)"-gradsq-file", argc, argv)) > 0) {
            strcpy(save_gradsq_file, argv[i + 1]);
            save_gradsq = 1;
        }
        else if (save_gradsq > 0) strcpy(save_gradsq_file, (char *)"gradsq");
        if ((i = find_arg((char *)"-input-file", argc, argv)) > 0) strcpy(input_file, argv[i + 1]);
        else strcpy(input_file, (char *)"secure.shuf.bin");
        if ((i = find_arg((char *)"-checkpoint-every", argc, argv)) > 0) checkpoint_every = atoi(argv[i + 1]);
        if ((i = find_arg((char *)"-init-param-file", argc, argv)) > 0) strcpy(init_param_file, argv[i + 1]);
        else strcpy(init_param_file, (char *)"vectors.000.bin");
        if ((i = find_arg((char *)"-load-init-param", argc, argv)) > 0) load_init_param = atoi(argv[i + 1]);
        if ((i = find_arg((char *)"-save-init-param", argc, argv)) > 0) save_init_param = atoi(argv[i + 1]);
        if ((i = find_arg((char *)"-init-gradsq-file", argc, argv)) > 0) strcpy(init_gradsq_file, argv[i + 1]);
        else strcpy(init_gradsq_file, (char *)"gradsq.000.bin");
        if ((i = find_arg((char *)"-load-init-gradsq", argc, argv)) > 0) load_init_gradsq = atoi(argv[i + 1]);
        if ((i = find_arg((char *)"-seed", argc, argv)) > 0) seed = atoi(argv[i + 1]);
        
        vocab_size = 0;
        fid = fopen(vocab_file, "r");
        if (fid == NULL) {log_file_loading_error("vocab file", vocab_file); free(cost); return 1;}
        while ((i = getc(fid)) != EOF) if (i == '\n') vocab_size++; 
        fclose(fid);
        if (vocab_size == 0) {fprintf(stderr, "Unable to find any vocab entries in vocab file %s.\n", vocab_file); free(cost); return 1;}
        result = train_glove();
        free(cost);
    }
    free(W);
    free(gradsq);

    return result;
}
