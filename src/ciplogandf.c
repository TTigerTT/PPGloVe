#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include<stdbool.h>
#include <pthread.h>
#include<math.h>
#include <time.h>
typedef double real;


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

typedef struct tpc_seccooccur{
    int word1;
    int word2;
    __uint128_t pt1;
    __uint128_t pt2;
}Sec_CREC;

typedef struct Beaver_triple{
    __uint128_t a1;
    __uint128_t a2;
    __uint128_t b1;
    __uint128_t b2;
    __uint128_t c1;
    __uint128_t c2;
}TRIPLE;

typedef struct cooccur_rec {
    int word1;
    int word2;
    real val;
} CREC;
const __uint128_t T=((__uint128_t)1)<<22; //T is the length of fraction part
const __uint128_t T_=(((__uint128_t)1)<<106); //t_ is the truncation param
__uint128_t *w1,*w2,*s1,*s2,*f1,*f2;
int frac_num=27;
TRIPLE* tri;
int MAX_TRIPLE=1000000;
int triple_count=0;

__uint128_t base[30],inverse[30];
long long num_lines,num_threads=50,communication=0;
long long *lines_per_thread;
int global=0;
    __uint128_t a;
    __uint128_t b;
    __uint128_t c;
    __uint128_t a1;
    __uint128_t a2;

    __uint128_t b1;
    __uint128_t b2;

    __uint128_t c1;
    __uint128_t c2; 


__uint128_t llrand() {
    __uint128_t r ;
    // srand((unsigned)time(NULL));
  
    for (int i = 0; i < 9; ++i) {
        r = (r << 15) | (rand() & 0x7FFF);
    }
    // printf("%llu\n",r);

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

double realf(double a){
    if(a>=10)
        return 1.0;
    else{
        return a/10;
    }

}

void mpcMultiBool(bool x1, bool y1,bool x2, bool y2, bool *result0, bool *result1) {

    bool a1, b1, c1, a2, b2, c2, e, f;
    c1 = true;
    b1 = true;
    a1 = false;
    c2 = true;
    b2 = true;
    a2 = false;

    e = x1  ^ a1 ^ x2  ^ a2;
    f = y1  ^ b1 ^ y2  ^ b2;

    (*result0) = f & a1 ^ e & b1 ^ c1;
    (*result1) = e & f ^f & a2 ^ e & b2 ^ c2;

    // communication+=4;


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

void MSB(__uint128_t x1, __uint128_t x2, bool* msb){
    bool x1_bits[128]={false}, x2_bits[128]={false},g_party1[128]={false},g_party2[128]={false},p_party1[128]={false},p_party2[128]={false},bit_x[128]={false};
    bool res0=false,res1=false;
    int length=128,num=1;
    Bit_X(x1,x1_bits);
    Bit_X(x2,x2_bits);
    for(int i=0;i<128;i++){

        mpcMultiBool(x1_bits[i],bit_x[i],bit_x[i],x2_bits[i],&res0,&res1);
        g_party1[i]=res0;
        g_party2[i]=res1;
        p_party1[i]=x1_bits[i];
        p_party2[i]=x2_bits[i];
    }

    for(int j=1;j<127;j+=2){
        mpcMultiBool(g_party1[j],p_party1[j+1],g_party2[j],p_party2[j+1],&res0,&res1);
        
        g_party1[num]=res0^g_party1[j+1];
        g_party2[num]=res1^g_party2[j+1];

        mpcMultiBool(p_party1[j],p_party1[j+1],p_party2[j],p_party2[j+1],&res0,&res1);
        p_party1[num]=res0;
        p_party2[num]=res1;
        
        num++;
    }
    length=length/2;

    for(int i=2;i<=7;i++){

        num=0;

        for(int k=0;k<length;k+=2){
            mpcMultiBool(g_party1[k],p_party1[k+1],g_party2[k],p_party2[k+1],&res0,&res1);
            g_party1[num]=res0^g_party1[k+1];
            g_party2[num]=res1^g_party2[k+1];

             mpcMultiBool(p_party1[k],p_party1[k+1],p_party2[k],p_party2[k+1],&res0,&res1);
             p_party1[num]=res0;
             p_party2[num]=res1;
        
             num++;
        }
        length=length/2;

    }
    *(msb+0)=p_party1[127]^g_party1[0];
    *(msb+1)=true^p_party2[127]^g_party2[0];


}
void mpcMulti(__uint128_t x1, __uint128_t y1 ,__uint128_t x2, __uint128_t y2,__uint128_t* res0,__uint128_t* res1 ){
  

    __uint128_t e1;
    __uint128_t e2;

    __uint128_t f1;
    __uint128_t f2;

    __uint128_t f;
    __uint128_t e;
    int i=triple_count;
    triple_count+=1;
    if(triple_count>=MAX_TRIPLE)
        triple_count=0;
    

     e1=x1-tri[i].a1;
     e2=x2-tri[i].a2;

     f1=y1-tri[i].b1;
     f2=y2-tri[i].b2;

     f =f1 +f2;
     e =e1 +e2;

    *(res0)=(f*x1+e*y1+tri[i].c1);
    *(res1)=(f*x2+e*y2+tri[i].c2-e*f);

    // communication+=4*128;
 
}

void B2A(bool x1,bool x2,__uint128_t *res0,__uint128_t *res1){

    __uint128_t e0,f0,e1,f1,t0,t1;
    if(x1){
        e0=1;
    }
    else{
        e0=0;
    }
    f0=0;
    e1=0;
    if(x2){
        f1=1;
    }
    else{
        f1=0;
    }
    mpcMulti(e0,f0,e1,f1,&t0,&t1);
    *(res0)=e0+f0-2*t0;
    *(res1)=e1+f1-2*t1;
}

void div2n(__uint128_t x0,__uint128_t x1,__uint128_t* base,__uint128_t * inverse,
int frac_num,__uint128_t *res0,__uint128_t *res1,__uint128_t *n0,__uint128_t *n1){

    __uint128_t res=0,x0_d,x1_d,b2a0,b2a1,acheck0[30],acheck1[30],tmp0,tmp1,t_res0=0,t_res1=0,t_n0=0,t_n1=0,num=-4;
    bool msb[2]={false};
    bool check0[30]={false};
    bool check1[30]={false};
    for(int i=0;i<frac_num-1;i++){
        x0_d=x0;
        x1_d=x1-(*(base+i));
        MSB(x0_d,x1_d,msb);
        check0[i]=msb[0];
        check1[i]=msb[1];
    }
    // check0[0]=false;
    // check1[0]=true;
    check0[frac_num-1]=false;
    check1[frac_num-1]=false;

    for(int i=frac_num-1;i>0;i--){
        check0[i]=check0[i]^check0[i-1];
        check1[i]=check1[i]^check1[i-1];
    }
    check0[0]=check0[0]^true;

    for(int i=0;i<frac_num;i++){
        B2A(check0[i],check1[i],&b2a0,&b2a1);
        acheck0[i]=b2a0;
        acheck1[i]=b2a1;
    }
    for(int i=0;i<frac_num;i++){
        t_res0+=acheck0[i]*(*(inverse+i));
        t_res1+=acheck1[i]*(*(inverse+i));
        t_n0+=acheck0[i]*num;
        t_n1+=acheck1[i]*num;
        num+=1;
    }

    mpcMulti(x0,t_res0,x1,t_res1,&tmp0,&tmp1);
    tmp0=truncate(tmp0);
    tmp1=truncate(tmp1);
    *(res0)=tmp0;
    *(res1)=tmp1;
    *(n0)=t_n0;
    *(n1)=t_n1;
}
void logarithm_app(__uint128_t x0, __uint128_t x1,__uint128_t n0,__uint128_t n1,__uint128_t* res0, __uint128_t* res1){
    double ln2=log(2);
    __uint128_t u_ln2=(__uint128_t )(ln2*T);
    __uint128_t x0_sub=x0-T;
    __uint128_t x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,x80,x81,x90,x91,x100,x101,x110,x111,x120,x121,x130,x131,x140,x141,x150,x151,x160,x161;
    __uint128_t inverse12=(__uint128_t) ((double)1/840*T);
    __uint128_t a=(__uint128_t )(1*T),b=(__uint128_t )((double)1/2*T),c=(__uint128_t )((double)1/3*T),d=(__uint128_t )((double)1/4*T),e=(__uint128_t )((double)1/5*T) ;
    __uint128_t f=(__uint128_t )((double)1/6*T),g=(__uint128_t )((double)1/7*T),h=(__uint128_t )((double)1/8*T);
    __uint128_t h9=(__uint128_t )((double)1/9*T),h10=(__uint128_t )((double)1/10*T),h11=(__uint128_t )((double)1/11*T),h12=(__uint128_t )((double)1/12*T),h13=(__uint128_t )((double)1/13*T),h14=(__uint128_t )((double)1/14*T),h15=(__uint128_t )((double)1/15*T),h16=(__uint128_t )((double)1/16*T);
    // printf("%llu,%f\n",x0+x1,(double)(x0_sub+x1)/T);
    // printf("%llu\n",n0+n1);
    mpcMulti(x0_sub,x0_sub,x1,x1,&x20,&x21);
    x20=truncate(x20);
    x21=truncate(x21);
    // printf("%llu,%f,%f",x0_sub+x1,(double)(x20+x21)/T);
    mpcMulti(x20,x0_sub,x21,x1,&x30,&x31);

    x30=truncate(x30);
    x31=truncate(x31);
    mpcMulti(x20,x20,x21,x21,&x40,&x41);
    x40=truncate(x40);
    x41=truncate(x41);
    mpcMulti(x40,x0_sub,x41,x1,&x50,&x51);
    mpcMulti(x40,x20,x41,x21,&x60,&x61);
    x50=truncate(x50);
    x51=truncate(x51);
    x60=truncate(x60);
    x61=truncate(x61);
     mpcMulti(x40,x30,x41,x31,&x70,&x71);
    mpcMulti(x40,x40,x41,x41,&x80,&x81);
    x70=truncate(x70);
    x71=truncate(x71);
 
    
    x80=truncate(x80);
    x81=truncate(x81);

    mpcMulti(x0_sub,x80,x1,x81,&x90,&x91);
    mpcMulti(x80,x20,x81,x21,&x100,&x101);
    mpcMulti(x80,x30,x81,x31,&x110,&x111);
    mpcMulti(x80,x40,x81,x41,&x120,&x121);
    mpcMulti(x80,x50,x81,x51,&x130,&x131);
    mpcMulti(x80,x60,x81,x61,&x140,&x141);
    mpcMulti(x80,x70,x81,x71,&x150,&x151);
    mpcMulti(x80,x80,x81,x81,&x160,&x161);

    x90=truncate(x90);
    x91=truncate(x91);
    x100=truncate(x100);
    x101=truncate(x101);
    x110=truncate(x110);
    x111=truncate(x111);
    x120=truncate(x120);
    x121=truncate(x121);

    x130=truncate(x130);
    x131=truncate(x131);
    x140=truncate(x140);
    x141=truncate(x141);
    x150=truncate(x150);
    x151=truncate(x151);
    x160=truncate(x160);
    x161=truncate(x161);




    *(res0)=n0*u_ln2+truncate(a*x0_sub)-truncate(b*x20)+truncate(c*x30)-truncate(d*x40)+truncate(e*x50)-truncate(f*x60)+truncate(g*x70)-truncate(h*x80)+truncate(h9*x90)-truncate(h10*x100)+truncate(h11*x110)-truncate(h12*x120)+truncate(h13*x130)-truncate(h14*x140)+truncate(h15*x150)-truncate(h16*x160);
    *(res1)=n1*u_ln2+truncate(a*x1)-truncate(b*x21)+truncate(c*x31)-truncate(d*x41)+truncate(e*x51)-truncate(f*x61)+truncate(g*x71)-truncate(h*x81)+truncate(h9*x91)-truncate(h10*x101)+truncate(h11*x111)-truncate(h12*x121)+truncate(h13*x131)-truncate(h14*x141)+truncate(h15*x151)-truncate(h16*x161);

    // *(res0)=n0*u_ln2+truncate(a*x0_sub)-truncate(b*x20)+truncate(c*x30)-truncate(d*x40)+truncate(e*x50)-truncate(f*x60)+truncate(g*x70)-truncate(h*x80);
    // *(res1)=n1*u_ln2+truncate(a*x1)-truncate(b*x21)+truncate(c*x31)-truncate(d*x41)+truncate(e*x51)-truncate(f*x61)+truncate(g*x71)-truncate(h*x81);
// *(res0)=n0*u_ln2+truncate(a*x0_sub)-truncate(b*x20)+truncate(c*x30)-truncate(d*x40)+truncate(e*x50)-truncate(f*x60)+truncate(g*x70)-truncate(h*x80)+truncate(h9*x90)-truncate(h10*x100)+truncate(h11*x110)-truncate(h12*x120);
//     *(res1)=n1*u_ln2+truncate(a*x1)-truncate(b*x21)+truncate(c*x31)-truncate(d*x41)+truncate(e*x51)-truncate(f*x61)+truncate(g*x71)-truncate(h*x81)+truncate(h9*x91)-truncate(h10*x101)+truncate(h11*x111)-truncate(h12*x121);

}
void *fun(__uint128_t d1,__uint128_t d2,__uint128_t* res0,__uint128_t* res1){
        __uint128_t d1plus=d1-(__uint128_t)10*T;
        __uint128_t a1,a2,b1,b2,c1,c2,f1,f2,g1,g2;
        bool msb[2]={false,false};
        __uint128_t xmax=(__uint128_t )((double)1/10*T);
        __uint128_t part2=T;
        a1=truncate(d1*xmax);
        a2=truncate(d2*xmax);
        MSB(d1plus,d2,msb);
        B2A(msb[0],msb[1],&b1,&b2);
        B2A(msb[0],true^msb[1],&c1,&c2);
        f1=part2*b1;
        f2=part2*b2;
        mpcMulti(a1,c1,a2,c2,&g1,&g2);
        *(res0)=f1+g1;
        *(res1)=f2+g2;        

}


void *log_traindata(void* vid){
 
    __uint128_t res0,res1,n0,n1;
    __uint128_t t;
    long long  numas;
    long long a;
     long long id = *(long long*)vid;
 

    numas=num_lines/num_threads*id;
    

     for (a = 0; a < lines_per_thread[id]; a++) {
        global+=1;
        if(global%10000000==0){
            fprintf(stderr,"get\n");
        }

             div2n(s1[numas+a],s2[numas+a],base,inverse,frac_num,&res0,&res1,&n0,&n1);
             logarithm_app(res0,res1,n0,n1,&res0,&res1);

            w1[numas+a]=res0;
            w2[numas+a]=res1;
         
            fun(s1[numas+a],s2[numas+a],&res0,&res1);
           
            f1[numas+a]=res0;
            f2[numas+a]=res1;
}
     pthread_exit(NULL);
}

int main(){
     a=llrand();
     b=llrand();
     c=a*b;
     a1=llrand();
     a2=a-a1;

     b1=llrand();
     b2=b-b1;

     c1=llrand();
     c2=c-c1;  

    FILE*fp;
  tri=(TRIPLE*)malloc(MAX_TRIPLE * sizeof(TRIPLE));

  char triplefile[50]="GenTriples";

   fp = fopen(triplefile, "rb");
     fread(tri, sizeof(TRIPLE), MAX_TRIPLE, fp);




    bool res33[128],msb[2];
        FILE *fin,*fout;
       for(int i=0;i<frac_num;i++){
        base[i]=(__uint128_t)1<<(i+22-4);
        inverse[i]=(__uint128_t)1<<(22-i+4);
    }

     int start,end;
    __uint128_t t,res0,res1,n0,n1,check;


    fin = fopen("mergeshuff_20mb.bin","rb");
    Sec_CREC cr;
    SecureCREC scr;

    fseeko(fin, 0, SEEK_END);
    long long file_size = ftello(fin);
    num_lines = file_size/((sizeof(Sec_CREC))); 
    printf("%lld\n",num_lines);
    long long i=0;
    long long a;
    fclose(fin);
    
    w1=(__uint128_t *) malloc(((int) num_lines)*sizeof(__uint128_t));
    w2=(__uint128_t *) malloc(((int) num_lines)*sizeof(__uint128_t));
    s1=(__uint128_t *) malloc(((int) num_lines)*sizeof(__uint128_t));
    s2=(__uint128_t *) malloc(((int) num_lines)*sizeof(__uint128_t));
    f1=(__uint128_t *) malloc(((int) num_lines)*sizeof(__uint128_t));
    f2=(__uint128_t *) malloc(((int) num_lines)*sizeof(__uint128_t));
   
    fin = fopen("mergeshuff_20mb.bin","rb");
    printf("begin");
    for(i=0;i<num_lines;i++){
         fread(&cr, sizeof(Sec_CREC), 1, fin);
            if (feof(fin)) break;
             s1[i]=cr.pt1;
             s2[i]=cr.pt2;
            // printf("%f %f\n",cr.val,cr.logval);
    if(i%100000000==0){
                   fprintf(stderr,"done.\n");

    }
    }
    fclose(fin);


    
   
        pthread_t *pt = (pthread_t *)malloc(num_threads * sizeof(pthread_t));
    lines_per_thread = (long long *) malloc(num_threads * sizeof(long long));
    fprintf(stderr,"begin the numlines %lld\n",num_lines);

        start=time(NULL);
        for (a = 0; a < num_threads - 1; a++) lines_per_thread[a] = num_lines / num_threads;
        lines_per_thread[a] = num_lines / num_threads + num_lines % num_threads;
        long long *thread_ids = (long long*)malloc(sizeof(long long) * num_threads);
        for (a = 0; a < num_threads; a++) thread_ids[a] = a;
        for (a = 0; a < num_threads; a++) pthread_create(&pt[a], NULL, log_traindata, (void *)&thread_ids[a]);
        for (a = 0; a < num_threads; a++) pthread_join(pt[a], NULL);
        free(thread_ids);

        end=time(NULL);

    printf("calculate use %f mins\n",(double)(end-start)/(60));

    int t0_num=0;
     double avger=0;
     for(int i=0;i<num_lines;i++){
        double resu=(double)(long long  )(w1[i]+w2[i])/T;
        double real_log=log((double)(long long  )(s1[i]+s2[i])/T);
      if(real_log!=0){
        avger+=fabs(resu-real_log)/fabs(real_log);
      }
      else{
        
        t0_num+=1;
      }

     }

     avger=avger/num_lines;
    printf("The average error is %f",avger);
    printf("The 0_num is %d",t0_num);
     fin = fopen("mergeshuff_20mb.bin","rb");
     fout=fopen("secmergeshuf_20series.bin","wb");
     for(i=0;i<num_lines;i++){
             fread(&cr, sizeof(Sec_CREC), 1, fin);
            if (feof(fin)) break;
             scr.word1=cr.word1;
             scr.word2=cr.word2;
             scr.pt1=cr.pt1;
             scr.pt2=cr.pt2;
             scr.lg1=w1[i];
             scr.lg2=w2[i];
             scr.f1=f1[i];
             scr.f2=f2[i];
            // printf("%d %d ",(int)(scr.pt1+scr.pt2),(int)(scr.f1+scr.f2));
            fwrite(&scr,sizeof(SecureCREC),1,fout);
        if(i%10000000==0)
            printf("%lld\n",i);
     }

    fclose(fout);

    



    }

