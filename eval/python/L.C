#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
 void myitoa(__uint128_t v, char* s)
 {
      char temp;
      int i=0, j;
  
     while(v >0) {
       
         s[i++] = v % 10 + '0';
         v /= 10;    }
     s[i] = '\0'; 
     j=0;
     i--;
     while(j < i) {        temp = s[j];
         s[j] = s[i];
         s[i] = temp;
         j++;
         i--;
     }
 }
unsigned long long llrand() {
    unsigned long long r = 0;

    for (int i = 0; i < 5; ++i) {
        r = (r << 15) | (rand() & 0x7FFF);
    }

    return r & 0xFFFFFFFFFFFFFFFFULL;
}
int main(){
    
//     __uint128_t base=(__uint128_t)1<<24;
//    __uint128_t M=((__uint128_t)1<<59);
//    __uint128_t  k=(__uint128_t)2<<48;
//   __uint128_t  x=base;
//   char s[40];
//     for(int i=0;i<20;i++){
      
//         x=(x*(k-((M*x)>>48)))>>48;
//         myitoa(x,s);
//       printf("%s\n",s);
//     }
//     x=x>>36;
//     myitoa(x,s);
//     printf("%s\n",s);

   time_t t;
   unsigned long long div[24];
   unsigned long long n[24];
   unsigned long long ran_num;

srand((unsigned) time(&t));

for( int i=0; i< 24;i++){
    div[i]=((unsigned long long)1)<<(30-i-1);
    n[i]=(unsigned long long )0;
}
n[3]=(unsigned long long ) 1;
unsigned long long p1=0, p2=0;
for( int i=0; i<24;i++){
    ran_num=llrand();
    unsigned long long t1= n[i]-ran_num;
    unsigned long long t2= ran_num;
    p1+=t1*div[i];
    p2+=t2*div[i];
    
}
unsigned long long  k=(unsigned long long) 19;

unsigned long long k1=llrand();
unsigned long long k2=k-k1;

unsigned long long a=llrand();
unsigned long long b=llrand();
unsigned long long c=a*b;
unsigned long long a1=llrand();
unsigned long long a2=a-a1;

unsigned long long b1=llrand();
unsigned long long b2=b-b1;

unsigned long long c1=llrand();
unsigned long long c2=c-c1;

unsigned long long e1=k1-a1;
unsigned long long e2=k2-a2;

unsigned long long f1=p1-b1;
unsigned long long f2=p2-b2;

unsigned long long f =f1 +f2;
unsigned long long e =e1 +e2;

unsigned long long result1=(f*k1+e*p1+c1);
unsigned long long result2=(f*k2+e*p2+c2-e*f);



printf("%llu\n",result1);
printf("%llu\n",result2);
printf("%llu\n",(result1+result2));


}