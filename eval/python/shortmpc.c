#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include<stdbool.h>
#include<math.h>

const unsigned short T= (unsigned short )1<<8; //T is the length of fraction part
const unsigned short T_=(unsigned short )1<<8; //t_ is the truncation param

    unsigned short a;
    unsigned short b;
    unsigned short c;
    unsigned short a1;
    unsigned short a2;

    unsigned short b1;
    unsigned short b2;

    unsigned short c1;
    unsigned short c2;     

    unsigned short e1;
    unsigned short e2;

    unsigned short f1;
    unsigned short f2;

    unsigned short f;
    unsigned short e;

unsigned short llrand() {
    unsigned short r ;
    for (int i = 0; i < 2; ++i) {
        r = (r << 15) | (rand() & 0x7FFF);
    }
     return r & (unsigned short)0xFFFF;
}

unsigned short truncate(unsigned short num)
{
    unsigned short res=0;
    res =(unsigned short )((double)num/T);
    if((short int)num<0)
    {
        res =res - T_;
    }
    return res;
}
void mpcMulti(unsigned short x1, unsigned short y1 ,unsigned short x2, unsigned short y2,unsigned short* res0,unsigned short* res1 ){
    
   a=llrand();
     b=llrand();
     c=a*b;
     a1=llrand();
    a2=a-a1;

   b1=llrand();
     b2=b-b1;

    c1=llrand();
   c2=c-c1;     

  e1=x1-a1;
    e2=x2-a2;

 f1=y1-b1;
    f2=y2-b2;

    f =f1 +f2;
    e =e1 +e2;

    *(res0)=(f*x1+e*y1+c1);
    *(res1)=(f*x2+e*y2+c2-e*f);
}

int main(){
    unsigned short res0,res1;
     unsigned short m=(unsigned short)(0.8*T);
    unsigned short n=(unsigned short)(0.8*T);
   printf("%hx\n",m);
   unsigned short t=llrand();
   printf("%hx\n",t);
  
   unsigned short k0=m-t,k1=t;
    unsigned short z0=n-t,z1=t;
    unsigned short z=k0+k1;
      printf("%u\n",k0);
      printf("%u\n",k1);
        printf("%hx\n",z);
       mpcMulti(k0,z0,k1,z1,&res0,&res1);

    res0=truncate(res0);
    res1=truncate(res1);


        for(int i=0;i<1000;i++){ 
        mpcMulti(k0,k0,k1,k1,&res0,&res1);
        printf("%d %d \n",res0,res1);
        res0=truncate(res0);
        res1=truncate(res1);
        printf("%hx %hx %hx %hx",59511,47641,65512,65466);
        // if ((unsigned short)(res0+res1)!=41616&&(unsigned short)(res0+res1)!=42025){
        //     printf("%d %d\n ",res0,res1);
        //     printf("%hx %hx %hx %hx %hx %hx %hx %hx\n",k0,k1,a1,a2,b1,b2,c1,c2);
        //     printf("%u %u %u %u\n",(unsigned short)(a1+a2),(unsigned short)(b1+b2),(unsigned short)((a1+a2)*(b1+b2)),(unsigned short)(c1+c2));
        //     printf("%d %d %d %d %d %d %d %d\n",k0,k1,a1,a2,b1,b2,c1,c2);
        //      printf("%d\n",(unsigned short)(k0-a1));
        //       printf("%d\n",(unsigned short)(k1-a2));
        //        printf("%d\n",(unsigned short)(k0-b1));
        //         printf("%d\n",(unsigned short)(k1-b2));printf("%d %d ",res0,res1);;
        //          printf("%d\n",(unsigned short)((k0-a1)+k1-a2));
        //           printf("%d\n",(unsigned short)(k0-b1+k1-b2));


        //     }

    }
}