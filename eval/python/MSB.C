#include <iostream>
#include <cmath>
#include <random>
#include <bitset>

using namespace std;
#define DATALENGTH 64 //modular
const unsigned long T= pow(2,16); //T is the length of fraction part
const unsigned long T_=pow(2,48); //t_ is the truncation param
mt19937_64 gen (random_device{}()); //gen is the 64-bit random generator


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

}

void mpcMsb(long long x1, long long x2, bool *msb) {
    bitset<DATALENGTH>  Bitx1(x1), Bitx2(x2);
    bool result0=false,result1=false;
    bool gParty1[DATALENGTH]={false}, pParty1[DATALENGTH]={false}, gParty2[DATALENGTH]={false}, pParty2[DATALENGTH]={false},
            x1Bits[DATALENGTH]={false}, x2Bits[DATALENGTH]={false}, bitx[DATALENGTH - 1]={false}, res[2 * DATALENGTH]={false};
    for (int i = 0; i < DATALENGTH; ++i) {
        x1Bits[i] = Bitx1[i];
        x2Bits[i] = Bitx2[i];
    }

    for (int i = 0; i < DATALENGTH - 1; ++i) {
        mpcMultiBool(x1Bits[i], bitx[i], bitx[i], x2Bits[i],&result0,&result1);
        res[i]=result0;
        res[DATALENGTH - 1 + i]=result1;

    }


    gParty1[0] = false;
    gParty2[0] = false;

    for (int i = 0; i < DATALENGTH - 1; ++i) {
        gParty1[i + 1] = res[i];
        gParty2[i + 1] = res[DATALENGTH - 1 + i];

        pParty1[i + 1] = x1Bits[i];
        pParty2[i + 1] = x2Bits[i];
    }

    int length = DATALENGTH;

    for (int i = 0; i < log2(DATALENGTH); ++i) {
        mpcMultiBool(pParty1[1], gParty1[0], pParty2[1], gParty2[0],&result0,&result1);
        gParty1[0] = result0 ^ gParty1[1];
        gParty2[0] = result1 ^ gParty2[1];

        int num = 1;
        for (int j = 2; j < length; j +=2) {
            mpcMultiBool(pParty1 [j+1], gParty1 [j], pParty2 [j+1], gParty2[j],&result0,&result1);

            gParty1[num] = result0 ^ gParty1[j + 1];
            gParty2[num] = result1 ^ gParty2[j + 1];

            mpcMultiBool(pParty1 [j+1], pParty1 [j], pParty2 [j+1], pParty2 [j],&result0,&result1);

            pParty1[num] = result0;
            pParty2[num++] = result1;
        }
        length /= 2;
    }

    msb[0] = x1Bits[DATALENGTH - 1] ^ gParty1[0];
    msb[1] = x2Bits[DATALENGTH - 1] ^ gParty2[0];
}


//truncate number after multiple
unsigned long truncate(unsigned long num)
{
    unsigned long res=0;
    res =round(num/T);
    if((long)num<0)
    {
        res =res - T_;
    }
    return res;
}

int main() {
    // test for calculate msb(x)
    double x=-5.6321;
    /6 x
    unsigned long ul_x= (unsigned long)(x*T); unsigned long x1=gen(); unsigned long x2=ul_x-x1;
    //test for calculate msb(x)
    bool msb[2];
    mpcMsb((long)x1,(long)x2,msb);
    cout<<(msb[0]^msb[1])<<endl;
}