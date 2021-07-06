

float eta=4.01,lmd=1.0;   //theory parameters
const int d=2, Ns=10, Nt=10;    //lattice size and dimensions
float dr=0.01,inf=10; //simulation parameters
const int configs=4000,gaps=10,equil=2000;    //simulation parameters
const int mu_n=20, int_val=3000;
long double rho, I_val[int_val], mu_min=0.0, mu_max=1.0, mu;

float n_avg;
int del, v, sign;

int x[d], x_[d], xx[d], px[4][d], x0[d]
    , k_[Nt][Ns][d]={0}, k1[Nt][Ns][d]={0}
    , k[Nt][Ns][d]={0}, a[Nt][Ns][d]={0}
    , a_[Nt][Ns][d]={0};
