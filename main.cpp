#include <iostream>
#include <math.h>
#include <random>
#include <time.h>
#include <chrono>
#include <fstream>
#include <unistd.h>
#include <string.h>
#include "parameters.h"

using namespace std::chrono;
using namespace std;

double Rand() {
    thread_local std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> dis(0, 1);
    return dis(gen);
}

int mod(int a,int b){
    return (a%b + b)%b;
}
void shiftx(int b[], int a[], int sign, int v){
    for (int i=0; i<d; i++){
        if (i==v){
            if (v==0){ b[i]=mod(a[i]+sign,Nt);  }
            else { b[i]=mod(a[i]+sign,Ns);  }
        }
        else { b[i]=a[i]; }
    }
}
void random(int a[]){
    for (int i=1; i<d; i++){
        a[i]=rand()%Ns;
    }
    a[0]=rand()%Nt;
}

long double I(int s){
    long double a=0,r=0;
    while(r<inf){
        a+=dr*pow(r,s+1)*exp(-eta*pow(r,2)-lmd*pow(r,4));
        r+=dr;
    }
    return a;
}

bool eql(int x1[],int x2[]){
    bool flag=true;
    for (int i=0;i<d;i++){
        if(x1[i]!=x2[i]){
            flag=false;
            break;
        }
    }
    return flag;
}

int ksum(int k[][Ns][d]){
    int sum=0;
    for (int i=0;i<Nt;i++){
        for (int j=0;j<Ns;j++){
            sum+=k[i][j][0];
        }
    }
    return sum;
}
long int sx(int x[], int k[][Ns][d], int a[][Ns][d]){
    int sum=0;
    int v[d]={0};
    for (int i=0;i<d;i++){
        v[i]=1;

        sum+=abs(k[x[0]][x[1]][i])
        +abs(k[mod(x[0]-v[0],Nt)][mod(x[1]-v[1],Ns)][i])
        +2*(a[x[0]][x[1]][i]
        +a[mod(x[0]-v[0],Nt)][mod(x[1]-v[1],Ns)][i]);

        v[i]=0;
    }
    return sum;
}

void update(){
    int y;
    for (int i=0;i<Nt;i++){
        for (int j=0;j<Ns;j++){
            for (int t=0;t<d;t++){
                
                y=a_[i][j][t];
                a_[i][j][t]=a[i][j][t]+2*(rand()%2)-1;
                
                if (a_[i][j][t]<0){
                    a_[i][j][t]=y;
                    continue;
                }
                int temp[d], temp1[d];
                temp[0]=i; temp[1]=j;
                shiftx(temp1, temp, 1, t);
                if (a_[i][j][t]>a[i][j][t]){
                    
                    rho=1.0/(abs(k[i][j][t])+a_[i][j][t])
                    /a_[i][j][t]
                    *I_val[sx(temp,k,a_)]/I_val[sx(temp,k,a)]
                    *I_val[sx(temp1,k,a_)]/I_val[sx(temp1,k,a)];
                }
                else {
                    
                    rho=1.0*(abs(k[i][j][t])+a[i][j][t])
                    *a[i][j][t]
                    *I_val[sx(temp,k,a_)]/I_val[sx(temp,k,a)]
                    *I_val[sx(temp1,k,a_)]/I_val[sx(temp1,k,a)];
                }
                if (Rand()<rho){
                    a[i][j][t]=a_[i][j][t];
                }
                else{
                    a_[i][j][t]=y;
                }
            }
        }
    }
    int x2[2], x1[2], x12[2];
    for (int i=0; i<Nt; i++){
        for (int j=0; j<Ns; j++){
            
            del=2*(rand()%2)-1;
            x[0]=i; x[1]=j;
            rho=1.0;
            shiftx(x2,x,1,1);
            shiftx(x1,x,1,0);
            shiftx(x12,x2,1,0);
            if (del>0){
                if (k[x[0]][x[1]][1]>=0) rho*=1.0/(k[x[0]][x[1]][1]+1+a[x[0]][x[1]][1]);
                else rho*=abs(k[x[0]][x[1]][1])+a[x[0]][x[1]][1];
                if (k[x2[0]][x2[1]][0]>=0) rho*=1.0/(k[x2[0]][x2[1]][0]+1+a[x2[0]][x2[1]][0]);
                else rho*=abs(k[x2[0]][x2[1]][0])+a[x2[0]][x2[1]][0];
                if (k[x[0]][x[1]][0]>0) rho*=k[x[0]][x[1]][0]+a[x[0]][x[1]][0];
                else rho*=1.0/(abs(k[x[0]][x[1]][0])+1+a[x[0]][x[1]][0]);
                if (k[x1[0]][x1[1]][1]>0) rho*=k[x1[0]][x1[1]][1]+a[x1[0]][x1[1]][1];
                else rho*=1.0/(abs(k[x1[0]][x1[1]][1])+1+a[x1[0]][x1[1]][1]);
                
                if (k[x[0]][x[1]][1]>=0 && k[x[0]][x[1]][0]<=0) rho*=I_val[sx(x,k,a)+2]/I_val[sx(x,k,a)];
                else if (k[x[0]][x[1]][1]<0 && k[x[0]][x[1]][0]>0) rho*=I_val[sx(x,k,a)-2]/I_val[sx(x,k,a)]; 
                if (k[x2[0]][x2[1]][0]>=0 && k[x1[0]][x1[1]][1]<=0) rho*=I_val[sx(x12,k,a)+2]/I_val[sx(x12,k,a)];
                else if (k[x2[0]][x2[1]][0]<0 && k[x1[0]][x1[1]][1]>0) rho*=I_val[sx(x12,k,a)-2]/I_val[sx(x12,k,a)];
                if (k[x2[0]][x2[1]][0]>=0 && k[x[0]][x[1]][1]>=0) rho*=I_val[sx(x2,k,a)+2]/I_val[sx(x2,k,a)];
                else if (k[x2[0]][x2[1]][0]<0 && k[x[0]][x[1]][1]<0) rho*=I_val[sx(x2,k,a)-2]/I_val[sx(x2,k,a)];
                if (k[x1[0]][x1[1]][1]<=0 && k[x[0]][x[1]][0]<=0) rho*=I_val[sx(x1,k,a)+2]/I_val[sx(x1,k,a)];
                else if(k[x1[0]][x1[1]][1]>0 && k[x[0]][x[1]][0]>0) rho*=I_val[sx(x1,k,a)-2]/I_val[sx(x1,k,a)];
            }
            else {
                if (k[x[0]][x[1]][1]<=0) rho*=1.0/(abs(k[x[0]][x[1]][1])+1+a[x[0]][x[1]][1]);
                else rho*=k[x[0]][x[1]][1]+a[x[0]][x[1]][1];
                if (k[x2[0]][x2[1]][0]<=0) rho*=1.0/(abs(k[x2[0]][x2[1]][0])+1+a[x2[0]][x2[1]][0]);
                else rho*=k[x2[0]][x2[1]][0]+a[x2[0]][x2[1]][0];
                if (k[x[0]][x[1]][0]<0) rho*=abs(k[x[0]][x[1]][0])+a[x[0]][x[1]][0];
                else rho*=1.0/(k[x[0]][x[1]][0]+1+a[x[0]][x[1]][0]);
                if (k[x1[0]][x1[1]][1]<0) rho*=abs(k[x1[0]][x1[1]][1])+a[x1[0]][x1[1]][1];
                else rho*=1.0/(k[x1[0]][x1[1]][1]+1+a[x1[0]][x1[1]][1]);
                
                if (k[x[0]][x[1]][1]<=0 && k[x[0]][x[1]][0]>=0) rho*=I_val[sx(x,k,a)+2]/I_val[sx(x,k,a)];
                else if (k[x[0]][x[1]][1]>0 && k[x[0]][x[1]][0]<0) rho*=I_val[sx(x,k,a)-2]/I_val[sx(x,k,a)]; 
                if (k[x2[0]][x2[1]][0]<=0 && k[x1[0]][x1[1]][1]>=0) rho*=I_val[sx(x12,k,a)+2]/I_val[sx(x12,k,a)];
                else if (k[x2[0]][x2[1]][0]>0 && k[x1[0]][x1[1]][1]<0) rho*=I_val[sx(x12,k,a)-2]/I_val[sx(x12,k,a)];
                if (k[x2[0]][x2[1]][0]<=0 && k[x[0]][x[1]][1]<=0) rho*=I_val[sx(x2,k,a)+2]/I_val[sx(x2,k,a)];
                else if (k[x2[0]][x2[1]][0]>0 && k[x[0]][x[1]][1]>0) rho*=I_val[sx(x2,k,a)-2]/I_val[sx(x2,k,a)];
                if (k[x1[0]][x1[1]][1]>=0 && k[x[0]][x[1]][0]>=0) rho*=I_val[sx(x1,k,a)+2]/I_val[sx(x1,k,a)];
                else if(k[x1[0]][x1[1]][1]<0 && k[x[0]][x[1]][0]<0) rho*=I_val[sx(x1,k,a)-2]/I_val[sx(x1,k,a)];
            }
            
            if (Rand()<rho){
                k[x[0]][x[1]][1]+=del;
                k[x2[0]][x2[1]][0]+=del;
                k[x1[0]][x1[1]][1]-=del;
                k[x[0]][x[1]][0]-=del;
                
            }
        }
    }
    
    for (int i=0; i<Nt; i++){
        del=2*(rand()%2)-1;
        rho=1.0;
        
        for (int j=0; j<Ns; j++){
            x[0]=i; x[1]=j;
            shiftx(x_,x,-1,1);
            if (del>0){
                if (k[i][j][1]>=0) rho*=1.0/(k[i][j][1]+1+a[i][j][1]);
                else rho*=abs(k[i][j][1])+a[i][j][1];
                if (k[x[0]][x[1]][1]>=0 && k[x_[0]][x_[1]][1]>=0) rho*=I_val[sx(x,k,a)+2]/I_val[sx(x,k,a)];
                else if (k[x[0]][x[1]][1]<0 && k[x_[0]][x_[1]][1]<0) rho*=I_val[sx(x,k,a)-2]/I_val[sx(x,k,a)];
            }
            else {
                if (k[i][j][1]>0) rho*=k[i][j][1]+a[i][j][1];
                else rho*=1.0/(abs(k[i][j][1])+a[i][j][1]+1);
                if (k[x[0]][x[1]][1]<=0 && k[x_[0]][x_[1]][1]<=0) rho*=I_val[sx(x,k,a)+2]/I_val[sx(x,k,a)];
                else if (k[x[0]][x[1]][1]>0 && k[x_[0]][x_[1]][1]>0) rho*=I_val[sx(x,k,a)-2]/I_val[sx(x,k,a)];
            }
        }
        
        if (Rand()<rho){
            for (int j=0; j<Ns; j++){
                k[i][j][1]+=del;
            }
        }
    }
    
    for (int j=0; j<Ns; j++){
        del=2*(rand()%2)-1;
        rho=1.0;
        
        for (int i=0; i<Nt; i++){
            x[0]=i; x[1]=j;
            shiftx(x_,x,-1,0);
            if (del>0){
                if (k[i][j][0]>=0) rho*=1.0/(k[i][j][0]+1+a[i][j][0]);
                else rho*=abs(k[i][j][0])+a[i][j][0];
                
                if (k[x[0]][x[1]][0]>=0 && k[x_[0]][x_[1]][0]>=0) rho*=I_val[sx(x,k,a)+2]/I_val[sx(x,k,a)];
                else if (k[x[0]][x[1]][0]<0 && k[x_[0]][x_[1]][0]<0) rho*=I_val[sx(x,k,a)-2]/I_val[sx(x,k,a)];
            }
            else {
                if (k[i][j][0]>0) rho*=k[i][j][0]+a[i][j][0];
                else rho*=1.0/(abs(k[i][j][0])+a[i][j][0]+1);
                
                if (k[x[0]][x[1]][0]<=0 && k[x_[0]][x_[1]][0]<=0) rho*=I_val[sx(x,k,a)+2]/I_val[sx(x,k,a)];
                else if (k[x[0]][x[1]][0]>0 && k[x_[0]][x_[1]][0]>0) rho*=I_val[sx(x,k,a)-2]/I_val[sx(x,k,a)];
            }
            rho*=exp(mu*del);
        }
        
        if (Rand()<rho){
            for (int i=0; i<Nt; i++){
                k[i][j][0]+=del;
            }
        }
    }
    
}
float phi2(int k[][Ns][d], int a[][Ns][d], long double I[int_val]){
    long double sum=0;
    int x[d];
    for (int i=0;i<Nt;i++){
        for (int j=0;j<Ns;j++){
            x[0]=i; x[1]=j;
            sum+=I[sx(x, k, a)+2]/I[sx(x, k, a)];
        }
    }
    return sum;
}
float phi4(int k[][Ns][d], int a[][Ns][d], long double I[int_val]){
    long double sum=0;
    int x[d];
    for (int i=0;i<Nt;i++){
        for (int j=0;j<Ns;j++){
            x[0]=i; x[1]=j;
            sum+=I[sx(x, k, a)+4]/I[sx(x, k, a)];
        }
    }
    return sum;
}
float errorjack(float xi[configs]){
    float x_i[configs], x_=0, stddev=0;
    for (int i=0; i<configs; i++){
        for (int j=0; j<configs; j++){
            x_i[i]=(1-(i==j))*xi[j];
        }
        x_i[i]=x_i[i]/(configs-1);
        x_+=x_i[i];
    }
    x_=x_/configs;
    for (int i=0; i<configs; i++){
        stddev+=(x_i[i]-x_)*(x_i[i]-x_);
    }
    stddev=sqrt(stddev*(configs-1)/configs);
    return stddev;
}
int main(int argc, char **argv)
{
    
    auto begin=high_resolution_clock::now();
    srand(time(NULL));
    for (int i=0; i<int_val; i++){
        I_val[i]=I(i);
    }
    float phi2_avg, phi4_avg;
    long double dmu=(mu_max-mu_min)/mu_n;
    float xi[configs],phi2i[configs];
    mu=mu_min;
    
    ofstream data, data1;
    string filename="mu_n_phi2_phi4_Nt"+to_string(Nt)+"_Ns"+to_string(Ns)+"_eta"
        +to_string(eta)+"_lmd"+to_string(lmd)+"_"+to_string(configs)+"_"+to_string(gaps)+"_"
        +to_string(equil)+".txt";
    //data.open(filename);
    data.open("mu_vs_rho.txt");
    data1.open("mu_vs_phi2.txt");
    //data2.open("mu_vs_phi4.txt");
    
    for (int g=0; g<mu_n; g++){
        
        for (int i=0; i<equil; i++){
            update();
            
        }
        phi2_avg=0;
        n_avg=0;
        //phi4_avg=0;
        for (int i=0; i<configs; i++){
            for (int j=0; j<gaps; j++){
                update();
                
            }
            update();
            
            xi[i]=1.0*ksum(k)/Nt/Ns;
            phi2i[i]=phi2(k,a,I_val);
            
            n_avg+=xi[i];
            phi2_avg+=phi2i[i];
            //phi4_avg+=phi4(k,a,I_val);
        }
        n_avg=n_avg/configs;
        phi2_avg=phi2_avg/configs;
        //phi4_avg=phi4_avg/configs;
        
        data<<mu<<"\t"<<n_avg<<"\t"<<errorjack(xi);
        data1<<"\t"<<phi2_avg<<"\t"<<errorjack(phi2i);
        
        mu+=dmu;
        cout<<g<<endl;
        
    }
    
    data.close();
    data1.close();
    //data2.close();
    auto stop=high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop-begin);
    cout<<duration.count()<<endl;
    
    return 0;
}