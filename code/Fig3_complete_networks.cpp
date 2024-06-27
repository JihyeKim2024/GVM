#define _CRT_SECURE_NO_WARNINGS
#include "Random.h"
#include<iostream>
#include<vector>
#include<cmath>
#include<set>
#include<algorithm>
using namespace std;
random_device rd;
splitmix gen(rd);
int N = 10000, q = 2, t_ensemble = 1;
double t_max=1000000000.0, dt=1.0/N;
uniform_real_distribution<double> urd(0,1);
uniform_int_distribution<int> random_node(0,N-1);
vector<double> aver_ep(11,0.0);
vector<double> exit_count_p(11,0.0);
FILE*fp;
char fname[100];

int main()
{
    for(int hyper_en=0; hyper_en<10000000; ++hyper_en)
    {
        for(int w=0; w<11; ++w)
        {
            vector<bool> state(N,1);
            double p=0.4+w*0.02; //rho(0)
            int sum_n=N;
            for(int i=(int)N*p; i<N; ++i)
            {
                state[i]=0;
                --sum_n;
            }//initialization
            if(sum_n==N)
            {
                ++exit_count_p[w];
                aver_ep[w]=(double)exit_count_p[w]/t_ensemble;
            }
            else if(sum_n==0)
            {
                aver_ep[w]=(double)exit_count_p[w]/t_ensemble;
            }
            else
            {
                for(double t=dt; t<=t_max; t+=dt)
                {
                    int voter=random_node(gen);
                    bool state_of_voter=state[voter];
                    int count=0;
                    for(int pick_count=0; pick_count<q; ++pick_count)
                    {
                        int selected_nei=random_node(gen);
                        while (selected_nei == voter)
                            selected_nei = random_node(gen);
                        if (state[selected_nei] != state_of_voter)
                            ++count;
                        else
                            break;
                    }
                    if(count==q)
                    {
                        state[voter]=1-state[voter];
                        int del_state=2*state[voter]-1;
                        sum_n+=del_state;
                    }//update step
                    if(sum_n==N or sum_n==0)
                    {
                        break;
                    }
                }// t
                if(sum_n==N)
                {
                    ++exit_count_p[w];
                    aver_ep[w]=(double)exit_count_p[w]/t_ensemble;
                }
                else if(sum_n==0)
                {
                    aver_ep[w]=(double)exit_count_p[w]/t_ensemble;
                }
            }
        }// rho(0)
        if(t_ensemble==1 or t_ensemble==50 or t_ensemble==100 or t_ensemble==1000 or t_ensemble==5000 or t_ensemble==10000 or t_ensemble==50000 or t_ensemble==100000 or t_ensemble==500000 or t_ensemble==1000000)
        {
            sprintf(fname,"Fig3_sN_q%d_en%d.txt",q,t_ensemble);
            fp=fopen(fname,"w");
            for(int i=0; i<11; ++i)
            {
                fprintf(fp,"%.4f\t%.9f\n",0.4+i*0.02,aver_ep[i]);
            }
            fclose(fp);
        }
        ++t_ensemble;
    }
}

