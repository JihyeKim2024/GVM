#define _CRT_SECURE_NO_WARNINGS
#include "Random.h"
#include<iostream>
#include<vector>
#include<cmath>
#include<set>
#include<algorithm>
#include <numeric>
using namespace std;
random_device rd;
splitmix gen(rd);
uniform_real_distribution<double> urd(0, 1);
uniform_int_distribution<int> random_nei(0, s - 2);
int s = 7, q = 2, t_ensemble = 1;
double t_max=1000000000.0, p=0.65;
vector<double> aver_et(3,0.0);
vector<double> aver_ep(3,0.0);
vector<double> exit_time(3,0.0); 
vector<double> exit_count(3,0.0); 
vector<double> exit_count_p(3,0.0);
FILE*fp;
char fname[100];

int main()
{
    for(int hyper_en=0; hyper_en<100000000; ++hyper_en)
    {
        for(int w=0; w<3; ++w)
        {
            int N=pow(10,w+2);
            double dt=1.0/N;
            uniform_int_distribution<int> random_node(0,N-1);
            vector<bool> state(N,1);
            int sum_n=N;
            for(int i=(int)N*p; i<N; ++i)
            {
                state[i]=0;
                --sum_n;
            }//initialization
            double t_save=0.0;
            for(double t=dt; t<=t_max; t+=dt)
            {
                int voter=random_node(gen);
                bool state_of_voter=state[voter];
                set<int> hyperedge;
                vector<int> hyperedge_vec;
                hyperedge.emplace(voter);
                while(hyperedge.size()<s)
                {
                    int nei=random_node(gen);
                    if(hyperedge.emplace(nei).second)
                    {
                        hyperedge_vec.emplace_back(nei);
                    }
                }
                int count=0;
                for(int pick_count=0; pick_count<q; ++pick_count)
                {
                    int selected_nei=hyperedge_vec[random_nei(gen)];
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
                    t_save=t;
                    break;
                }
            }//t
            if(sum_n==N)
            {
                exit_time[w]+=t_save;
                ++exit_count_p[w];
                ++exit_count[w];
                aver_et[w]=(double)exit_time[w]/exit_count[w];
                aver_ep[w]=(double)exit_count_p[w]/t_ensemble;
            }
            else if(sum_n==0)
            {
                exit_time[w]+=t_save;
                ++exit_count[w];
                aver_et[w]=(double)exit_time[w]/exit_count[w];
                aver_ep[w]=(double)exit_count_p[w]/t_ensemble;
            }
        }// N
        if(t_ensemble==10 or t_ensemble==50 or t_ensemble==100 or t_ensemble==1000 or t_ensemble==5000 or t_ensemble==10000 or t_ensemble==50000 or t_ensemble==100000 or t_ensemble==1000000 or t_ensemble==10000000 or t_ensemble==100000000 or t_ensemble==1000000000)
        {
            sprintf(fname,"FigS2_s%d_q%d_en%d.txt",s,q,t_ensemble);
            fp=fopen(fname,"w");
            for(int i=0; i<3; ++i)
            {
                unsigned int NN=pow(10,i+2);
                fprintf(fp,"%d\t%.9f\t%.9f\n",NN,aver_ep[i],aver_et[i]);
            }
            fclose(fp);
        }
        ++t_ensemble;
    }// ensemble
}
