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
int N = 10000, s = 7, t_ensemble = 1, iq = 2, fq = 6, dq = fq - iq;
double t_max = 1000000000.0, dt = 1.0 / N;
uniform_int_distribution<int> random_node(0,N-1);
uniform_real_distribution<double> urd(0,1);
uniform_int_distribution<int> random_nei(0,s-2);
vector<double> aver_et(dq,0.0);
vector<double> aver_ep(dq,0.0);
vector<double> exit_time(dq,0.0); 
vector<double> exit_count(dq,0.0);
vector<double> exit_count_p(dq,0.0);
FILE*fp;
char fname[100];

int main()
{
    for(int hyper_en=0; hyper_en<100000000; ++hyper_en)
    {
        for(int q=iq; q<fq; ++q)
        {
            vector<bool> state(N,0);
            unsigned int sum_n=0;
            for(int i=(int)N*0.5; i<N; ++i)
            {
                    state[i]=1;
                    ++sum_n;
            }//initialization
            double t_save=0.0;
            for(double t=dt; t<=t_max; t+=dt)
            {
                unsigned int voter=random_node(gen);
                bool state_of_voter=state[voter];
                set<int> hyperedge;
                vector<int> hyperedge_vec;
                hyperedge.emplace(voter);
                while(hyperedge.size()<s)
                {
                    unsigned int nei=random_node(gen);
                    if(hyperedge.emplace(nei).second)
                    {
                        hyperedge_vec.emplace_back(nei);
                    }    
                }
                int count=0;
                for(int pick_count=0; pick_count<q; ++pick_count)
                {
                        unsigned int selected_nei=hyperedge_vec[random_nei(gen)];
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
            }// t
            if(sum_n==N)
            {
                exit_time[q-2]+=t_save;
                ++exit_count_p[q-2];
                ++exit_count[q-2];
                aver_et[q-2]=(double)exit_time[q-2]/exit_count[q-2];
                aver_ep[q-2]=(double)exit_count_p[q-2]/t_ensemble;
            }
            else if(sum_n==0)
            {
                exit_time[q-2]+=t_save;
                ++exit_count[q-2];
                aver_et[q-2]=(double)exit_time[q-2]/exit_count[q-2];
                aver_ep[q-2]=(double)exit_count_p[q-2]/t_ensemble;
            }

        }// q
        if(t_ensemble==1 or t_ensemble==50 or t_ensemble==100 or t_ensemble==1000 or t_ensemble==5000 or t_ensemble==10000 or t_ensemble==50000 or t_ensemble==100000 or t_ensemble==1000000 or t_ensemble==10000000 or t_ensemble==100000000 or t_ensemble==1000000000)
        {
            sprintf(fname,"qdep_completes%d_en%d.txt",s,t_ensemble);
            fp=fopen(fname,"w");
            for(int i=0; i<dq; ++i)
            {
                fprintf(fp,"%d\t%.9f\t%.9f\n",i+2,aver_ep[i],aver_et[i]);
            }
            fclose(fp);
        }
        ++t_ensemble;
    }// ensemble
}

