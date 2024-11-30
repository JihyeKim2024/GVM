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
unsigned int s=5,q=2;
double t_max=1000000000.0;
double p=0.5;
uniform_real_distribution<double> urd(0,1);
uniform_int_distribution<int> random_nei(0,s-1);
vector<double> aver_et(3,0.0);
vector<double> aver_ep(3,0.0);
vector<double> exit_time(3,0.0); //total sum of exit time
vector<double> exit_count(3,0.0); //total sum of exit ensemble
vector<double> exit_count_p(3,0.0);
vector<vector<double>> exit_times(3);
FILE*fp;

char fname[100];
int main()
{
    unsigned int t_ensemble=1;
    for(int hyper_en=0; hyper_en<100000000; ++hyper_en)
    {
        for(int w=2; w<5; ++w)
        {
            unsigned int N=pow(10,w+3);
            double dt=1.0/N;
            uniform_int_distribution<int> random_node(0,N-1);
            vector<bool> state(N);
            unsigned int sum_n=0;
            for(int i=(int)N*0.5; i<N; ++i)
            {
                state[i]=1;
                ++sum_n;
                
            }//initialization
            if(sum_n==N)
            {
                exit_times[w-2].emplace_back(0);
                ++exit_count_p[w-2];
                ++exit_count[w-2];
                aver_ep[w-2]=(double)exit_count_p[w-2]/t_ensemble;
                aver_et[w-2]=(double)exit_time[w-2]/exit_count[w-2];
            }
            else if(sum_n==0)
            {
                exit_times[w-2].emplace_back(0);
                ++exit_count[w-2];
                aver_et[w-2]=(double)exit_time[w-2]/exit_count[w-2];
                aver_ep[w-2]=(double)exit_count_p[w-2]/t_ensemble;
            }
            else
            {
                double t_save=0.0;
                for(double t=dt; t<=t_max; t+=dt)
                {
                    set<int> hyperedge;
                    vector<int> hyperedge_vec;
                    while(hyperedge.size()<s)
                    {
                        unsigned int nei=random_node(gen);
                        if(hyperedge.emplace(nei).second)
                        {
                            hyperedge_vec.emplace_back(nei);
                        }
                        
                    }
                    int count=0;
                    set<int> members;
                    while(members.size()<q)
                    {
                        unsigned int selected_nei=hyperedge_vec[random_nei(gen)];
                        if(members.emplace(selected_nei).second)
                        {
                            count+=state[selected_nei];
                            
                        }
                    }
                    if(count==q or count==0)
                    {
                        int del_state=0;
                        for(auto & mem: hyperedge_vec)
                        {
                            bool pre_state=state[mem];
                            state[mem]=count/q;
                            del_state+=state[mem]-pre_state;
                        }
                        sum_n+=del_state;
                        
                    }
                    //update step
                    if(sum_n==N or sum_n==0)
                    {
                        t_save=t;
                        break;
                    }
                }//t
                if(sum_n==N)
                {
                    exit_times[w-2].emplace_back(t_save);
                    exit_time[w-2]+=t_save;
                    ++exit_count_p[w-2];
                    ++exit_count[w-2];
                    aver_et[w-2]=(double)exit_time[w-2]/exit_count[w-2];
                    aver_ep[w-2]=(double)exit_count_p[w-2]/t_ensemble;
                }
                else if(sum_n==0)
                {
                    exit_times[w-2].emplace_back(t_save);
                    exit_time[w-2]+=t_save;
                    ++exit_count[w-2];
                    aver_et[w-2]=(double)exit_time[w-2]/exit_count[w-2];
                    aver_ep[w-2]=(double)exit_count_p[w-2]/t_ensemble;
                }
            }
        }//w
        if(t_ensemble==10 or t_ensemble==50 or t_ensemble==100 or t_ensemble==1000 or t_ensemble==5000 or t_ensemble==10000 or t_ensemble==50000 or t_ensemble==100000 or t_ensemble==1000000 or t_ensemble==10000000 or t_ensemble==100000000 or t_ensemble==1000000000)
        {
            sprintf(fname,"FigS6_s5_q%d_en%d.txt",q,t_ensemble);
            fp=fopen(fname,"w");
            for(int i=0; i<3; ++i)
            {
                unsigned int NN=pow(10,i+5);
                double sum_Ti=accumulate(exit_times[i].begin(),exit_times[i].end(),0.0);
                double sq_sum=inner_product(exit_times[i].begin(),exit_times[i].end(),exit_times[i].begin(),0.0);
                double stdev=sqrt((sq_sum-2*sum_Ti*aver_et[i]+t_ensemble*aver_et[i]*aver_et[i])/(t_ensemble-1));
                fprintf(fp,"%d\t%.9f\t%.9f\t%.9f\n",NN,aver_ep[i],aver_et[i],stdev);
            }
            fclose(fp);
        }
        ++t_ensemble;
    }//making a hypergraph
}

