#define _CRT_SECURE_NO_WARNINGS
#include "Random.h"
#include<iostream>
#include<vector>
#include<cmath>
#include<set>
#include<algorithm>
#include <numeric>
#include<map>
using namespace std;
random_device rd;
splitmix gen(rd);
unsigned int s=5,q=2;
double N=10000000;
double dt=1.0/N, gam=3.0;
double t_max=1000000000.0;
double p=0.5;
uniform_real_distribution<double> urd(0,1);
uniform_int_distribution<int> random_member(0,s-2);
uniform_int_distribution<int> random_initial(0,N-1);
double aver_et=0.0;
double aver_ep=0.0;
double exit_time=0.0; //total sum of exit time
double exit_count=0.0; //total sum of exit ensemble
double exit_count_p=0.0;
vector<double> exit_times;
vector<double> q_vector;
vector<double> updated_qvector;
vector<vector<int>> y_vector(N);
vector<int> poor_vec;
vector<int> rich_vec;
double denom=0;
FILE*fp;
char fname[100];
int main()
{
    unsigned int t_ensemble=1;
    denom=0.0;
    for(double i=N; i>0; i-=1)
    {
        denom+=pow(i,-gam);
    }
    for(int i=1; i<N+1; ++i)
    {
        double Pk=pow(i,-gam)/denom;
        q_vector.emplace_back(N*Pk);
        updated_qvector.emplace_back(N*Pk);
        if(q_vector[i-1]>1){rich_vec.emplace_back(i);}
        else if(q_vector[i-1]<1){poor_vec.emplace_back(i);}
    }
    while(poor_vec.size()!=0)
    {
        int pooper=poor_vec[0];
        int richer=rich_vec[0];
        y_vector[pooper-1].emplace_back(richer);
        updated_qvector[richer-1]-=(1-q_vector[pooper-1]);
        poor_vec.erase(poor_vec.begin());
        if(updated_qvector[richer-1]<1)
        {
            rich_vec.erase(rich_vec.begin());
            poor_vec.emplace_back(richer);
        }
    }
    for(int hyper_en=0; hyper_en<100000000; ++hyper_en)
    {
        double sum_degree=0;
        vector<int> node_pocket;
        vector<int> degree_vector(N);
        for(int i=0; i<N; ++i)
        {
            double r=urd(gen);
            double xx=r*N+1;
            int nn=floor(xx);
            double dd=xx-nn;
            int k;
            if(dd<q_vector[nn-1]){k=nn;}
            else{k=y_vector[nn-1][0];}
            degree_vector[i]=k;
        }
            //after fixing node degree
            for(int i=0; i<N; ++i)
            {
                for(int dd=0; dd<degree_vector[i]; ++dd)
                {
                    node_pocket.emplace_back(i);
                    ++sum_degree;
                }
            }
            uniform_int_distribution<int> random_node(0,sum_degree-1);
            vector<bool> state(N,0);
            unsigned int sum_n=N/2;
            set<int> active_mem;
            while(active_mem.size()<N/2)
            {
                int mem=random_initial(gen);
                if(active_mem.emplace(mem).second){state[mem]=1;}
            }
            double t_save=0.0;
            for(double t=dt; t<=t_max; t+=dt)
            {
                unsigned int voter=random_initial(gen);
                bool state_of_voter=state[voter];
                set<int> hyperedge;
                vector<int> hyperedge_vec;
                hyperedge.emplace(voter);
                int count=0;
                while(hyperedge.size()<s)
                {
                    unsigned int nei=node_pocket[random_node(gen)];
                    if(hyperedge.emplace(nei).second)
                    {
                        hyperedge_vec.emplace_back(nei);
                    }
                }
                for(int pick_count=0; pick_count<q; ++pick_count)
                {
                    int amem=hyperedge_vec[random_member(gen)];
                    if(state[amem]!=state_of_voter){++count;}
                    if(state[amem]==state_of_voter){break;}
                }
                if (count == q)
                {
                    state[voter] = 1 - state[voter];
                    int del_state = 2 * state[voter] - 1;
                    sum_n += del_state;
                }
                //update step
                if(sum_n==N or sum_n==0)
                {
                    t_save=t;
                    break;
                }
            }
            if(sum_n==N)
            {
                exit_times.emplace_back(t_save);
                exit_time+=t_save;
                ++exit_count_p;
                ++exit_count;
                aver_et=(double)exit_time/exit_count;
                aver_ep=(double)exit_count_p/t_ensemble;
            }
            else if(sum_n==0)
            {
                exit_times.emplace_back(t_save);
                exit_time+=t_save;
                ++exit_count;
                aver_et=(double)exit_time/exit_count;
                aver_ep=(double)exit_count_p/t_ensemble;
            }
        if(t_ensemble%10==0)
        {
            sprintf(fname,"FigS4_gam%.1f_q%d_N%.1f_en%d_revision.txt",gam,q,N,t_ensemble);
            fp=fopen(fname,"w");
            double sum_Ti=accumulate(exit_times.begin(),exit_times.end(),0.0);
            double sq_sum=inner_product(exit_times.begin(),exit_times.end(),exit_times.begin(),0.0);
            double stdev=sqrt((sq_sum-2*sum_Ti*aver_et+t_ensemble*aver_et*aver_et)/(t_ensemble-1));
            fprintf(fp,"%.1f\t%.9f\t%.9f\t%.9f\n",N,aver_ep,aver_et,stdev);
            fclose(fp);
        }
        ++t_ensemble;
    }

}



