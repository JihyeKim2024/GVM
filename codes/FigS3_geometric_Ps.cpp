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
int N = 10000, n_size = N - 1, t_ensemble = 1;
int q=2;
double t_max = 1000000000.0, dt = 1.0 / N;
double smin=3;
int num_data = 6;
uniform_int_distribution<int> random_node(0,N-1);
uniform_real_distribution<double> urd(0,1);
vector<double> aver_et(num_data,0.0);
vector<double> aver_ep(num_data,0.0);
vector<double> exit_time(num_data,0.0);
vector<double> exit_count(num_data,0.0);
vector<double> exit_count_p(num_data,0.0);
vector<vector<double>> exit_times(num_data);
vector<vector<double>> q_vector(num_data);
vector<vector<double>> updated_qvector(num_data);
vector<vector<int>> y_vector(n_size);
vector<vector<int>> poor_vec(num_data);
vector<vector<int>> rich_vec(num_data);
vector<double> denom(num_data);
FILE*fp;
char fname[100];

int main()
{
    //making a table for Walker algorithm
    for(int i=0; i< num_data; ++i)
    {
        double mm = i + smin; // <s>
        denom[i]=0;
        
        for(int j=N; j>1; --j)
        {
            denom[i]+=j*pow((mm-2)/(mm-1),j-2)/(mm-1);
        } 
        for(int ss=0; ss<N-1; ++ss)
        {
            int cs=ss+2;
            double p_cs=cs*pow((mm-2)/(mm-1),ss)/(mm-1)/denom[i]; //geometric P(s)
            q_vector[i].emplace_back(n_size*p_cs);
            updated_qvector[i].emplace_back(n_size*p_cs);
            if(q_vector[i][ss]>1)
                rich_vec[i].emplace_back(ss);
            else if(q_vector[i][ss]<1)
                poor_vec[i].emplace_back(ss);
        }
        while(poor_vec[i].size()!=0)
        {
            int poorer=poor_vec[i][0];
            int richer=rich_vec[i][0];
            y_vector[poorer].emplace_back(richer+2);
            updated_qvector[i][richer]-=(1-q_vector[i][poorer]);
            poor_vec[i].erase(poor_vec[i].begin());
            if(updated_qvector[i][richer]<1)
            {
                rich_vec[i].erase(rich_vec[i].begin());
                poor_vec[i].emplace_back(richer);
            }
        }
        for(int ss=0; ss<N-1; ++ss)
        {
            if(y_vector[ss].size()!=i+1)
                y_vector[ss].emplace_back(-1);
        }
    } //making a table for Walker algorithm

    for(int hyper_en=0; hyper_en<100000000; ++hyper_en)
    {
        for(int m=0; m< num_data; ++m)
        {
            vector<bool> state(N);
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
                double r=urd(gen);
                double x=r*n_size+1;
                int nn=floor(x);
                double dd=x-nn;
                double s;
                if(dd<q_vector[m][nn-1])
                    s=nn+1;
                else
                    s=y_vector[nn-1][m];
                int count=0;
                while(hyperedge.size()<s)
                {
                    unsigned int nei=random_node(gen);
                    if(hyperedge.emplace(nei).second)
                    {
                        hyperedge_vec.emplace_back(nei);
                    }
                }
                uniform_int_distribution<int> random_member(0,s-2);
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
                }//update step
                if(sum_n==N or sum_n==0)
                {
                    t_save=t;
                    break;
                }
            }//t
            if(sum_n==N)
            {
                exit_times[m].emplace_back(t_save);
                exit_time[m]+=t_save;
                ++exit_count_p[m];
                ++exit_count[m];
                aver_et[m]=(double)exit_time[m]/exit_count[m];
                aver_ep[m]=(double)exit_count_p[m]/t_ensemble;
            }
            else if(sum_n==0)
            {
                exit_times[m].emplace_back(t_save);
                exit_time[m]+=t_save;
                ++exit_count[m];
                aver_et[m]=(double)exit_time[m]/exit_count[m];
                aver_ep[m]=(double)exit_count_p[m]/t_ensemble;
            }
            
        }// <s>
        if(t_ensemble==10 or t_ensemble==100 or t_ensemble==1000 or t_ensemble==5000 or t_ensemble==10000 or t_ensemble==50000 or t_ensemble==100000 or t_ensemble==1000000 or t_ensemble==10000000 or t_ensemble==100000000 or t_ensemble==1000000000)
        {
            sprintf(fname,"FigS3_q%d_geo_en%d.txt",q,t_ensemble);
            fp=fopen(fname,"w");
            for(int i=0; i<5; ++i)
            {
                double sum_Ti=accumulate(exit_times[i].begin(),exit_times[i].end(),0.0);
                double sq_sum=inner_product(exit_times[i].begin(),exit_times[i].end(),exit_times[i].begin(),0.0);
                double stdev=sqrt((sq_sum-2*sum_Ti*aver_et[i]+t_ensemble*aver_et[i]*aver_et[i])/(t_ensemble-1));
                fprintf(fp,"%.9f\t%.9f\t%.9f\t%.9f\n",  i + smin,aver_ep[i],aver_et[i],stdev);
            }
            fclose(fp);
        }
        ++t_ensemble;
    }// ensemble 
}
