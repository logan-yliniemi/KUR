/* 
 * File:   KUR.cpp
 * Author: ylinieml
 */

/// Evoluationary Algorithm Parameters
/// <PARAM>
#define POPULATION 100
#define ELIMINATE 50
#define GENERATIONS 2500
#define STEPS 1
#define STAT_RUNS 1
#define BETA 0.5


#define DO_LC 1
#define DO_PACCET 0
#define DO_NSGA 0
#define DO_SPEA 0

#define USE_ANCHORS 0

#include <cstdlib>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <numeric>

#include <stdio.h>
#include <iostream>
using namespace std;

#ifndef VECTOR_INCLUDE
#define VECTOR_INCLUDE
#include <vector>
#endif

#define PI 3.1415

/// Small Functions/Macros
#define LYRAND (double)rand()/RAND_MAX
double LYrand_norm(double a){
    double theta=LYRAND*2*PI;
    double rsq=-1*a*log(LYRAND);
    double x=rsq*cos(theta);   
    return(x);
}
#define SMALL 0.0001
#define HUNDI_ROUND(x) (double)floor(x*100)/100

vector< vector<double> > Anchors;

bool pretty_print = true;

using namespace std;

double vector_median(vector<double> fit){
    double median;
    /// sort vector
    sort(fit.begin(),fit.end());
    /// even or odd size
    int even = (fit.size()+1)%2;
    int odd = fit.size()%2;
    if(even){ /// even case
        int mid = fit.size()/2;
        double m1 = fit.at(mid);
        double m2 = fit.at(mid+1);
        median = (m1+m2)/2;
    }
    if(odd){ /// odd case
        int mid = fit.size()/2;
        median = fit.at(mid);
    }
    return median;
}
double vector_mean(vector<double> fit){
    double sum=accumulate(fit.begin(), fit.end(), 0.0);
    double mean = sum / fit.size();
    return mean;
}

void report(FILE* pFILE, double value) { /// report to text file
    fprintf(pFILE, "%.5f\t", value);
}

void newline(FILE* pFILE) { /// report to text file
    fprintf(pFILE, "\b \b\n");
}


#include "Evo_Agent_KUR.h"
#include "Procedural_Transformation.h"
#include "NSGAheader.h"
#include "SPEAheader.h"

void grid_visualize(Procedural_Transformation* pT);
void contour_visualize(Procedural_Transformation* pT);

class KURclass{
public:
    double f1(double x1,double x2,double x3);
    double f2(double x1,double x2,double x3);
    void start();
};

void KURclass::start(){
    
}

double KURclass::f1(double x1, double x2, double x3){
    double val=0;
    val += -10*exp(-0.2 * sqrt(pow(x1,2)+pow(x2,2)));
    val += -10*exp(-0.2 * sqrt(pow(x2,2)+pow(x3,2)));
    return val;
}
double KURclass::f2(double x1, double x2, double x3){
    double val=0;
    val+=pow(fabs(x1),0.8) + 5*sin(pow(x1,3));
    val+=pow(fabs(x2),0.8) + 5*sin(pow(x2,3));
    val+=pow(fabs(x3),0.8) + 5*sin(pow(x3,3));
    return val;
}

int main(){
    srand(time(NULL));
    
    //srand(14);
    //Pro_Pareto_Filter_Testing();
    //Procedural_Testing();
    //int iii;
    //cin >> iii;
    
    FILE* pFILE_fit;
    FILE* pFILE_time;
    FILE* pFILE_treasure;
    FILE* pFILE_pareto_number;
    FILE* pFILE_pareto_discovery;
    FILE* pFILE_pareto_front;
    pFILE_fit=fopen("fitness.txt","w");
    pFILE_time=fopen("time.txt","w");
    pFILE_treasure=fopen("treasure.txt","w");
    pFILE_pareto_number=fopen("pareto_size.txt","w");
    pFILE_pareto_discovery=fopen("pareto_discovery.txt","w");
    
    Procedural_Transformation T;
    SPEA_2 SPEA;
    
    if(USE_ANCHORS){
    vector<double> one;
    vector<double> two; 
    
    one.push_back(20); 
    one.push_back(0); 
    
    two.push_back(14.44);
    two.push_back(11.61);
    
    Anchors.push_back(one);
    Anchors.push_back(two);
    }
    
    for(int stat_run=0; stat_run < STAT_RUNS; stat_run++) {
            T.Pareto_Reset();
            vector<double> badness;
            badness.push_back(1000.0);
            badness.push_back(1000.0);
            T.Pareto_Check(badness);
            
            KURclass environment;
            environment.start();
            KURclass* pE = &environment; /// pointer to Environment

            vector<Evo_Agent_KUR> Agents;
            vector<Evo_Agent_KUR>* pVA = &Agents; /// pointer to Vector of Agents

            for (int i = 0; i < POPULATION; i++) {
                Evo_Agent_KUR EA;
                EA.id=i;
                EA.start();
                pVA->push_back(EA);
            }

            for (int gen = 0; gen < GENERATIONS; gen++){
                if (gen % (GENERATIONS / 10) == 0) {
                    cout << "Run No." << stat_run << " is " << (double) gen / GENERATIONS * 100 << " % Complete!" << endl;
                }
                
                /// For each population member in pA, execute 1 round of the DST domain:
                for (int mem=0; mem<POPULATION; mem++) {
                    Evo_Agent_KUR* pA = &pVA->at(mem);
                    pA->reset();
                    /// KUR is a bi-minimization problem. This implementation of PaCcET is for maximization.
                    /// Objective reversal is done here and here only.
                    double f1hold = pE->f1(pA->get_action(0),pA->get_action(1),pA->get_action(2));
                    double f2hold = pE->f2(pA->get_action(0),pA->get_action(1),pA->get_action(2));
                    pA->set_f1(f1hold);
                    pA->set_f2(f2hold);
               }
                
                /// we set up the pareto front:
                for (int a = 0; a < pVA->size(); a++) {
                    vector<double> MO;
                    MO.push_back(pVA->at(a).get_f1());
                    MO.push_back(pVA->at(a).get_f2());
                    T.Pareto_Check(MO);
                }
                
                /// and now we do comparisons.
                for (int a = 0; a < pVA->size(); a++) {
                    /// <PARAM>
                    /// <Linear Combination of Objectives>
                    double tr=0.5;
                    pVA->at(a).fitness = (pVA->at(a).get_f1())*tr + pVA->at(a).get_f2()*(1-tr);
                    /// <F1 Only>
                    // pVA->at(a).fitness = pVA->at(a).get_f1();
                    /// <F2 Only>
                    //pVA->at(a).fitness = pVA->at(a).get_f2();
                    /// <Procedural Transformation>
                    vector<double> MO;
                    vector<double>* pMO;
                    vector<double> OMO;
                    pMO = &MO;
                    MO.push_back(pVA->at(a).get_f1());
                    MO.push_back(pVA->at(a).get_f2());
                    OMO=MO; /// Copy for pareto check after transforming
                    //cout << "KUR MO1 " << MO.at(0) << endl;
                    //cout << "KUR MO2 " << MO.at(1) << endl;
                    
                    //if(DO_NSGA){
                    //NSGA.vector_input(MO,a);
                    //}
                    if(DO_SPEA){
                    SPEA.vector_input(MO,a);
                    }
                    if(DO_PACCET){
                    T.execute_N_transform(pMO);
                    //T.Pareto_Check(OMO);
                    
                    pVA->at(a).transformed_fitness = MO.at(0)*BETA + MO.at(1)*(1-BETA);
                    pVA->at(a).fitness = pVA->at(a).transformed_fitness;
                    }
                
                    
            }
                
                if(DO_NSGA){
                    NSGA_2 NSGA;
                    NSGA.declare_NSGA_dimension(2);
                    NSGA.NSGA_reset();
                    for (int a = 0; a < pVA->size(); a++) {
                        vector<double> afit;
                        afit.push_back(pVA->at(a).get_f1());
                        afit.push_back(pVA->at(a).get_f2());
                        //afit.push_back(pVA->at(a).get_fxn(2));
                        NSGA.vector_input(afit,a);
                    }
                    NSGA.execute();
                    for (int a = 0; a < pVA->size(); a++) {
                        pVA->at(a).fitness=-NSGA.NSGA_member_fitness(a);
                    }
                }
                
                if(DO_SPEA){
                    for(int a=0; a<pVA->size(); a++){
                    vector<double> MO;
                    MO.push_back(pVA->at(a).get_f1());
                    MO.push_back(pVA->at(a).get_f2());
                    SPEA.vector_input(MO,a);
                    SPEA.take_agent(pVA->at(a),a);
                    }
                    
                vector<int> survivors;
                vector<int>* pS=&survivors;
                SPEA.execute(pS);
                
                pVA->clear();
                for(int i=0; i< pS->size(); i++){
                    int el = pS->at(i);
                    pVA->push_back(SPEA.archive.at(el).agent);
                    pVA->back().mutate();
                }
                }
                
                vector<double> fit;
                vector<double> times;
                vector<double> treasures;
                
                for (int a = 0; a < pVA->size(); a++) {
                    fit.push_back(-pVA->at(a).get_fitness());
                    times.push_back(pVA->at(a).get_f1());
                    treasures.push_back(pVA->at(a).get_f2());
                   report(pFILE_fit,fit.back()); /// Report every result
                   report(pFILE_time,times.back()); // Report every result
                   report(pFILE_treasure,treasures.back()); // Report every result
                   report(pFILE_pareto_number,T.get_pareto_size());
                }
                /// determine mean, median for reporting
                double generation_median = vector_median(fit);
                double generation_mean = vector_mean(fit);
                double median_time = vector_median(times);
                double median_treasure = vector_median(treasures);
                if(gen % 100 == 0){
                cout << "Generation\t" << gen << "\tf1:\t" << median_time << "\tf2:\t" << median_treasure << "\tfitness:\t" << generation_median << endl;
                }
                
                /// always eliminate worst-performing solutions.
                if(!DO_SPEA){
                for(int e=0; e<ELIMINATE; e++){
                    double minfit = *min_element(fit.begin(), fit.end());
                    int spot = min_element(fit.begin(), fit.end()) - fit.begin();
                    // kill this fitness;
                    fit.erase(fit.begin() + spot);
                    // kill this agent
                    pVA->erase(pVA->begin() + spot);
                }
                
                /// duplicate best-performing agents
                for(int d=0; d<ELIMINATE; d++){
                    double maxfit = *max_element(fit.begin(), fit.end());
                    int spot = max_element(fit.begin(), fit.end()) - fit.begin();
                    if(rand()%10)
                    {
                        int az=rand();
                        int sz=pVA->size();
                        int spt=az%sz;
                        spot=spt;
                    }
                    
                    /// <PARAM>
                    /// to reduce the rate of convergence, we select the best one (spot) 50% of the time,
                    /// and the other 50% of the time we select a random survivor.
                    
              
                    /// create exact copy
                    pVA->push_back(pVA->at(spot));
                    /// mutate
                    pVA->back().mutate();
                }
                }
                
                if (pretty_print) {
                   report(pFILE_fit,generation_mean); /// report every result
                   report(pFILE_time,median_time); // Report every result
                   report(pFILE_treasure,median_treasure); // Report every result
                   report(pFILE_pareto_number,T.get_pareto_size());
               } else {                
                   //For Coarse Results
                   if (gen % (GENERATIONS / 100) == 0) {
                       report(pFILE_fit,generation_mean); // Report only occasionally
                       report(pFILE_time,median_time); // Report only occasionally
                       report(pFILE_treasure,median_treasure); // Report only occasionally
                       report(pFILE_pareto_number,T.get_pareto_size());
                   }
               }
                
                
            
    }
            //Start a new line in output file for next run
            fprintf(pFILE_fit,"\n");
            fprintf(pFILE_time,"\n");
            fprintf(pFILE_treasure,"\n");
            fprintf(pFILE_pareto_number,"\n");
            T.print_pareto(pFILE_pareto_front);
}      
    fclose(pFILE_time);
    fclose(pFILE_treasure);
    fclose(pFILE_fit);
    fclose(pFILE_pareto_number);
    T.cout_pareto();
    
    Procedural_Transformation* pT=&T;
    //grid_visualize(pT);
    //contour_visualize(pT);
}

void grid_visualize(Procedural_Transformation* pT){
    cout << "Grid Visualize Starting" << endl;
    
    FILE* pFILE_coord;
    FILE* pFILE_dom;
    FILE* pFILE_trans;
    
    pFILE_coord=fopen("before_coord.txt","w");
    pFILE_dom=fopen("dom.txt","w");
    pFILE_trans=fopen("after_coord.txt","w");
    
    vector<vector< double > > line;
    vector<double>* pcoord;
    
    /// Decide on grid line spacing
    int grid_lines,spaces;
    int ppl;
    double xmin,xmax,ymin,ymax;
    
    /// <PARAM>
    grid_lines=10;
    spaces=grid_lines-1;
    xmax = 19.9;
    xmin = 14.44;
    ymax = 11.42;
    ymin = -0.1;
    ppl = 200;
    
    double x_const_spacing = (xmax-xmin) / spaces;
    double y_const_spacing = (ymax-ymin) / spaces;
    double x_dot_spacing = (xmax-xmin) / ppl;
    double y_dot_spacing = (ymax-ymin) / ppl;
    
    double xconst;
    double xvar;
    double yconst;
    double yvar;
    
    /// x constant lines
    for(int ln=0; ln<grid_lines; ln++){
        line.clear();
    /// Make Single Grid line
         xconst= xmin + x_const_spacing*ln;
         for(int pt=0; pt<ppl; pt++){
         yvar= ymin + y_dot_spacing * pt;
         vector<double> point;
         point.push_back(xconst);
         point.push_back(yvar);
         line.push_back(point);
         } 
    for(int i=0; i<line.size(); i++){
    /// Print out before
        for(int j=0; j<line.at(i).size(); j++){
        report(pFILE_coord,line.at(i).at(j));
        }
        report(pFILE_dom,pT->Dominated_Check(line.at(i)));
        newline(pFILE_coord);
        newline(pFILE_dom);
    /// Take one coords at a time, feed it through transformation
        pcoord = &line.at(i);
        pT->execute_N_transform(pcoord);
    /// Print out after
        for(int j=0; j<line.at(i).size(); j++){
        report(pFILE_trans,line.at(i).at(j));
        }
        newline(pFILE_trans);
    }
    
    }
    
    /// y constant lines
    for(int ln=0; ln<grid_lines; ln++){
        line.clear();
    /// Make Single Grid line
        yconst = ymin + y_const_spacing*ln;
        for(int pt=0; pt<ppl; pt++){
            xvar= xmin + x_dot_spacing*pt;
            vector<double> point;
            point.push_back(xvar);
            point.push_back(yconst);
            line.push_back(point);
        }
        
    for(int i=0; i<line.size(); i++){
    /// Print out before
        for(int j=0; j<line.at(i).size(); j++){
        report(pFILE_coord,line.at(i).at(j));
        }
        report(pFILE_dom,pT->Dominated_Check(line.at(i)));
        newline(pFILE_coord);
        newline(pFILE_dom);
    /// Take one coords at a time, feed it through transformation
        pcoord = &line.at(i);
        pT->execute_N_transform(pcoord);
    /// Print out after
        for(int j=0; j<line.at(i).size(); j++){
        report(pFILE_trans,line.at(i).at(j));
        }
        newline(pFILE_trans);
    
    }
    }
    fclose(pFILE_coord);
    fclose(pFILE_dom);
    fclose(pFILE_trans);
    
    
    vector<double> held;
    vector<double>* pH;
    
    FILE* pFILE_orig_pareto=fopen("orig_pareto.txt","w");
    FILE* pFILE_trans_pareto=fopen("trans_pareto.txt","w");
    for(int i=0; i<pT->get_pareto_size(); i++){
        /// print out original Pareto front
        held=pT->get_ith_pareto_approximate_member(i);
        pH=&held;
        for(int obj=0; obj<OBJECTIVES; obj++){
        report(pFILE_orig_pareto,held.at(obj));
        }
        /// print out transformed Pareto front
        pT->execute_N_transform(pH);
        for(int obj=0; obj<OBJECTIVES; obj++){
        report(pFILE_trans_pareto,held.at(obj));
        }
        newline(pFILE_orig_pareto);
        newline(pFILE_trans_pareto);
    }
    fclose(pFILE_orig_pareto);
    fclose(pFILE_trans_pareto);
    
    cout << "Grid Visualize Ending" << endl;
}

void contour_visualize(Procedural_Transformation* pT){
    cout << "Contour Visualize Starting" << endl;
    
    int EQUAL_CONTOUR=1;
    int STEEPEST_GROWTH=0;
    
    FILE* pFILE_coord;
    FILE* pFILE_dom;
    FILE* pFILE_trans;
    
    pFILE_coord=fopen("regular_space_coords.txt","w");
    pFILE_dom=fopen("dom.txt","w");
    pFILE_trans=fopen("transformed_space_coords.txt","w");
    
    vector<vector< double > > line;
    vector<double>* pcoord;
    
    /// Decide on grid line spacing
    int grid_lines,spaces;
    int ppl;
    double min_contour,max_contour,width;
    
    /// <PARAM>
    grid_lines=60;
    spaces=grid_lines-1;
    min_contour = 0;
    max_contour = 1;
    width = 1.4;
    ppl = 200;
    
    double line_spacing = (max_contour - min_contour) / spaces;
    double dot_spacing = width/ppl;
    
    
    double contour = min_contour;
    for(int ln=0; ln<grid_lines; ln++){
        line.clear();
    /// Make Single Contour line
    /// Contour Value
        contour += line_spacing;
        double x,y;
        if(EQUAL_CONTOUR){
            x = contour*sqrt(2) - width*sqrt(2)/2; 
            y = contour*sqrt(2) + width*sqrt(2)/2;
        }
        if(STEEPEST_GROWTH){
            x = contour*sqrt(2) - width*sqrt(2)/2; 
            y = - contour*sqrt(2) + width*sqrt(2)/2;
        }
         for(int pt=0; pt<ppl; pt++){
         vector<double> point;
         point.push_back(x);
         point.push_back(y);
         
         
         line.push_back(point);
         if(EQUAL_CONTOUR){
         y-=dot_spacing*sqrt(2);
         x+=dot_spacing*sqrt(2);
         }
         if(STEEPEST_GROWTH){
         y+=dot_spacing*sqrt(2);
         x+=dot_spacing*sqrt(2);
         }
         
         }
        
        /// Print out transformed space
    for(int i=0; i<line.size(); i++){
        for(int j=0; j<line.at(i).size(); j++){
        report(pFILE_trans,line.at(i).at(j));
        }
        newline(pFILE_trans);
        /// ^V.W.A.I. Feb 25, 2014.
        
        /// Reverse Transformation
        /// Take one coords at a time, feed it through transformation
        pcoord = &line.at(i);
        pT->execute_N_reverse_transform(pcoord);

    /// Print out untransformed
        for(int j=0; j<line.at(i).size(); j++){
        report(pFILE_coord,line.at(i).at(j));
        }
        report(pFILE_dom,pT->Dominated_Check(line.at(i)));
        newline(pFILE_coord);
        newline(pFILE_dom);
    }
    
    }
    
    fclose(pFILE_coord);
    fclose(pFILE_dom);
    fclose(pFILE_trans);
    
    
    vector<double> held;
    vector<double>* pH;
    
    FILE* pFILE_orig_pareto=fopen("orig_pareto.txt","w");
    FILE* pFILE_trans_pareto=fopen("trans_pareto.txt","w");
    for(int i=0; i<pT->get_pareto_size(); i++){
        /// print out original Pareto front
        held=pT->get_ith_pareto_approximate_member(i);
        pH=&held;
        for(int obj=0; obj<OBJECTIVES; obj++){
        report(pFILE_orig_pareto,held.at(obj));
        }
        /// print out transformed Pareto front
        pT->execute_N_transform(pH);
        for(int obj=0; obj<OBJECTIVES; obj++){
        report(pFILE_trans_pareto,held.at(obj));
        }
        newline(pFILE_orig_pareto);
        newline(pFILE_trans_pareto);
    }
    fclose(pFILE_orig_pareto);
    fclose(pFILE_trans_pareto);
    
    cout << "Contour Visualize Ending" << endl;
}