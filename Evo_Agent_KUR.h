/* 
 * File:   Evo_Agent_KUR_KUR.h
 * Author: ylinieml
 *
 * Created on December 7, 2013, 12:50 PM
 */

#ifndef Evo_Agent_KUR_KUR_H
#define	Evo_Agent_KUR_KUR_H

/* 
 * File:   Evo_Agent_KUR_DST.h
 * Author: ylinieml
 *
 * Created on October 15, 2013, 12:27 PM
 */

#ifndef Evo_Agent_KUR_H
#define	Evo_Agent_KUR_H

#ifndef VECTOR_INCLUDE
#define VECTOR_INCLUDE
#include <vector.h>
#endif

#define PI 3.1415

class Evo_Agent_KUR {
    double raw_fitness; /// raw fitness calculated by treasure, time
    double f1;
    double f2;
    
    double transformed_treasure;
    double transformed_time;

    void console_1vector(vector<double>);
    void console_2vector(vector<vector<double> >);

    double Initial_Q_Value;
    
    
    
    void create_action_vector();

public:
    vector<double> actions;
    double transformed_fitness; /// transformed fitness after Pareto Transformation
    double fitness; /// used for actual evolutionary method
    double get_action(int time);
    double get_fitness();
    void mutate();
    double action;
    void show_action_vector();
    int id;
    
    void start();       //Used before first run of a statistical run / repeat
    void reset();       //Used at start of episode
    
    void set_f1(double val);
    void set_f2(double val);
    double get_f1();
    double get_f2();
    
    void set_raw_fitness(double raw);
    void set_fitness(double fit);
};

void Evo_Agent_KUR::start() {
    /// id is set outside of start, in main();
    create_action_vector();
    fitness=-1;
}

void Evo_Agent_KUR::reset() {
    f2=0;
    f1=0;
}

void Evo_Agent_KUR::create_action_vector(){
    cout << "action vector in! " << id << endl;
    double r1;
    actions.clear();
    for(int i=0; i<3; i++){
        r1=(double)rand()/RAND_MAX;
        r1=-5+r1*10;
        actions.push_back(r1);
    }
}

double Evo_Agent_KUR::get_fitness(){
    return fitness;
}

double Evo_Agent_KUR::get_action(int spot){
    return actions.at(spot);
}

void Evo_Agent_KUR::set_f1(double val){
    f1=val;
}

double Evo_Agent_KUR::get_f1(){
    return f1;
}

void Evo_Agent_KUR::set_f2(double val){
    f2=val;
}

double Evo_Agent_KUR::get_f2(){
    return f2;
}

void Evo_Agent_KUR::set_raw_fitness(double raw){
    raw_fitness=raw;
}

void Evo_Agent_KUR::set_fitness(double fit){
    fitness=fit;
}

void Evo_Agent_KUR::console_2vector(vector< vector<double> > a) {
    for (int i = 0; i < a.size(); i++) {
        console_1vector(a.at(i));
        cout << endl;
    }
}

void Evo_Agent_KUR::console_1vector(vector<double> a) {
    for (int i = 0; i < a.size(); i++) {
        cout << a.at(i);
        cout << "\t";
    }
}

void Evo_Agent_KUR::show_action_vector() {
    int A=actions.size();
    for(int a=0; a<A; a++){
        cout << actions.at(a) << "\t";
    }
    cout << endl;
}

void Evo_Agent_KUR::mutate() {
    /// <PARAM>
    actions.resize(3);
    actions.at(0)+=LYrand_norm(0.3);
    actions.at(1)+=LYrand_norm(0.3);
    actions.at(2)+=LYrand_norm(0.3);
    for(int i=0; i<3; i++){
    if(actions.at(i) < -5){actions.at(i) = -5;}
    if(actions.at(i) > 5 ){actions.at(i) = 5;}
    }
}

#endif	/* Evo_Agent_KUR_H */




#endif	/* Evo_Agent_KUR_KUR_H */

