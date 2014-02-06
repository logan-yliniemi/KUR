/* 
 * File:   Procedural_Transformation.h
 * Author: ylinieml
 *
 * Created on October 22, 2013, 11:58 AM
 */


#ifndef PROCEDURAL_TRANSFORMATION_H
#define	PROCEDURAL_TRANSFORMATION_H

#define PFRONT_THRESHOLD 250
#define PFRONT_BUFFER 2
#define K_FOR_KNN 20

#define OBJECTIVES 2
#define PI 3.1415
#define TIME_DEX 0
#define TREASURE_DEX 1

#define DISCOVERY_THRESHOLD 200

#ifndef VECTOR_INCLUDE
#define VECTOR_INCLUDE
#include <vector.h>
#endif

double vector_dist(vector<double> a, vector<double> b){
    int A=a.size();
    int B=b.size();
    if (A!=B){
        cout << endl;
        cout << "vector_dist sizes don't match" << endl;
    }
    double dx,dxsq,sum=0,d;
    for(int i=0; i<A; i++){
        dx=a.at(i)-b.at(i);
        dxsq=dx*dx;
        sum+=dxsq;
    }
    d=sqrt(sum);
}

class Procedural_Transformation{
    vector<double> utopia;
    vector<double> nadir;
    vector< vector<double> > single_obj_max;
    vector< vector<double> > PFront;
    vector< vector<double> > Domination_Corners;
    vector< vector<double> > scPFront;
    vector< vector<double> > scCorners;

    vector<double> input;
    void take_input(vector<double>* coords);
    
    void scale();
    void Pro_transform();
    void N_Pro_transform(int indicator);
    void N_Dummy_transform();
    
    void calculate_corners();
    void calculate_scaled_corners();
    void calculate_scaled_pareto();
    
    vector<double> output;
    void give_output(vector<double>* coords);

    
public:
    void nad_ut();
    void train();
    
    bool Pareto_Check(vector<double>);
    bool Discrete_Pareto_Check(vector<double>);
    bool Dominated_Check(vector<double> coords);
    int is_Pareto(vector<double> coords);
    
    void Pareto_Reset();
    
    void calc_utopia();
    void calc_nadir();
    void take_objective_max(vector<double> coords);
    void execute_transform(vector<double>* pinputs);
    void execute_N_transform(vector<double>* pinputs,int indicator);
    int get_pareto_size();
    void cout_pareto();
    void print_pareto(FILE*);
    void cout_scaled_pareto();
    vector<double> get_ith_pareto_approximate_member(int i);
};

vector<double> Procedural_Transformation::get_ith_pareto_approximate_member(int i){
    return PFront.at(i);
}

int Procedural_Transformation::get_pareto_size(){
    return PFront.size();
}

void Procedural_Transformation::Pareto_Reset(){
    PFront.clear();
    Domination_Corners.clear();
    scPFront.clear();
    scCorners.clear();
    single_obj_max.clear();
    utopia.clear();
    nadir.clear();
    input.clear();
    output.clear();
}

void Procedural_Transformation::cout_pareto(){
    cout << "Current Non-Dominated Set:" << endl;
    int P = PFront.size();
    for(int p=0; p<P; p++){
        for(int q=0; q<PFront.at(p).size(); q++){
            cout << PFront.at(p).at(q) << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

void Procedural_Transformation::print_pareto(FILE* pFILE){
    cout << "Printing pareto to file" << endl;
    pFILE=fopen("pareto_front.txt","a");
    
    /// Identifier line for split between statistical Pareto front runs.
    for(int q=0; q<PFront.at(0).size(); q++){
        fprintf(pFILE, "%.5f\t", -1.98700);
    }
    fprintf(pFILE, "\n");
    
    int P = PFront.size();
    for(int p=0; p<P; p++){
        for(int q=0; q<PFront.at(p).size(); q++){
            fprintf(pFILE, "%.5f\t", PFront.at(p).at(q));
        }
        fprintf(pFILE, "\n");
    }
    fclose(pFILE);
    cout << "done printing pareto to file" << endl;        
}

void Procedural_Transformation::cout_scaled_pareto(){
    cout << "Current Non-Dominated Set (NORM):" << endl;
    int P = scPFront.size();
    for(int p=0; p<P; p++){
        for(int q=0; q<PFront.at(p).size(); q++){
            cout << scPFront.at(p).at(q) << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

void Procedural_Transformation::calculate_corners(){
    /// Calculate corners of dominated hyperspace. 
    /// These include Pareto Optimal points and joints where Pareto Optimal point shadows intersect.
    
    /// First, Pareto-Optimal Points:
    for(int i=0; i<PFront.size(); i++){
        Domination_Corners.push_back(PFront.at(i));
    }
    /// Then the shadow corners:
    for(int i=0; i<PFront.size()-1; i++){
        vector<double> V;
        V.push_back(PFront.at(i).at(0));
        V.push_back(PFront.at(i+1).at(1));
        Domination_Corners.push_back(V);
    }
}

bool Procedural_Transformation::Dominated_Check(vector<double> coords){
    /// <Is it dominated by any point in the Pareto front??>
    /// For each Pareto point
    for(int pt=0; pt<PFront.size(); pt++){
        int counter=0;
    for(int obj=0; obj<OBJECTIVES; obj++){
        /// If the pareto point scores higher on a criteria, increment counter
        if(coords.at(obj) <= PFront.at(pt).at(obj)){
            counter++;
        }
    }
        if(counter==(OBJECTIVES)){
            /// If the pareto point scored higher or equal on all criteria... 
            /// coords is dominated, and is not a Pareto point.
            /// It is dominated. Return True.
            return true;
        }
    }
    /// It is not dominated. Return False.
    return false;
}

int Procedural_Transformation::is_Pareto(vector<double> coords){
    for(int i=0; i<PFront.size(); i++){
        int yes=1;
        for(int j=0; j<PFront.at(i).size(); j++){
            if(PFront.at(i).at(j)!=coords.at(j)){yes=0;}
        }
        if(yes){return i;}
    }
    return -1;
}

bool Procedural_Transformation::Discrete_Pareto_Check(vector<double> coords){
    /// if we've seen this point X amount of times before, then we pareto check
    static vector< vector<double> > purgatory;
    int not_yet_in_purgatory=1;
    for(int i=0; i<purgatory.size(); i++){
        int inc=1;
        for(int obj=0; obj<OBJECTIVES; obj++){
        if(purgatory.at(i).at(obj) != coords.at(obj)){inc=0; break;}
        }
        if(inc==1){
            not_yet_in_purgatory=0;
            purgatory.at(i).at(OBJECTIVES)+=1.0;
            if(purgatory.at(i).at(OBJECTIVES)>DISCOVERY_THRESHOLD){
                Pareto_Check(coords);
            }
        }
    }
    if(not_yet_in_purgatory==1){
        purgatory.push_back(coords);
        purgatory.back().push_back(0.0);
    }
    
    return false;
    return true;
}

bool Procedural_Transformation::Pareto_Check(vector<double> coords){
    /// <Display the point in question>
    //cout << "Pareto Checking Point: ";
    //for(int i=0; i<coords.size(); i++){
    //    cout << coords.at(i) << "\t";
    //}
    //cout << endl;
    
    /// <Is it dominated by any point in the Pareto front??>
    /// For each Pareto point
    for(int pt=0; pt<PFront.size(); pt++){
        int counter=0;
    for(int obj=0; obj<OBJECTIVES; obj++){
        /// If the pareto point scores higher on a criteria, increment counter
        if(coords.at(obj) <= PFront.at(pt).at(obj)){
            counter++;
        }
    }
        if(counter==(OBJECTIVES)){
            /// If the pareto point scored higher or equal on all criteria... 
            /// coords is dominated, and is not a Pareto point.
            return false;
        }
    }
    
    /// <Does it dominate any points on the Pareto front?>
    vector<int> eliminate;
    for(int pt=0; pt<PFront.size(); pt++){
        int counter=0;
    for(int obj=0; obj<OBJECTIVES; obj++){
        /// If the new point scores higher than or equal to on a criteria, increment counter
        if(coords.at(obj) >= PFront.at(pt).at(obj)){
            counter++;
        }
    }
        if(counter==(OBJECTIVES)){
            /// If the new point scored higher or equal on all criteria
            /// The "Pareto" point is dominated, and should be eliminated.
            eliminate.push_back(pt); 
        }
    }
    
    /// <Eliminate dominated points on the Pareto Front>
    //for(int e=0; e<eliminate.size(); e++){
    //    int spot = eliminate.at(e);
    //    PFront.erase(PFront.begin()+spot);
    //}
    
    /// <Eliminate dominated points on the Pareto Front>
    for(int e=eliminate.size()-1; e>=0; e--){
        /// We eliminate from end -> beginning so that the indices we calculated remain valid.
        int spot = eliminate.at(e);
        PFront.erase(PFront.begin()+spot);
    }
    
    /// <Add new point in correct spot of Pareto Front>
    PFront.push_back(coords);
    /// TODO PFront.emplace()
//    int X=0;
//    /// We want to insert it in the sorted position, so we must find where to insert it.
//    for(int pt=0; pt<PFront.size(); pt++){
//        if(PFront.at(pt).at(0) < coords.at(0)){
//            X++;
//        }
//    }
//    /// X marks the spot. 
//    
//    //Push back dummy vector to make space
//    vector<double> dummy;
//    dummy.push_back(0.1); dummy.push_back(0.2); dummy.push_back(0.3);
//    PFront.push_back(dummy);
//    
//    /// Move all PFront entries above X to spot X+1.
//    for(int pt=PFront.size()-2; pt>X-1; pt--){
//        for(int obj=0; obj<OBJECTIVES; obj++){
//        PFront.at(pt+1).at(obj)=PFront.at(pt).at(obj);
//        }
//    }
//    
//    /// Insert new point into PFront.
//    for(int obj=0; obj<OBJECTIVES; obj++){
//        PFront.at(X).at(obj) = coords.at(obj);
//    }
    
    /// If PFront is over the threshold
    if(PFront.size()>=PFRONT_THRESHOLD+PFRONT_BUFFER){
        /// We find the KNN distance of each.
        vector<double> PFront_KNN_distance;
        int P = PFront.size();
        PFront_KNN_distance.resize(PFront.size());
        for(int piq=0; piq<P; piq++){
            /// find distances from piq to every PFront point.
            vector<double> all_distances;
            for(int i=0; i<P; i++){
            all_distances.push_back(vector_dist(PFront.at(piq),PFront.at(i)));
            }
            /// eliminate K smallest distances
            int K=K_FOR_KNN;
            for(int i=0; i<K; i++){
            int redux = min_element(all_distances.begin(),all_distances.end()) - all_distances.begin();
            all_distances.erase (all_distances.begin()+redux);
            }
            /// remaining minimum is KNN dist.
            int spot = min_element(all_distances.begin(),all_distances.end()) - all_distances.begin();
            PFront_KNN_distance.at(piq)=all_distances.at(spot);
        }
        
        /// PFront_KNN_distance.at(i) now is the ith point on PFront's KNN distance.
        /// Let's eliminate PFRONT_BUFFER amount of PFront and of PFront_KNN_distance.
    /// We eliminate the ones with the lowest KNN distance.
        for(int i=0; i<PFRONT_BUFFER; i++){
            int redux = min_element(PFront_KNN_distance.begin(),PFront_KNN_distance.end()) - PFront_KNN_distance.begin();
            PFront_KNN_distance.erase(PFront_KNN_distance.begin()+redux);
            PFront.erase(PFront.begin()+redux);
        }
        
        /// We ensure that our anchor points still exist.
        bool match=true; 
        bool add=true;
        for(int i=0; i<Anchors.size(); i++){
            add=true;
            for(int p=0; p<PFront.size(); p++){
                match=true;
                for(int obj=0; obj<OBJECTIVES; obj++){
                    if(PFront.at(p).at(obj) != Anchors.at(i).at(obj)){
                        match=false;
                        break;
                    }
                }
                if(match){
                    add=false;
                    break;
                /// do not add 
                }
            }
            if(add){
                PFront.push_back(Anchors.at(i));
            }
        }
        /// And so it is.
    }
    
    /// Since Pareto Front has changed, we recalculate the dominated hyperspace.
    //LYLY //calculate_corners();
    nad_ut();
    calculate_scaled_pareto();
    //calculate_scaled_corners();
    return true;
}

void Procedural_Transformation::calculate_scaled_pareto(){
    /// update vector< vector<double> > scPFront;
    //cout << "Calculating Scaled Pareto\t";
    static int counter;
    if(counter%100==0){
    cout << PFront.size() << "\t" ;
    }
    counter++;
    if(PFront.size() <= 1){return;}
    scPFront.clear();
    for(int i=0; i < PFront.size(); i++){
        vector<double> dual;
        for(int j=0; j < PFront.at(i).size(); j++){
            double val = PFront.at(i).at(j);
            double min = nadir.at(j);
            double range = utopia.at(j)-nadir.at(j);
            dual.push_back((val - min) / range);
        }
        scPFront.push_back(dual);
    }
    
}

void Procedural_Transformation::calculate_scaled_corners(){
    /// update vector< vector<double> > scCorners;
    if(PFront.size() <= 1){return;}
    scCorners.clear();
    
    for(int i=0; i < Domination_Corners.size(); i++){
        vector<double> dual;
        for(int j=0; j<Domination_Corners.at(i).size(); j++){
            double val = Domination_Corners.at(i).at(j);
            double min = nadir.at(j);
            double range = utopia.at(j)-nadir.at(j);
            dual.push_back((val - min) / range);
        }
        scCorners.push_back(dual);
    }

}

void Procedural_Transformation::take_objective_max(vector<double> coords){
    single_obj_max.push_back(coords);
}

void Procedural_Transformation::calc_utopia(){
    /// Calculate Utopia point based on single_obj_max.
    for(int o=0; o<OBJECTIVES; o++){
        double best=-999999.0;
        for(int i=0; i<single_obj_max.size(); i++){
            double val = single_obj_max.at(i).at(o);
            if(val > best){
                best = val;
            }
        }
        utopia.push_back(best);
    }
}

void Procedural_Transformation::nad_ut(){
    /// 1) calculate nadir from PFront
    /// 2) calculate utopia from PFront
    
    // 1)
    nadir.clear();
    double min;
    double mindex;
    
    for(int o=0; o<OBJECTIVES; o++){
        min=99999999999;
        mindex=-1;
        for(int p=0; p<PFront.size(); p++){
        if(PFront.at(p).at(o) < min){
            min=PFront.at(p).at(o);
            mindex=p;
        }
        }
        nadir.push_back(PFront.at(mindex).at(o)); 
    }
    //nadir.push_back(PFront.at(0).at(0)); /// time
    //nadir.push_back(PFront.back().at(1)); /// treasure
    
    // 2)
    utopia.clear();
    double max;
    double maxdex;
    for(int o=0; o<OBJECTIVES; o++){
        max=-999999999999;
        maxdex=-1;
        for(int p=0; p<PFront.size(); p++){
        if(PFront.at(p).at(o)>max){
            max=PFront.at(p).at(o);
            maxdex=p;
        }
        }
        utopia.push_back(PFront.at(maxdex).at(o)+0.00001);   
    }
    /// Added constant avoids atan2 problems at the edges of the PFront.
    //utopia.push_back(PFront.back().at(0)+0.000001); /// time
    //utopia.push_back(PFront.at(0).at(1)+0.000001); /// treasure
    
    //cout << "Nadir: " << nadir.at(0) << "\t" << nadir.at(1) << endl;
    //cout << "Utopia: " << utopia.at(0) << "\t" << utopia.at(1) << endl;
}

void Procedural_Transformation::N_Dummy_transform() {
    output.clear();
    for(int i=0; i<input.size(); i++){
        output.push_back(input.at(i));
    }
}

void Procedural_Transformation::N_Pro_transform(int indicator) {
    /// <LYLY Experimental, Deactivated>
    static vector<vector<double> > popcorn;
    if(indicator==0){
        popcorn.clear();
        for(int i=0; i<scPFront.size(); i++){
            vector<double> pop;
            for(int o=0; o<OBJECTIVES; o++){
                //pop.push_back(LYRAND*0.05);
                pop.push_back(0.0);
            }
            popcorn.push_back(pop);
        }
    }
    
    /// TODO
    /// N Dimensional
    //cout_scaled_pareto();
    //cout << "PRO Transform input 0: " << input.at(0) << endl;
    //cout << "PRO Transform input 1: " << input.at(1) << endl;
    //cout << "PRO Transform input 2: " << input.at(2) << endl;
    vector<double> deltas;
    int I=input.size();
    for(int i=0; i<I; i++){
        deltas.push_back(input.at(i)-1);
        //cout << "Delta " << i << " = " << deltas.at(i) << endl;
        //if(input.at(i)>=1){cout << "!" << i << " " << input.at(i) << endl;}
    }
    
    /// paper notation december 13, 2013.
    /// <FIND d>
    double d=0;
    for(int i=0; i<I; i++){
        d+=deltas.at(i)*deltas.at(i);
    }
    d=sqrt(d);
    //cout << "d = " << d << endl;
    
    /// <FIND D>
    double D;
    /// find distance along utopia->input vector which is first dominated.
    double lowerbound = 0;
    double upperbound = 10;
    double candidate;
    double margin=upperbound-lowerbound;
    bool dominated;
    vector<double> td;
    vector<double> directional_ratios;
    
    double sumdeltassq=0.0;
    for(int i=0; i<I; i++){
        //directional_ratios.push_back(inputs.at(i)/inputs.at(0));
        sumdeltassq+=deltas.at(i)*deltas.at(i);
        //sumdeltassq+=deltas.at(i);
    }
    sumdeltassq=sqrt(sumdeltassq);
    for(int i=0; i<I; i++){
        /// More truly directional cosines.
    directional_ratios.push_back(deltas.at(i)/sumdeltassq);
    //cout << "DIRECTIONAL RATIO " << i << " :: " << directional_ratios.back()<< endl;
    }
    td.resize(I);
    
    /// TODO exception handing for extending the PFront on the extremes.
    /// TODO calculate max upper bound based on number of objectives.
    
    while(margin>0.0001){
        dominated=false;
        //cout << "DOMCHECK: " << dominated << endl;
        candidate=(upperbound+lowerbound)/2;
        //cout << "candidate: " << candidate << endl;
        for(int i=0; i<I; i++){
        td.at(i) = 1 + candidate*directional_ratios.at(i);
        //cout << "point: " << i << " = " << td.at(i) << endl;
        }
        
        for(int p=0; p<scPFront.size(); p++){
            /// for each point in the scaled pareto front
            int counter=0;
            /// if all objectives score better than the td vector
            for(int i=0; i<I; i++){
                //cout << "sc" << scPFront.at(p).at(i) << " " << scPFront.size() << " " << scPFront.at(p).size() << endl;
                //cout << "pop" << popcorn.at(p).at(i) << " " << popcorn.size() << " " << popcorn.at(p).size() << endl;
                //cout << "td" << td.at(i) << " " << td.size() <<  endl;
                
                if(scPFront.at(p).at(i)+popcorn.at(p).at(i)>=td.at(i)){
                    counter++;
                }
                //cout << "COUNTER: " << counter << endl;
            }
            if(counter==I){dominated=true;}
        }
        //cout << "DOMCHECK2: " << dominated << endl;
        if(dominated==true){upperbound = candidate;}
        if(dominated==false){lowerbound = candidate;}
        margin=upperbound-lowerbound;
    }
    D=(upperbound+lowerbound)/2;
    //cout << "D = " << D << endl;
    
    /// <FIND D'>
    double Dprime;
    
    lowerbound = 0;
    upperbound = 10;
    margin=upperbound-lowerbound;
    while(margin>0.0001){
        dominated=false;
        candidate=(upperbound+lowerbound)/2;
        for(int i=0; i<I; i++){
        td.at(i) = 1 + candidate*directional_ratios.at(i);
        //cout << "point: " << i << " = " << td.at(i) << endl;
        }
        //cout << "LC CHECKERS " << (double)accumulate(td.begin(),td.end(),0.0) << endl;
        if(accumulate(td.begin(),td.end(),0.0) <= 1){dominated = true;}
        if(dominated==true){upperbound = candidate;}
        if(dominated==false){lowerbound = candidate;}
        margin=upperbound-lowerbound;
    }
    Dprime = (upperbound+lowerbound)/2;
    //cout << "DPrime = " << Dprime << endl;
    
    /// <CALCULATE d'>
    double dprime = d*Dprime/D;
    //cout << "dprime = " << dprime << endl;
    
    output.clear();
    double checksum=0;
    int stopper;
    for(int i=0; i<I; i++){
        output.push_back(1+dprime*directional_ratios.at(i));
    }
    
    
    
}

void Procedural_Transformation::Pro_transform() { 
    /// ONLY VALID FOR TWO DIMENSIONS
    // Theoretically valid for "n" dimensions with changes (N_Pro_transform())
    // cout << "PRO Transform input 0: " << input.at(0) << endl;
    // cout << "PRO Transform input 1: " << input.at(1) << endl;
    //cout << "input size: " << input.size() << endl;
    double dx = input.at(0)-1;
    double dy = input.at(1)-1;
    // cout << "dy,dx: " << dy << " , " << dx << endl;
    double scaled_theta = atan2(dy,dx);
    double rad = sqrt(dx*dx+dy*dy);

    // cout << "rad: " << rad << "\ttheta: " << scaled_theta << endl;
    
    double newrad;
    double coef;
    coef = 1;
    /// determine coef based on PFront.
    double lowerbound = 0;
    double upperbound = 2;
    double candidate;
    double tx;
    double ty;
    double delta = 2;
    bool dominated;
    while(delta>0.0001){
        dominated=false;
        candidate=(upperbound+lowerbound)/2;
        tx = 1 + candidate*cos(scaled_theta);
        ty = 1 + candidate*sin(scaled_theta);
        for(int i=0; i<PFront.size(); i++){
           // cout << "tx: " << tx << "\tty: " << ty << endl;
           // cout << "scPFront at " << i << " : " << scPFront.at(i).at(0) << "\t" << scPFront.at(i).at(1) << endl;
           if(scPFront.at(i).at(0) >= tx && scPFront.at(i).at(1)>=ty){
               dominated=true;
               break;
           }
        }
        if(dominated==true){upperbound = candidate;}
        if(dominated==false){lowerbound = candidate;}
        delta=upperbound-lowerbound;
    }
    double dominated_radius = (upperbound+lowerbound)/2;
    
    lowerbound = 0;
    upperbound = 2;
    delta=2;
    while(delta>0.0001){
        dominated=false;
        candidate=(upperbound+lowerbound)/2;
        tx = 1 + candidate*cos(scaled_theta);
        ty = 1 + candidate*sin(scaled_theta);
        if(tx+ty <= 1){dominated = true;}
        if(dominated==true){upperbound = candidate;}
        if(dominated==false){lowerbound = candidate;}
        delta=upperbound-lowerbound;
    }
    double ideal_radius = (upperbound+lowerbound)/2;
    
    coef = dominated_radius/ideal_radius;
    
    // cout << "dominated_radius: " << dominated_radius << endl;
    // cout << "ideal radius: " << ideal_radius << endl;
    // cout << "Coef: " << coef << endl;
    
    newrad = rad/coef;
    
    // cout << "newrad: " << newrad << endl;
    
    /// We need to determine the x,y values to give the transform output.
    double x = 1 + newrad * cos(scaled_theta);
    double y = 1 + newrad * sin(scaled_theta);
    
    // cout << "Newx: " << x << "\tNewy: " << y << endl;
    
    output.clear();
    //for (int a = 0; a < 1; a++) {
        output.push_back(x);
        output.push_back(y);
    //}
}


void Procedural_Transformation::calc_nadir(){
    for(int o=0; o<OBJECTIVES; o++){
        double worst=999999.0;
        for(int i=0; i<single_obj_max.size(); i++){
            double val = single_obj_max.at(i).at(o);
            if(val < worst){
                worst = val;
            }
        }
        nadir.push_back(worst);
    }
}

void Procedural_Transformation::take_input(vector<double>* pcoords){
    /// Assign input to class variable.
    // cout << "inside take_input\t" <<  pcoords->size() << endl;
    if(pcoords->size() != OBJECTIVES){cout << "Are we doing a " << OBJECTIVES << " objective problem or not?" << endl;}
    input.clear();
    for(int obj=0; obj<OBJECTIVES; obj++){
        input.push_back(pcoords->at(obj));
    }
    // cout << "inside take_input2\t" <<  input.size() << endl;
}

void Procedural_Transformation::give_output(vector<double>* pcoords){
    /// Assign output to external variable.
    for(int obj=0; obj<OBJECTIVES; obj++){
        pcoords->at(obj) = output.at(obj);
    }
}

void Procedural_Transformation::scale(){
    /// Scale values of objectives to be less than one, with the nadir point taking on (0,0).
    for(int obj=0; obj<OBJECTIVES; obj++){
        double range = utopia.at(obj)-nadir.at(obj);
        input.at(obj) = (input.at(obj) - nadir.at(obj)) / range;
    }
}

void Procedural_Transformation::execute_N_transform(vector<double>* pinputs, int indicator){
    int I=pinputs->size();
    take_input(pinputs);
    // cout << "INPUTS TAKEN" << endl;
    // for(int i=0; i<I; i++){
    //     cout << i << ":" << pinputs->at(i) << "\t";
    // }
    // cout << endl;
    scale();
    // cout << "INPUTS SCALED" << endl;
    // for(int i=0; i<I; i++){
    //     cout << i << ":" << pinputs->at(i) << "\t";
    // }
    // cout << endl;
    N_Pro_transform(indicator);
    //N_Dummy_transform();
    // cout << "TRANSFORMED!!!" << endl;
    // for(int i=0; i<I; i++){
    //     cout << i << ":" << pinputs->at(i) << "\t";
    // }
    // cout << endl;
    give_output(pinputs);
    // cout << "OUTPUTS:" << endl;
    // for(int i=0; i<I; i++){
    //     cout << i << ":" << pinputs->at(i) << "\t";
    // }
    // cout << endl;
}

void Procedural_Transformation::execute_transform(vector<double>* pinputs){
    /// Public method for performing the whole transformation.
    //cout << "pinputs size: " << pinputs->size() << endl;
    //cout << "pinputs: " << pinputs->at(0) << "\t" << pinputs->at(1) << endl;
    take_input(pinputs);
    //cout << "Inputs: " << input.at(0) << "\t" << input.at(1) << endl;
    scale(); 
    //cout << "Scaled Inputs: " << input.at(0) << "\t" << input.at(1) << endl;
    Pro_transform();
    give_output(pinputs);
    //cout << "Outputs: " << output.at(0) << "\t" << output.at(1) << endl;
}

void Pro_Pareto_testing(){
    Procedural_Transformation T;
    vector<double> coords;
    for(int i=0; i<1000; i++){
    cout << "Trial " << i << endl;
    coords.push_back(rand()%100000);
    coords.push_back(rand()%100000);
    cout << "coords 1\t" << coords.at(0) << endl;
    cout << "coords 2\t" << coords.at(1) << endl;
    T.Pareto_Check(coords);
    coords.clear();
    T.cout_pareto();
    }
    cout << "Placement" << endl;
    coords.push_back(100000);
    coords.push_back(100000);
    T.Pareto_Check(coords);
    T.cout_pareto();
    coords.clear();
}

#endif	/* PROCEDURAL_TRANSFORMATION_H */

