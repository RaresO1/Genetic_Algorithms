#include <iostream>
#include <fstream>
#include <Windows.h>

#define _USE_MATH_DEFINES
#include <cmath>
#include <chrono>
#include <random>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <sys/types.h> 
#include <unistd.h>

using namespace std;

ofstream g("C:/Users/Rares/Documents/Facultate/An2/Sem1/AG/H_2p/res.txt", std::ios_base::app);
//bool v[901], v2[901];


vector<vector<bool>> pop, ofs;
vector<vector<double>> values, ofs_values;
vector<double>eval(10000), q_eval(10000), p(10000), ofs_eval(10000);

double best_of_all = 999999999, best_of_gen;
int nr_dim, nr_biti, nr_bits_in_dim, nr_cand, pop_size, max_gen, gen_nr;
double pc, pm;
int norm_val;
chrono::high_resolution_clock::time_point t1,t2;

unsigned t = chrono::high_resolution_clock::now().time_since_epoch().count();

std::uniform_real_distribution<> distribution(0.0, 1.0);



void getValues(vector<vector<bool>> &v, vector<vector<double>> &val, int cand, double lim_dr, double lim_st){
    if(val.size() <= cand){
        val.push_back(vector<double> (nr_dim));
    }
    //decodificare din gray in binar
    for(int i=0; i<nr_dim; ++i){

        val[cand][i]=v[cand][i*nr_bits_in_dim];
        for(int j=0; j<nr_bits_in_dim-1; ++j){
            if(v[cand][i*nr_bits_in_dim + j] == v[cand][i*nr_bits_in_dim + j+1]){
                val[cand][i] = val[cand][i] * 2;
            }
            else {
                val[cand][i] = val[cand][i] * 2 + 1;
            }
        }
        // val[cand][i] = 0;
        // for(int j=0; j<nr_bits_in_dim; ++j){
        //     val[cand][i] = val[cand][i] * 2 + v[cand][i*nr_bits_in_dim + j];
        // }

        val[cand][i] = (val[cand][i] / norm_val) * (lim_dr - lim_st) + lim_st;
    }
}

double calcDeJong(vector<vector<bool>> &v, vector<vector<double>> &val, int cand){
    
    double limDr = 5.12;
    double limSt = -5.12;

    getValues(v, val, cand, limDr, limSt);
    double rez = 0;
    for(int i=0; i<nr_dim; ++i){
        rez += val[cand][i] * val[cand][i];
    }
    return rez;
}

double calcDeJongPartial(vector<vector<bool>> &v, vector<vector<double>> &val, vector<double> &eval, 
                            int cand, int poz, double &newVal){

    double limDr = 5.12;
    double limSt = -5.12;

    int dimPoz = poz / nr_bits_in_dim;
    double rez = eval[cand];

    double vall=v[cand][dimPoz*nr_bits_in_dim];
        for(int j=0; j<nr_bits_in_dim-1; ++j){
            if(v[cand][dimPoz*nr_bits_in_dim + j] == v[cand][dimPoz*nr_bits_in_dim + j+1]){
                vall = vall * 2;
            }
            else {
                vall = vall * 2 + 1;
            }
        }

    // double vall = 0;
    // for(int j=0; j<nr_bits_in_dim; ++j){
    //     vall = vall*2 + v[cand][dimPoz*nr_bits_in_dim+j];
    // }
    vall = (vall / norm_val) * (limDr - limSt) + limSt;

    newVal = vall;
    
    double newRez = rez - (val[cand][dimPoz]* val[cand][dimPoz]) 
                + (vall * vall);

    return newRez;
}


double calcSchwefel(vector<vector<bool>> &v, vector<vector<double>> &val, int cand){
    
    double limDr = 500.0;
    double limSt = -500.0;

    getValues(v, val, cand, limDr, limSt);
    double rez = nr_dim * 418.9829;
    for(int i=0; i<nr_dim; ++i){
        rez -= val[cand][i] * sin(sqrt(abs(val[cand][i])));
    }
    return rez;
}

double calcSchwefelPartial(vector<vector<bool>> &v, vector<vector<double>> &val, vector<double> &eval, 
                            int cand, int poz, double &newVal){
    
    double limDr = 500.0;
    double limSt = -500.0;

    int dimPoz = poz / nr_bits_in_dim;
    double rez = eval[cand];

    double vall=v[cand][dimPoz*nr_bits_in_dim];
        for(int j=0; j<nr_bits_in_dim-1; ++j){
            if(v[cand][dimPoz*nr_bits_in_dim + j] == v[cand][dimPoz*nr_bits_in_dim + j+1]){
                vall = vall * 2;
            }
            else {
                vall = vall * 2 + 1;
            }
        }
    // double vall = 0;
    // for(int j=0; j<nr_bits_in_dim; ++j){
    //     vall = vall*2 + v[cand][dimPoz*nr_bits_in_dim+j];
    // }
    vall = (vall / norm_val) * (limDr - limSt) + limSt;

    newVal = vall;

    double newRez = rez + (val[cand][dimPoz] * sin(sqrt(abs(val[cand][dimPoz]))))
                    - (vall * sin(sqrt(abs(vall))));
    return newRez;

}

double calcMichalewicz(vector<vector<bool>> &v, vector<vector<double>> &val, int cand){
    
    double limDr = M_PI;
    double limSt = 0;

    getValues(v, val, cand, limDr, limSt);
    double rez = 0;
    for(int i=0; i<nr_dim; ++i){
        rez -= (sin(val[cand][i]))*pow((sin(i*(val[cand][i]*val[cand][i])/M_PI)),20);
    }
    return rez;
    
}


double calcMichalewiczPartial(vector<vector<bool>> &v, vector<vector<double>> &val, vector<double> &eval, 
                            int cand, int poz, double &newVal){
    
    double limDr = M_PI;
    double limSt = 0;

    int dimPoz = poz / nr_bits_in_dim;
    double rez = eval[cand];

    double vall=v[cand][dimPoz*nr_bits_in_dim];
        for(int j=0; j<nr_bits_in_dim-1; ++j){
            if(v[cand][dimPoz*nr_bits_in_dim + j] == v[cand][dimPoz*nr_bits_in_dim + j+1]){
                vall = vall * 2;
            }
            else {
                vall = vall * 2 + 1;
            }
        }
    // double vall = 0;

    // for(int j=0; j<nr_bits_in_dim; ++j){
    //     vall = vall*2 + v[cand][dimPoz*nr_bits_in_dim+j];
    // }

    vall = (vall / norm_val) * (limDr - limSt) + limSt;

    newVal = vall;

    double newRez = rez + ((sin(val[cand][dimPoz]))*pow((sin(dimPoz*(val[cand][dimPoz]*val[cand][dimPoz])/M_PI)),20))
                    - ((sin(vall))*pow((sin(dimPoz*(vall*vall)/M_PI)),20));

    return newRez;
}


double calcRastrigin(vector<vector<bool>> &v, vector<vector<double>> &val, int cand){

    double lim_dr = 5.12;
    double lim_st = -5.12;
    double rez;
    getValues(v, val, cand, lim_dr, lim_st);
    //cout << "ok5\n";

    rez = 10*nr_dim;
    for(int i=0; i<nr_dim; ++i){
        rez += val[cand][i]* val[cand][i] - 10*cos(2*M_PI*val[cand][i]);
    }
    return rez;
}

double calcRastriginPartial(vector<vector<bool>> &v, vector<vector<double>> &val, vector<double> &eval, 
                            int cand, int poz, double &newVal){

    double limDr = 5.12;
    double limSt = -5.12;

    int dimPoz = poz / nr_bits_in_dim;
    double rez = eval[cand];
    //rez -= values[cand][dimPoz]* values[cand][dimPoz] - 10*cos(2*M_PI*values[cand][dimPoz]);
    
    double vall=v[cand][dimPoz*nr_bits_in_dim];
        for(int j=0; j<nr_bits_in_dim-1; ++j){
            if(v[cand][dimPoz*nr_bits_in_dim + j] == v[cand][dimPoz*nr_bits_in_dim + j+1]){
                vall = vall * 2;
            }
            else {
                vall = vall * 2 + 1;
            }
        }

    // double vall = 0;

    // //values[cand][dimPoz] = 0;
    // for(int j=0; j<nr_bits_in_dim; ++j){
    //     //values[cand][dimPoz] = values[cand][dimPoz]*2 + pop[cand][dimPoz*20+j];
    //     vall = vall*2 + v[cand][dimPoz*nr_bits_in_dim+j];
    // }
    //values[cand][dimPoz] = (values[cand][dimPoz] / norm_val) * (limDr - limSt) + limSt;
    vall = vall / norm_val * (limDr - limSt) + limSt;

    newVal = vall;

    double newRez = rez - (val[cand][dimPoz]* val[cand][dimPoz] - 10*cos(2*M_PI*val[cand][dimPoz])) 
                + (vall * vall - 10 * cos(2*M_PI*vall));

    //rez += values[cand][dimPoz]* values[cand][dimPoz] - 10*cos(2*M_PI*values[cand][dimPoz]);
    //rez += val * val - 10 * cos(2*M_PI*val);
    return newRez;
}


void genRandCand(vector<bool> &v, std::mt19937 &generator){

    double r;
    for(int i=0; i<nr_biti; ++i){
        
        r = distribution(generator);
        v.push_back( (r > 0.5) );
    }

}


void hillClimb(vector<vector<bool>> &v , vector<vector<double>> &val, vector<double> &eval, int poz, 
        double (*f2) (vector<vector<bool>> &, vector<vector<double>> &, vector<double> &eval, int, int, double &)){
    while(1){
        int fpoz = -1;
        double nv, fnv;
        double ev = eval[poz];
        for(int i=0; i<nr_biti; ++i){
            v[poz][i] = !v[poz][i];
            double vall = f2(v, val, eval, poz, i, nv);
            if(vall < ev){
                ev = vall;
                fpoz = i;
                fnv = nv;
            }
            v[poz][i] = !v[poz][i];
        }
        if(fpoz != -1){
            eval[poz] = ev;
            v[poz][fpoz] = !v[poz][fpoz];
            val[poz][fpoz/nr_bits_in_dim] = fnv;
        }
        else{
            break;
        }
    }
}


// void heapify(int n, int i){

//     int largest = i;
//     int l = 2 * i + 1;
//     int r = 2 * i + 2;

//     if (l < n && eval[l] > eval[largest])
//         largest = l;

//     if (r < n && eval[r] > eval[largest])
//         largest = r;

//     if (largest != i) {
//         swap(eval[i], eval[largest]);
//         swap(pop[i], pop[largest]);
//         heapify(n, largest);
//     }
// }


//  void heapSort(){

//     for (int i = pop_size / 2 - 1; i >= 0; i--)
//         heapify(pop_size, i);

//     for (int i = pop_size - 1; i > 0; i--) {
//         swap(eval[0], eval[i]);
//         swap(pop[0], pop[i]);
//         heapify(i, 0);
//     }
// }



void GA(double (*f1) (vector<vector<bool>> &, vector<vector<double>> &, int), 
        double (*f_part) (vector<vector<bool>> &, vector<vector<double>> &, vector<double> &, int, int, double &)){  

    //#pragma omp parallel for reduction(min:bestVal)

    pid_t pid = getpid();
    std::mt19937 generator(t+pid*10000);
    
    // generate the first generation
    for(int i=0; i < pop_size; ++i){
        vector<bool> c;
        pop.push_back(c);
        genRandCand(pop[i], generator);
        eval[i] = f1(pop, values, i); 
        hillClimb(pop, values, eval, i, f_part);
        //cout << "ok4\n";
    }
    //return;
    //cout << "ok2\n";




    // for(int i=0; i<pop_size; ++i){
    //     cout << eval[i] << '\n';
    // }
    // return;

    // double to_positive = 0;
    // if(f1 == calcMichalewicz){
    //     to_positive = 30;
    // }

    int best_cand_pos;
    gen_nr = 0;
    while(gen_nr != max_gen){ // add break when 100 gen without change
        best_of_gen = 99999999;
        ++gen_nr;
        ofs.clear();
        ofs_values.clear();
        // evaluate the generation
        double sum_eval = 0;

         // once every 200 gen we remove a part (worst) of the population and replace it with new random cand
        // if(gen_nr % 200 == 0){
        //     heapSort();
            
        //     for(int i=pop_size/10; i<pop_size; ++i){
        //         //cout << eval[i] << ' ';
        //         pop[i].clear();
        //         genRandCand(pop[i], generator);
        //         eval[i] = f1(pop, values, i);
        //         hillClimb(pop, values, eval, i, f_part);
        //         //cout << eval[i] << '\n';
        //     }

        //     // double verif_eval;
        //     // vector<vector<double>> verif_values;
        //     // for(int i=0; i<pop_size; ++i){
        //     //     verif_eval = f1(pop, verif_values, i);
        //     //     cout << eval[i] << ' ' << verif_eval << '\n';
        //     // }
        //     // return;
        // }

        for(int i=0; i<pop_size; ++i){
            
            if(eval[i] < best_of_gen){
                best_of_gen = eval[i];
                best_cand_pos = i;
                //best_cand = pop[i];
            }
        }
        //cout << gen_nr << ' ' << best_of_gen << '\n';
        
        if(best_of_gen < best_of_all){
            best_of_all = best_of_gen;
        }
        //cout << best_of_gen << '\n';
        
        double selection_pressure=best_of_gen;
        double eps = 2;
        


        for(int i=0; i<pop_size; ++i){
            sum_eval += 1.0/(eval[i] - selection_pressure  + eps); 
        }
        //cout << "ok3\n";

        //double sum_p = 0;
        // calc the probabilities for each cand
        for(int i=0; i<pop_size; ++i){
            p[i] = (1.0/(eval[i] - selection_pressure  + eps))/sum_eval;
            //sum_p += p[i];
            //cout << eval[i] << ' ' << p[i] << '\n';
        }
        //cout << sum_p << '\n';

        // calc the cumulative probabilities
        q_eval[0] = p[0];
        for(int i=1; i<pop_size; ++i){
            q_eval[i] = q_eval[i-1] + p[i];
        }

        // create the next generation
        double r;
        for(int i=0; i<pop_size; ++i){
            r = distribution(generator);

            int st=0, dr=pop_size-1;
            while(st<dr){
                int mij = (st+dr)/2;
                if(q_eval[mij] > r){
                    dr = mij;
                }
                else{
                    st = mij+1;
                }
            }
            ofs.push_back(pop[dr]);
            ofs_values.push_back(values[dr]);
            ofs_eval[i] = eval[dr];

            // if(gen_nr >= 200)
            //     cout << ofs_eval[i] << '\n';
        }
        
        // if(gen_nr >= 200)
        //     cout << "\n\n";

        // x-over the next generation 
        vector<int> x_over;
        for(int i=0; i<pop_size; ++i){
            r = distribution(generator);
            if(r < pc){
                x_over.push_back(i);
            }
        }

        //cout << "Nr x-over: " << x_over.size() << '\n';
        // if(x_over.size() % 2 == 1){
        //     x_over.pop_back();
        // }
        int size = x_over.size();
        for(int i=0; i < size / 2; ++i){
            r = distribution(generator);
            int partner = r*(x_over.size()-1);

            r = distribution(generator);
            int pos = r * (nr_biti-1)+1;

            for(int j=pos; j<nr_biti; ++j){
                swap(ofs[x_over[x_over.size()-1]][j], ofs[x_over[partner]][j]);
            }
            ofs_eval[x_over[x_over.size()-1]] = f1(ofs, ofs_values, x_over[x_over.size()-1]);
            ofs_eval[x_over[partner]] = f1(ofs, ofs_values, x_over[partner]);

            hillClimb(ofs, ofs_values, ofs_eval, x_over[x_over.size()-1], f_part);
            hillClimb(ofs, ofs_values, ofs_eval, x_over[partner], f_part);

            x_over.pop_back();
            x_over.erase(x_over.begin() + partner);
        }

        // cout << "ofs values after x-over\n";
        // for(int i=0; i<pop_size; ++i){
        //     cout << ofs_eval[i] << '\n';
        // }

        // mutation on the next generation
        for(int i=0; i<nr_biti * pop_size; ++i){
            r = distribution(generator);
            if(r < pm){
                ofs[i/nr_biti][i%nr_biti] = !ofs[i/nr_biti][i%nr_biti];
                double dim_val = 0;
                ofs_eval[i/nr_biti] = f_part(ofs, ofs_values, ofs_eval, i/nr_biti, i%nr_biti, dim_val);
                ofs_values[i/nr_biti][(i%nr_biti)/nr_bits_in_dim] = dim_val;
            }
        }

        // for(int i=0; i<pop_size; ++i){
        //      hillClimb(ofs, ofs_values, ofs_eval, i, f_part);
        // }


        // cout << "ofs values after x-over and mutation\n";
        // for(int i=0; i<pop_size; ++i){
            

        //     cout << ofs_eval[i] << '\n';
        // }
         //return;
        

        pop = ofs;
        values = ofs_values;
        eval = ofs_eval;

        //cout << best_of_gen << '\n';
    }
    

}


int main(){

    nr_dim=5; // 30
    //nr_cand = 200;
    nr_bits_in_dim = 20;
    nr_biti = nr_dim * nr_bits_in_dim; // at least 27 for schwefel func - 30 for good measure :))
    norm_val = (1<<nr_bits_in_dim) - 1;

    pop_size = 400; // 50 - 100
    pc = .3;
    pm = .002;
    max_gen = 500;

    
    t1 = chrono::high_resolution_clock::now();
    double (*f1) (vector<vector<bool>> &, vector<vector<double>> &, int);
    double (*f_part) (vector<vector<bool>> &, vector<vector<double>> &, vector<double> &, int, int, double &);

    // f1 = calcMichalewicz;
    // f_part = calcMichalewiczPartial;

    // f1 = calcRastrigin;
    // f_part = calcRastriginPartial;

    // f1 = calcSchwefel;
    // f_part = calcSchwefelPartial;

    f1 = calcDeJong;
    f_part = calcDeJongPartial;
    
    //cout << "ok\n";
    GA(f1, f_part);

    
    t2 = chrono::high_resolution_clock::now();

    int best_val1, best_val2;

    best_val1 = (int)best_of_all;
    best_val2 = (int)(abs((best_of_all-best_val1)*100000));
    int c = best_val2;
    int maxi = 4;

    g << best_val1 << '.';

    while(c < 10000 && maxi){
        g << '0';
        c *= 10;
        maxi --;
    }
    g << best_val2 << ' ';
    
    chrono::milliseconds time_span = chrono::duration_cast<chrono::milliseconds>(t2-t1);
    g << time_span.count() << "\n";
}