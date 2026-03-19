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


vector<vector<bool>> pop, ofs, meta_pop, meta_ofs;
vector<vector<double>> values, ofs_values, meta_val, meta_ofs_val;
vector<double>eval(200), q_eval(200), p(200), ofs_eval(200), meta_eval(100), meta_q_eval(100), meta_p(100), meta_ofs_eval(100);

double best_of_all = 999999999, best_of_gen, meta_best=999999999, meta_best_of_gen;
int nr_dim, nr_biti, nr_bits_in_dim, nr_cand, pop_size, max_gen, gen_nr;
int meta_nr_dim, meta_nr_biti, meta_pop_size, meta_max_gen, meta_gen_nr, cnt;
double best_pc, best_pm, best_sel_press;
double pc, pm, sel_press;
double meta_pc, meta_pm;
int norm_val;
chrono::high_resolution_clock::time_point t1,t2;

unsigned t = chrono::high_resolution_clock::now().time_since_epoch().count();

std::uniform_real_distribution<> distribution(0.0, 1.0);


void getValues(vector<vector<bool>> &v, vector<vector<double>> &val, int cand, double lim_dr, double lim_st){
    if(val.size() <= cand){
        val.push_back(vector<double> (nr_dim));
    }
    for(int i=0; i<nr_dim; ++i){
        val[cand][i]=0;
        for(int j=0; j<nr_bits_in_dim; ++j){
            val[cand][i] = val[cand][i] * 2 + v[cand][i*nr_bits_in_dim + j];
        }
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

    double vall = 0;
    for(int j=0; j<nr_bits_in_dim; ++j){
        vall = vall*2 + v[cand][dimPoz*nr_bits_in_dim+j];
    }
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

    double vall = 0;
    for(int j=0; j<nr_bits_in_dim; ++j){
        vall = vall*2 + v[cand][dimPoz*nr_bits_in_dim+j];
    }
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

    double vall = 0;
    for(int j=0; j<nr_bits_in_dim; ++j){
        vall = vall*2 + v[cand][dimPoz*nr_bits_in_dim+j];
    }
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

    double vall = 0;

    //values[cand][dimPoz] = 0;
    for(int j=0; j<nr_bits_in_dim; ++j){
        //values[cand][dimPoz] = values[cand][dimPoz]*2 + pop[cand][dimPoz*20+j];
        vall = vall*2 + v[cand][dimPoz*nr_bits_in_dim+j];
    }
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


void genRandCand_GA(vector<bool> &v, std::mt19937 &generator){

    double r;
    for(int i=0; i<meta_nr_biti; ++i){
        
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


void GA(double (*f1) (vector<vector<bool>> &, vector<vector<double>> &, int), 
        double (*f_part) (vector<vector<bool>> &, vector<vector<double>> &, vector<double> &, int, int, double &), int rand){  

    //#pragma omp parallel for reduction(min:bestVal)

    pid_t pid = getpid();
    std::mt19937 generator(t+pid*10000 + rand * 100000);
    
    // generate the first generation
    //if(nr_dim == 30) cout <<"inceput in GA\n";
    for(int i=0; i < pop_size; ++i){
        vector<bool> c;
        pop.push_back(c);
        genRandCand(pop[i], generator);
        //if(nr_dim == 30)return;
        eval[i] = f1(pop, values, i); 
        hillClimb(pop, values, eval, i, f_part);
        //cout << "ok4\n";
    }
    //if(nr_dim == 30) cout <<"dupa pop initiala in GA\n";
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
    //if(nr_dim == 30) cout <<"intru in while in GA\n";
    while(gen_nr != max_gen){ // add break when 100 gen without change
        best_of_gen = 99999999;
        ++gen_nr;
        ofs.clear();
        ofs_values.clear();
        // evaluate the generation
        double sum_eval = 0;

       
        for(int i=0; i<pop_size; ++i){
            
            if(eval[i] < best_of_gen){
                best_of_gen = eval[i];
                best_cand_pos = i;
                //best_cand = pop[i];
            }
        }
        //if(nr_dim == 30) cout << gen_nr << ' ' << best_of_gen << '\n';
        
        if(best_of_gen < best_of_all){
            best_of_all = best_of_gen;
        }
        //cout << best_of_gen << '\n';
        
        double selection_pressure=(sel_press + 1)*best_of_gen;
         double eps = 1e-5;


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

            // hillClimb(ofs, ofs_values, ofs_eval, x_over[x_over.size()-1], f_part);
            // hillClimb(ofs, ofs_values, ofs_eval, x_over[partner], f_part);

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
                //hillClimb(ofs, ofs_values, ofs_eval, i/nr_biti, f_part);
            }
        }


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

    pop.clear();
    values.clear();
    ofs.clear();
    ofs_values.clear();
}



void getValues_GA(vector<vector<bool>> &v, vector<vector<double>> &val, int cand, 
    double dr1, double st1, double dr2, double st2, double dr3, double st3){
    if(val.size() <= cand){
        val.push_back(vector<double> (meta_nr_dim));
    }

    int i;
    int nr, bpos, meta_norm_val;
    i=0;
    val[cand][i]=0;
    nr = 5;
    bpos = 0;
    meta_norm_val = (1<<nr) - 1;
 
    for(int j=0; j<nr; ++j){
        val[cand][i] = val[cand][i] * 2 + v[cand][bpos + j];
    }
    val[cand][i] = (val[cand][i] / meta_norm_val) * (dr1 - st1) + st1;
    
    i=1;
    val[cand][i]=0;
    nr = 7;
    bpos = 5;
    meta_norm_val = (1<<nr) - 1;

    for(int j=0; j<nr; ++j){
        val[cand][i] = val[cand][i] * 2 + v[cand][bpos + j];
    }
    val[cand][i] = (val[cand][i] / meta_norm_val) * (dr2 - st2) + st2;

    i=2;
    val[cand][i]=0;
    nr = 4;
    bpos = 12;
    meta_norm_val = (1<<nr) - 1;

    for(int j=0; j<nr; ++j){
        val[cand][i] = val[cand][i] * 2 + v[cand][bpos + j];
    }
    val[cand][i] = (val[cand][i] / meta_norm_val) * (dr3 - st3) + st3;



}

double calc_GA(vector<vector<bool>> &v, vector<vector<double>> &val, int cand){

    double lim_1_dr = 1;
    double lim_1_st = 0.1;
    double lim_2_dr = 0.1;
    double lim_2_st = 0;
    double lim_3_dr = 1;
    double lim_3_st = 0;

    double rez;
    getValues_GA(v, val, cand, lim_1_dr, lim_1_st, lim_2_dr, lim_2_st, lim_3_dr, lim_3_st);
    //cout << "ok5\n";

    nr_dim=30; // 30
    //nr_cand = 200;
    nr_bits_in_dim = 20;
    nr_biti = nr_dim * nr_bits_in_dim; 
    norm_val = (1<<nr_bits_in_dim) - 1;

    pop_size = 10; // 50 - 100
    pc = meta_val[cand][0];
    pm = meta_val[cand][1];
    sel_press =  meta_val[cand][2];
    max_gen = 5;

    double (*f1) (vector<vector<bool>> &, vector<vector<double>> &, int);
    double (*f_part) (vector<vector<bool>> &, vector<vector<double>> &, vector<double> &, int, int, double &);
    f1 = calcMichalewicz;
    f_part = calcMichalewiczPartial;
    
    //cout << "ok\n";
    best_of_all = 0;

    ++cnt;
    //cout << "intram in GA\n";
    GA(f1, f_part, cnt);
    //cout << "am iesit din GA\n";
    rez = best_of_all;
    if(rez < meta_best){
        meta_best = rez;
        best_pc = meta_val[cand][0];
        best_pm = meta_val[cand][1];
        best_sel_press = meta_val[cand][2];
    }
    return rez;
}

void hillClimb_GA(vector<vector<bool>> &v , vector<vector<double>> &val, vector<double> &eval, int poz,
    double (*ga) (vector<vector<bool>> &, vector<vector<double>> &, int)){
        //cout << "intru in while in hill Climb GA\n";
        //cout << eval[poz] << '\n';
    while(1){
        int fpoz = -1;
        double ev = eval[poz];
        for(int i=0; i<meta_nr_biti; ++i){
            v[poz][i] = !v[poz][i];
            double vall = ga(v, val, poz);
            if(vall < ev){
                ev = vall;
                fpoz = i;
            }
            v[poz][i] = !v[poz][i];
        }
        if(fpoz != -1){
            v[poz][fpoz] = !v[poz][fpoz];
            eval[poz] = ev;
            //cout << ev << '\n';
            
        }
        else{
            break;
        }
    }
}

void meta_GA(double (*ga) (vector<vector<bool>> &, vector<vector<double>> &, int)){

    pid_t pid = getpid();
    std::mt19937 generator(t+pid*10000);
    
    // generate the first generation
    for(int i=0; i < meta_pop_size; ++i){
        vector<bool> c;
        meta_pop.push_back(c);
        genRandCand_GA(meta_pop[i], generator);
        meta_eval[i] = ga(meta_pop, meta_val, i); 
        //hillClimb_GA(meta_pop, meta_val, meta_eval, i, ga);
        //cout << '\n';
        //cout << "ok4\n";
    }

    int best_cand_pos;
    for(int i=0; i<meta_pop_size; ++i){
            
        if(meta_eval[i] < meta_best_of_gen){
            meta_best_of_gen = meta_eval[i];
            best_cand_pos = i;
            //best_cand = pop[i];
        }
    }
    cout << meta_gen_nr << ' ' << meta_best_of_gen << '\n';
    
    if(meta_best_of_gen < meta_best){
        meta_best = meta_best_of_gen;
    }

    meta_gen_nr = 0;
    while(meta_gen_nr != meta_max_gen){ // add break when 100 gen without change
        meta_best_of_gen = 99999999;
        ++meta_gen_nr;
        meta_ofs.clear();
        meta_ofs_val.clear();
        // evaluate the generation
        double meta_sum_eval = 0;

        
        //cout << best_of_gen << '\n';
        
        double selection_pressure= meta_best_of_gen;
         double eps = 1e-5;


        for(int i=0; i<meta_pop_size; ++i){
            meta_sum_eval += 1.0/(meta_eval[i] - selection_pressure  + eps); 
        }
        //cout << "ok3\n";

        //double sum_p = 0;
        // calc the probabilities for each cand
        for(int i=0; i<meta_pop_size; ++i){
            meta_p[i] = (1.0/(meta_eval[i] - selection_pressure  + eps))/meta_sum_eval;
            //sum_p += p[i];
            //cout << eval[i] << ' ' << p[i] << '\n';
        }
        //cout << sum_p << '\n';

        // calc the cumulative probabilities
        meta_q_eval[0] = meta_p[0];
        for(int i=1; i<meta_pop_size; ++i){
            meta_q_eval[i] = meta_q_eval[i-1] + meta_p[i];
        }

        // create the next generation
        double r;
        for(int i=0; i<meta_pop_size; ++i){
            r = distribution(generator);

            int st=0, dr=meta_pop_size-1;
            while(st<dr){
                int mij = (st+dr)/2;
                if(meta_q_eval[mij] > r){
                    dr = mij;
                }
                else{
                    st = mij+1;
                }
            }
            meta_ofs.push_back(meta_pop[dr]);
            meta_ofs_val.push_back(meta_val[dr]);
            meta_ofs_eval[i] = meta_eval[dr];

            // if(gen_nr >= 200)
            //     cout << ofs_eval[i] << '\n';
        }
        
        // if(gen_nr >= 200)
        //     cout << "\n\n";

        // x-over the next generation 
        vector<int> x_over;
        for(int i=0; i<meta_pop_size; ++i){
            r = distribution(generator);
            if(r < meta_pc){
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
            int pos = r * (meta_nr_biti-1)+1;

            for(int j=pos; j<meta_nr_biti; ++j){
                swap(meta_ofs[x_over[x_over.size()-1]][j], meta_ofs[x_over[partner]][j]);
            }
            meta_ofs_eval[x_over[x_over.size()-1]] = ga(meta_ofs, meta_ofs_val, x_over[x_over.size()-1]);
            meta_ofs_eval[x_over[partner]] = ga(meta_ofs, meta_ofs_val, x_over[partner]);

            // hillClimb_GA(meta_ofs, meta_ofs_val, meta_ofs_eval, x_over[x_over.size()-1], ga);
            // hillClimb_GA(meta_ofs, meta_ofs_val, meta_ofs_eval, x_over[partner], ga);

            x_over.pop_back();
            x_over.erase(x_over.begin() + partner);
        }

        // cout << "ofs values after x-over\n";
        // for(int i=0; i<pop_size; ++i){
        //     cout << ofs_eval[i] << '\n';
        // }

        // mutation on the next generation
        for(int i=0; i<meta_nr_biti * meta_pop_size; ++i){
            r = distribution(generator);
            if(r < meta_pm){
                meta_ofs[i/meta_nr_biti][i%meta_nr_biti] = !meta_ofs[i/meta_nr_biti][i%meta_nr_biti];
                double dim_val = 0;
                meta_ofs_eval[i/meta_nr_biti] = ga(meta_ofs, meta_ofs_val, i/meta_nr_biti);
            }
        }

        for(int i=0; i<meta_pop_size; ++i){
            
            if(meta_eval[i] < meta_best_of_gen){
                meta_best_of_gen = meta_eval[i];
                best_cand_pos = i;
                //best_cand = pop[i];
            }
        }
        cout << meta_gen_nr << ' ' << meta_best_of_gen << '\n';
        
        if(meta_best_of_gen < meta_best){
            meta_best = meta_best_of_gen;
        }


        // cout << "ofs values after x-over and mutation\n";
        // for(int i=0; i<pop_size; ++i){
            

        //     cout << ofs_eval[i] << '\n';
        // }
         //return;
        

        meta_pop = meta_ofs;
        meta_val = meta_ofs_val;
        meta_eval = meta_ofs_eval;
    }

    // try the GA with the best params found
    
    nr_dim=30; // 30
    //nr_cand = 200;
    nr_bits_in_dim = 20;
    nr_biti = nr_dim * nr_bits_in_dim; 
    norm_val = (1<<nr_bits_in_dim) - 1;
    
    pop_size = 200; // 50 - 100
    pc = best_pc;
    pm = best_pm ;//* 300 / nr_biti;
    sel_press =  best_sel_press;
    max_gen = 1000;
    
    double (*f1) (vector<vector<bool>> &, vector<vector<double>> &, int);
    double (*f_part) (vector<vector<bool>> &, vector<vector<double>> &, vector<double> &, int, int, double &);
    f1 = calcMichalewicz;
    f_part = calcMichalewiczPartial;
    
    //cout << "ok\n";
    best_of_all = 0;
    
    ++cnt;
    //g << 1 << '\n';
    cout << "started the long run\n";
    g << best_pc << ' ' << best_pm << ' ' << best_sel_press << ' ' << best_of_all << '\n';
    //cout << best_pc << ' ' << best_pm << ' ' << best_sel_press << ' ' << best_of_all << '\n';
    //return;
    GA(f1, f_part, cnt);
    cout << best_of_all << '\n';
    
}


int main(){

    meta_nr_dim = 3;
    meta_nr_biti = 16;
    meta_pop_size = 100; //100
    meta_max_gen = 10; // 500
    meta_pc = 0.5;
    meta_pm = 1.0 / meta_nr_biti;

    
    t1 = chrono::high_resolution_clock::now();
    
    double (*ga) (vector<vector<bool>> &, vector<vector<double>> &, int);
    ga = calc_GA;

    g << "start the meta_GA\n";
    meta_GA(ga);
    
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
