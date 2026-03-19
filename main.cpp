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

ifstream f1("C:/Users/Rares/Documents/Facultate/An2/Sem1/AG/H_3/ex1.txt");
ifstream f2("C:/Users/Rares/Documents/Facultate/An2/Sem1/AG/H_3/ex2.txt");
ifstream f3("C:/Users/Rares/Documents/Facultate/An2/Sem1/AG/H_3/ex3.txt");
ifstream f4("C:/Users/Rares/Documents/Facultate/An2/Sem1/AG/H_3/ex4.txt");
ofstream g("C:/Users/Rares/Documents/Facultate/An2/Sem1/AG/H_3/res.txt", std::ios_base::app);

#define MAX_L 22505
#define MAX_C 124005

vector<vector<int>> pop_cl_sat, ofs_cl_sat, sol_cl_sat(6, vector<int>(MAX_C)); // cati literali satisfac clauza
vector<vector<bool>> sol_cand(6, vector<bool>(MAX_L)); // sol_cand[i] - val lit i
vector<int>eval_sol_cand(6);
vector<vector<int>> cl_cu_lit(MAX_C); // literalii pt fiecare clauza
vector<int> cl_in_care_apare_lit[MAX_L*2]; // ficare var (l si -l) in ce clauze se afla
vector<int> pattern(MAX_L);

vector<vector<bool>> pop, ofs;
vector<double>q_eval(2000), p(2000); // eval = nr de clauze satisfacute
vector<int> pop_eval(2000), ofs_eval(2000);
int nr_cl_sat;
int nr_l, nr_c;

double best_of_all = 0, best_of_gen = 0;
int nr_biti, pop_size, max_gen, gen_nr;
double pc, pm;

int nr_input;

chrono::high_resolution_clock::time_point t1,t2;

unsigned t = chrono::high_resolution_clock::now().time_since_epoch().count();

std::uniform_real_distribution<> distribution(0.0, 1.0);


//functions
void genRandCand(vector<bool> &v, std::mt19937 &generator);
void hillClimb(vector<vector<bool>> &pop, vector<vector<int>> &cl_sat, vector<int> &eval, int poz);
int calc_cand(vector<vector<bool>> &pop, vector<vector<int>> &cl_sat, int cand);
int calc_cand_part(vector<vector<bool>> &pop, vector<vector<int>> &cl_sat, int cand, int lit, int rez);


void genRandCand(vector<bool> &v, std::mt19937 &generator){

    double r;
    v.push_back(0); // poz 0 in vector;
    for(int i=1; i<=nr_biti; ++i){
        if(pattern[i] == -1){
            r = distribution(generator);
            v.push_back( (r > 0.5) );
        }
        else{
            v.push_back(pattern[i]);
        }
    }

}


void hillClimb(vector<vector<bool>> &pop, vector<vector<int>> &cl_sat, vector<int> &eval, int poz, int step_nr){
    while(step_nr){
        --step_nr;
        int fpoz = -1;
        int ev = eval[poz];
        int check = ev;
        for(int i=1; i<=nr_biti; ++i){
            if(pattern[i] != -1){
                continue;
            }
            //v[poz][i] = !v[poz][i];
            eval[poz] = calc_cand_part(pop, cl_sat, poz, i, eval[poz]);
            int val2 = eval[poz];
            eval[poz] = calc_cand_part(pop, cl_sat, poz, i, eval[poz]);
            if(eval[poz] != check){
                cout << "Eroare la calc_cand_part() in hillClimb: " << ev << ' ' << val2 << '\n';
                return;
            }
            if(val2 > ev){
                ev = val2;
                fpoz = i;
            }
        }
        if(fpoz != -1){
            eval[poz] = calc_cand_part(pop, cl_sat, poz, fpoz, eval[poz]);

            if(eval[poz] != ev){
                cout << "Eroare in HillClimb\n";
                return;
            }
        }
        else{
            break;
        }
    }
}


int calc_cand(vector<vector<bool>> &pop, vector<vector<int>> &cl_sat, int cand){
    int rez = 0;
    while(cl_sat.size() <= cand){
        vector<int> v(nr_c+1);
        cl_sat.push_back(v);
    }
    for(int i=0; i<nr_c; ++i){
        cl_sat[cand][i] = 0;
    }

    for(int i=1; i<=nr_biti; ++i){
        int poz = ((int)pop[cand][i] * 2 - 1) * i + nr_l;
        for(int j=0; j<cl_in_care_apare_lit[poz].size(); ++j){
            cl_sat[cand][cl_in_care_apare_lit[poz][j]] ++;
            if(cl_sat[cand][cl_in_care_apare_lit[poz][j]] == 1){
                rez++;
            }
        }
    }
    if(rez > best_of_all){
        //sol_cand = pop[cand];
        best_of_all = rez;
    }
    if(eval_sol_cand[4] < rez){
        eval_sol_cand[4] = rez;
        sol_cand[4] = pop[cand];
        sol_cl_sat[4] = pop_cl_sat[cand];
    }
    for(int i=3; i>=0; --i){
        if(eval_sol_cand[i] < rez){
            swap(eval_sol_cand[i+1], eval_sol_cand[i]);
            swap(sol_cand[i+1], sol_cand[i]);
            swap(sol_cl_sat[i+1], sol_cl_sat[i]);
        }
    }
    return rez;
}

int calc_cand_part(vector<vector<bool>> &pop, vector<vector<int>> &cl_sat, int cand, int lit, int rez){
    
    int poz = (pop[cand][lit] * 2 - 1) * lit + nr_l;
    for(int j=0; j<cl_in_care_apare_lit[poz].size(); ++j){
        cl_sat[cand][cl_in_care_apare_lit[poz][j]] --;
        if(cl_sat[cand][cl_in_care_apare_lit[poz][j]] == 0){
            rez--;
        }
    }

    poz = ((1-pop[cand][lit]) * 2 - 1) * lit + nr_l;
    for(int j=0; j<cl_in_care_apare_lit[poz].size(); ++j){
        cl_sat[cand][cl_in_care_apare_lit[poz][j]] ++;
        if(cl_sat[cand][cl_in_care_apare_lit[poz][j]] == 1){
            rez++;
        }
    }
    pop[cand][lit] = 1-pop[cand][lit];

    if(rez > best_of_all){
        //sol_cand = pop[cand];
        best_of_all = rez;
    }
     if(eval_sol_cand[4] < rez){
        eval_sol_cand[4] = rez;
        sol_cand[4] = pop[cand];
        sol_cl_sat[4] = pop_cl_sat[cand];
    }
    for(int i=3; i>=0; --i){
        if(eval_sol_cand[i] < rez){
            swap(eval_sol_cand[i+1], eval_sol_cand[i]);
            swap(sol_cand[i+1], sol_cand[i]);
            swap(sol_cl_sat[i+1], sol_cl_sat[i]);
        }
    }

    return rez;
}

void GA(){  

    pid_t pid = getpid();
    std::mt19937 generator(t+pid*10000);
    int best_before_HC = 0;
    int best_after_HC = 0;
    // generate the first generation
    for(int i=0; i < pop_size; ++i){
        vector<bool> c;
        pop.push_back(c);
        genRandCand(pop[i], generator);
        pop_eval[i] = calc_cand(pop, pop_cl_sat, i); 
        // if(best_before_HC < pop_eval[i]){
        //     best_before_HC = pop_eval[i];
        // }
        //g << pop_eval[i] << ' ';
        hillClimb(pop, pop_cl_sat, pop_eval, i, 100);

        // if(best_after_HC < pop_eval[i]){
        //     best_after_HC = pop_eval[i];
        // }
        //g << pop_eval[i] << '\n';

        //cout << "ok4\n";
    }
    // cout << best_before_HC << ' ' << best_after_HC << '\n';
    // fflush(0);
    //return;

    // for(int i=0; i<pop_size; ++i){
    //     cout << eval[i] << '\n';
    // }
    // return;


    int best_cand_pos;
    int worst_of_gen;
    gen_nr = 0;
    while(gen_nr != max_gen){
        best_of_gen = 0;
        worst_of_gen = 9999999;
        ++gen_nr;
        ofs.clear();
        ofs_cl_sat.clear();
        // evaluate the generation
        double sum_eval = 0;

        for(int i=0; i<pop_size; ++i){
            
            if(pop_eval[i] > best_of_gen){
                best_of_gen = pop_eval[i];
                best_cand_pos = i;
            }
            if(pop_eval[i] < worst_of_gen){
                worst_of_gen = pop_eval[i];
            }
        }
        //cout << gen_nr << ' ' << best_of_gen << '\n';
        
        if(best_of_gen > best_of_all){
            best_of_all = best_of_gen;
        }
        if(best_of_all == nr_c){
            return;
        }
        
        double selection_pressure=worst_of_gen;
        double eps = 3 * (best_of_gen - worst_of_gen);
        // if(gen_nr == 1 || gen_nr == max_gen){
        //     cout << best_of_gen << ' ' << worst_of_gen << '\n';
        // }
        //double eps = 1000;

        for(int i=0; i<pop_size; ++i){
            sum_eval += pop_eval[i] - selection_pressure + eps; 
        }
        // calc the probabilities for each cand
        for(int i=0; i<pop_size; ++i){
            p[i] = (pop_eval[i] - selection_pressure  + eps)/sum_eval;
        }

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
            ofs_eval[i] = pop_eval[dr];
            ofs_cl_sat.push_back(pop_cl_sat[dr]);

        }
        // x-over the next generation 
        vector<int> x_over;
        for(int i=0; i<pop_size; ++i){
            r = distribution(generator);
            if(r < pc){
                x_over.push_back(i);
            }
        }

        int size = x_over.size();
        for(int i=0; i < size / 2; ++i){
            r = distribution(generator);
            int partner = r*(x_over.size()-1);

            r = distribution(generator);
            int pos = r * (nr_biti-1)+2;

            // g << pos << '\n';

            // for(int j=1; j<=nr_biti; ++j){
            //     g << ofs[x_over[x_over.size()-1]][j];
            // }
            // g << '\n';
            // for(int j=1; j<=nr_biti; ++j){
            //     g << ofs[x_over[partner]][j];
            // }
            // g << '\n';

            for(int j=pos; j<=nr_biti; ++j){
                swap(ofs[x_over[x_over.size()-1]][j], ofs[x_over[partner]][j]);
            }

            // for(int j=1; j<=nr_biti; ++j){
            //     g << ofs[x_over[x_over.size()-1]][j];
            // }
            // g << '\n';
            // for(int j=1; j<=nr_biti; ++j){
            //     g << ofs[x_over[partner]][j];
            // }
            // g << '\n';

            ofs_eval[x_over[x_over.size()-1]] = calc_cand(ofs, ofs_cl_sat, x_over[x_over.size()-1]);
            ofs_eval[x_over[partner]] = calc_cand(ofs, ofs_cl_sat, x_over[partner]);

            // cout << ofs_eval[x_over[x_over.size()-1]] << ' ' << ofs_eval[x_over[partner]] << '\n';
            // return;

            hillClimb(ofs, ofs_cl_sat, ofs_eval, x_over[x_over.size()-1], 50);
            hillClimb(ofs, ofs_cl_sat, ofs_eval, x_over[partner], 50);

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
                ofs_eval[i/nr_biti] = calc_cand_part(ofs, ofs_cl_sat, i/nr_biti, (i%nr_biti) + 1, ofs_eval[i/nr_biti]);
            }
        }


        // cout << "ofs values after x-over and mutation\n";
        // for(int i=0; i<pop_size; ++i){
        //     cout << ofs_eval[i] << '\n';
        // }
        // cout << '\n';
        // for(int i=0; i<pop_size; ++i){
        //     cout << pop_eval[i] << '\n';
        // }
         //return;
        

        pop = ofs;
        pop_eval = ofs_eval;
        pop_cl_sat = ofs_cl_sat;

        //cout << best_of_gen << '\n';
    }
    // cout << "intru in hill_climb cu top 5, best_of_all: " << best_of_all << "\n";
    // fflush(0);
    for(int i=0; i<5; ++i){
        pop[0] = sol_cand[i];
        pop_eval[0] = eval_sol_cand[i];
        pop_cl_sat[0] = sol_cl_sat[i];
        hillClimb(pop, pop_cl_sat, pop_eval, 0, 25000);
        if(pop_eval[0] > best_of_all){
            best_of_all = pop_eval[0];
        }
    }
}

void citire_input(){
    // cout << "Incep citirea\n";
    // fflush(0);
    int cnt = 0;
    if(nr_input == 1){
        f1 >> nr_l >> nr_c;
    }
    else if(nr_input == 2){
        f2 >> nr_l >> nr_c;
    }
    else if(nr_input == 3){
        f3 >> nr_l >> nr_c;
    }
    else{
        f4 >> nr_l >> nr_c;
    }
    //cout << nr_l << ' ' << nr_c << '\n';

    int a, b, c;
    for(int i=0; i<nr_c; ++i){
        vector<int> new_cl;
        if(nr_input == 1){
        f1 >> a;
    }
    else if(nr_input == 2){
        f2 >> a;
    }
    else if(nr_input == 3){
        f3 >> a;
    }
    else{
        f4 >> a;
    }
        b = a;
        while(a){
            new_cl.push_back(a);
            cl_in_care_apare_lit[a+nr_l].push_back(i);
            if(nr_input == 1){
                f1 >> a;
            }
            else if(nr_input == 2){
                f2 >> a;
            }
            else if(nr_input == 3){
                f3 >> a;
            }
            else{
                f4 >> a;
            }
        }

        if(new_cl.size() == 1){
            pattern[abs(b)] = (b + abs(b)) / (2*b);
            ++cnt;
        }

        cl_cu_lit[i] = new_cl;
        

    }
    // cout << cnt << ' ';

     
     
     for(int i=1; i<= nr_l; ++i){
         int x = cl_in_care_apare_lit[i + nr_l].size();
         int y = cl_in_care_apare_lit[-i + nr_l].size(); 
         if(x == 0){
            pattern[i] = 0;
            ++cnt;
        }
        else if(y == 0){
            pattern[i] = 1;
            ++cnt ;
        }
    }
    //cout << cnt << ' ';
    
    // int nr = 0;
    // for(int i=0; i<=nr_c; ++i){
    //     for(int j=0; j<cl_cu_lit[i].size(); ++j){
    //         if(pattern[abs(cl_cu_lit[i][j] != -1)]){
    //             if((pattern[abs(cl_cu_lit[i][j])]*2 - 1) * abs(cl_cu_lit[i][j]) != cl_cu_lit[i][j]){
    //                 cl_cu_lit[i].erase(cl_cu_lit[i].begin() + j);
    //                 j--;
    //                 if(cl_cu_lit[i].size() == 1){

    //                 }
    //             }
    //         }
    //     }
    // }

    // int nr = 0;
    // for(int i=1; i<= nr_l; ++i){
    //     if(pattern[i] != -1){
    //         int poz = (-(pattern[i] * 2 -1) *i + nr_l);
    //         for(int j = 0; j < cl_in_care_apare_lit[poz].size(); ++j){
    //             int cl = cl_in_care_apare_lit[poz][j];
    //             for(int k = 0; k < cl_cu_lit[cl].size(); ++k){
    //                 if(poz == cl_cu_lit[cl][k]){
    //                     cl_cu_lit[cl].erase(cl_cu_lit[cl].begin() + k);
    //                     if(cl_cu_lit.size() == 1){
    //                         pattern[abs(cl_cu_lit[cl][0])] = (cl_cu_lit[cl][0] + abs(cl_cu_lit[cl][0])) / (2*cl_cu_lit[cl][0]);
    //                         ++ cnt;
    //                     }
    //                     break;
    //                 }
    //             }
    //         }
    //     }
    // }


    //cout << cnt << " nr lit in pattern\n";
    //fflush(0);
    // cout << "Am terminat citirea\n";
    // fflush(0);

}

int main(int argc, char* argv[]){
    nr_input = atoi(argv[1]);
    // cout  << nr_input << '\n';
    // fflush(0);
    for(int i=0; i<MAX_L; ++i){
        pattern[i] = -1;
    }

    citire_input();

    nr_biti = nr_l;
    pop_size = 100; // 100 - 200
    pc = .4;
    pm = 1/nr_biti;
    max_gen = 200; // 1000 - 2000

    t1 = chrono::high_resolution_clock::now();

    GA();
    
    t2 = chrono::high_resolution_clock::now();


    g << best_of_all << ' ';
    
    chrono::milliseconds time_span = chrono::duration_cast<chrono::milliseconds>(t2-t1);
    g << time_span.count() << ' ' << nr_input <<"\n";
}