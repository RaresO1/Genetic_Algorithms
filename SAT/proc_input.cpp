#include <iostream>
#include <fstream>
#include <cstring>

using namespace std;

//ifstream f("C:/Users/Rares/Documents/Facultate/An2/Sem1/AG/H_3/Instante_extracted/aim/aim-200-6_0-yes1-1.cnf");
//ifstream f("C:/Users/Rares/Documents/Facultate/An2/Sem1/AG/H_3/Instante_extracted/Bejing/3blocks.cnf");

//ifstream f("C:/Users/Rares/Documents/Facultate/An2/Sem1/AG/H_3/Instante_extracted/uf20-91/uf20-01.cnf");
//ifstream f("C:/Users/Rares/Documents/Facultate/An2/Sem1/AG/H_3/Instante_extracted/uf75-325/uf75-01.cnf");
//ifstream f("C:/Users/Rares/Documents/Facultate/An2/Sem1/AG/H_3/Instante_extracted/uf175-753/uf175-01.cnf");

// S
ifstream f("C:/Users/Rares/Documents/Facultate/An2/Sem1/AG/H_3/Instante_extracted/uf50-218/uf50-01.cnf");
//ifstream f("C:/Users/Rares/Documents/Facultate/An2/Sem1/AG/H_3/Instante_extracted/uf100-430/uf100-01.cnf");
//ifstream f("C:/Users/Rares/Documents/Facultate/An2/Sem1/AG/H_3/Instante_extracted/uf150-645/uf150-01.cnf");

// M
//ifstream f("C:/Users/Rares/Documents/Facultate/An2/Sem1/AG/H_3/Instante_extracted/uf250-1065/uf250-01.cnf");
//ifstream f("C:/Users/Rares/Documents/Facultate/An2/Sem1/AG/H_3/Instante_extracted/f/f600.cnf");
//ifstream f("C:/Users/Rares/Documents/Facultate/An2/Sem1/AG/H_3/Instante_extracted/f/f1000.cnf");
//ifstream f("C:/Users/Rares/Documents/Facultate/An2/Sem1/AG/H_3/Instante_extracted/f/f2000.cnf");


// L
//ifstream f("C:/Users/Rares/Documents/Facultate/An2/Sem1/AG/H_3/Instante_extracted/Bejing/enddr2-10-by-5-8.cnf");
//ifstream f("C:/Users/Rares/Documents/Facultate/An2/Sem1/AG/H_3/Instante_extracted/Bejing/ewddr2-10-by-5-1.cnf");
//ifstream f("C:/Users/Rares/Documents/Facultate/An2/Sem1/AG/H_3/Instante_extracted/Bejing/ewddr2-10-by-5-8.cnf");

ofstream g("C:/Users/Rares/Documents/Facultate/An2/Sem1/AG/H_3/ex3.txt");

int main(){

    int nr_l=0, nr_c=0;
    int a,b,c,x;

    char buffer[1024];
    f.getline(buffer, 1024);
    while(buffer[0] != 'p'){
        f.getline(buffer, 1024);
    }
    if(buffer[0] == 'p'){
        cout << "ok\n";
        fflush(0);
    }
    else{
        cout << "nu am ajuns la enunt\n";
        return 1;
    }

    int nr=0;
    for(int i = 6; i<sizeof(buffer) && buffer[i] != '\n' && buffer[i] != '\0' ; i++){
        
        if(buffer[i] != ' '){
            nr = nr*10 + buffer[i] - '0';
        }
        else{
            if(!nr_l){
                nr_l = nr;
                nr=0;
            }
            while(buffer[i] == ' ')i++;
            i--;
        }
        //cout << buffer[i] << ' ' << nr << '\n';
    }
    nr_c = nr;
    cout << nr_l <<' ' << nr_c << '\n';
    fflush(0);
    g << nr_l << ' ' << nr_c << '\n';
    for(int i=0; i<nr_c; ++i){
        f >> a;
        while(a){
            g << a << ' ';
            f >> a;
        }
        g << 0 << '\n';
    }

    return 0;
}
