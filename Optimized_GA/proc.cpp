#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;
ifstream f("C:/Users/Rares/Documents/Facultate/An2/Sem1/AG/H_2p/res_for_proc.txt");
ofstream g("C:/Users/Rares/Documents/Facultate/An2/Sem1/AG/H_2p/fin_res.txt", std::ios_base::app);

int main(){

    double bestVal=999999;
    double worstVal=-bestVal;
    double avgVal = 0;

    double bestT=999999;
    double worstT=0;
    double avgT = 0;

    double val, time;
    int cnt = 0;
    while(f>>val>>time){
        cnt ++;

        if(val < bestVal)bestVal = val;
        if(val > worstVal)worstVal = val;
        avgVal += val;

        if(time < bestT)bestT = time;
        if(time > worstT)worstT = time;
        avgT += time;
    }

    avgT /= cnt;
    avgVal /= cnt;
    g << fixed << setprecision(5);
    g << "BV: " << bestVal << '\n';
    g << "AV: " << avgVal << '\n';
    g << "WV: " << worstVal << '\n';
    
    g << "BT: " << bestT << '\n';
    g << "AT: " << avgT << '\n';
    g << "WT: " << worstT << '\n';

    g << '\n';

    return 0;
}