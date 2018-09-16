#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <cmath>
#include <chrono>
using namespace std;

int buscarValorObjetivoPD_bottom_up( long long int &valorObjetivo, vector<int> &conj, int& sol, long long int &pasos);
int buscarValorObjetivoPD_top_down(vector<int>,  long long int, int& sol, long long int &pasos);
int aux_top_down(vector<int> conjActual, long long int valorObjetivo, int tamInicial, vector<vector<int> > & memoria, int desde,long long int &pasos);
int bruteforce(vector<int>& c, long long int valorObjetivo, int& sol,long long int &pasos);
int backtracking(vector<int> c, long long int valorObjetivo, int& sol,long long int &pasos);
void auxBacktracking(vector<int>& c, int desde, long long int valorObjetivo, int& sol, long long int sumaActual, long long int cantElem,long long int &pasos);


#define medir_tiempo(K, CODIGO) \
    [&] () -> double {\
        double tiempo_promedio = 0.0;\
        for (int i = 0; i < K; ++i)\
        {\
            auto tiempo_inicio = chrono::steady_clock::now();\
            CODIGO \
            auto tiempo_fin = chrono::steady_clock::now();\
            tiempo_promedio += chrono::duration<double, milli>(tiempo_fin - tiempo_inicio).count();\
        }\
        return tiempo_promedio / (double)K;\
    }();


int main(int argc, char** argv) {

    using namespace chrono;
    ofstream outFile;
    ifstream inFile;
    vector<int> c;
    long long int valorObjetivo, cantElem;
    int cantRep = 2;
    ////INICIO CARGA DATOS
    string nombreIn, nombreOut;
    if(strcmp(argv[1], "")==0) { //significa que la entrada es por teclado
        cin >> cantElem>>valorObjetivo;
        c.resize(cantElem);
        for (int i = 0; i < cantElem; ++i) {
            cin >> c[i];
        }
        int sbf,sbt,spd1,spd2;
        long long int pbf=0,pbt=0,ppd1=0,ppd2=0;
        auto tbf = medir_tiempo(cantRep,bruteforce(c,valorObjetivo,sbf,pbf););
        auto tbt = medir_tiempo(cantRep, backtracking(c,valorObjetivo,sbt,pbt););
        auto tpd1 = medir_tiempo(cantRep, buscarValorObjetivoPD_top_down(c,valorObjetivo,spd1,ppd1););
        auto tpd2 = medir_tiempo(cantRep, buscarValorObjetivoPD_bottom_up(valorObjetivo,c,spd2,ppd2););
        //cout <<"resultados: "<<sbf<<" "<<sbt<<" "<<spd1<<" "<<spd2<<endl;
        outFile<<cantElem<<","<<tbf<< ","<<tbt<<","<<tpd1<<","<<tpd2<<","<<spd1<<","<<valorObjetivo<<","<<pbf<<","<<pbt<<","<<ppd1<<","<<ppd2<<endl;
    }else {
        inFile.open(argv[1]);
        outFile.open("OUT.csv");//,ios_base::app);
        if(!inFile.is_open()) return 1;
        outFile<<"TAM"<<","<<"T-FUERZA-BRUTA"<<","<<"T-BACKTRACKING"<<","<<"T-TOP-DOWN"<<","<<"T-BOTTOM-UP"<<","<< "SOLUCION"<<","<<"V"<<","<<"Pasos-BF"<<","<<"Pasos-BT"<<","<<"Pasos-TD"<<","<<"Pasos-BU"<<endl;
        while (!inFile.eof()){
            c.clear();
            inFile >> cantElem >> valorObjetivo;
            //cout <<cantElem<<" "<<valorObjetivo<<endl;
            c.resize(cantElem);
            for (int j = 0; j < cantElem; ++j) {
                inFile >> c[j];
            }
            int sbf,sbt,spd1,spd2;
            long long int pbf=0,pbt=0,ppd1=0,ppd2=0;
            auto tbf = medir_tiempo(cantRep,bruteforce(c,valorObjetivo,sbf,pbf););
            auto tbt = medir_tiempo(cantRep, backtracking(c,valorObjetivo,sbt,pbt););
            auto tpd1 = medir_tiempo(cantRep, buscarValorObjetivoPD_top_down(c,valorObjetivo,spd1,ppd1););
            auto tpd2 = medir_tiempo(cantRep, buscarValorObjetivoPD_bottom_up(valorObjetivo,c,spd2,ppd2););
            //cout <<"resultados: "<<sbf<<" "<<sbt<<" "<<spd1<<" "<<spd2<<endl;
            outFile<<cantElem<<","<<tbf<< ","<<tbt<<","<<tpd1<<","<<tpd2<<","<<spd1<<","<<valorObjetivo<<","<<pbf<<","<<pbt<<","<<ppd1<<","<<ppd2<<endl;
            //cout<<spd1<<"\t"<<tbf<<"\t"<<tbt<<"\t"<<tpd1<<"\t"<<tpd2<<endl;
        }

        inFile.close();
        outFile.close();
    }

    

    return 0;
}



void auxBacktracking(vector<int>& c, int desde,  long long int valorObjetivo, int& sol, long long int sumaActual,  long long int cantElem,long long int &pasos){
    if(sol==-1 || (sol != -1 && cantElem < sol)){
        if(cantElem > 0) {
            if(sumaActual == valorObjetivo) {
                if(sol==-1){
                    sol = cantElem;
                }
                else {
                    if(sol > cantElem) {
                        sol = cantElem;
                    }
                }
            }
        }
    }
    for (int i = desde; i < c.size() && sumaActual < valorObjetivo; ++i) {
        pasos++;
        if( sol==-1 || sol!=-1 && sol > cantElem) {//cota por fact y optº

            desde++;
            sumaActual +=c[i];
            cantElem++;
            auxBacktracking(c,desde,valorObjetivo,sol,sumaActual, cantElem,pasos);
            cantElem--;
            sumaActual -=  c[i];
        }
    }
}


int backtracking(vector<int> c,  long long int valorObjetivo, int& sol,long long int &pasos) {
    sol=-1;
    sort(c.begin(),c.end()); //ordenamos de menor a mayor en aprox segun c++ (log |c| * |c|)
    auxBacktracking(c,0,valorObjetivo,sol,0,0,pasos);
    if(sol==-1 || sol >= c.size()) {
        return -1;
    }else{
        return sol;
    }
}

int bruteforce(vector<int>& c,  long long int valorObjetivo,int &sol,long long int &pasos){
    int n = c.size();
    long long int total = pow(2,n);
    sol=-1;
    for (long long int i=0; i< total; i++) {
        long card=0;
        long sum = 0;
        for (int j=0; j<n; j++){
            pasos++;
            if (i & (1 << j)) { //chequeo si el elem jesimo tiene el bit en 1 (es decir, debe estar en el subconjunto i), si es 0 no esta entonces lo salto.
                card++;
                sum += c[j];
            }

        }
        if(sum == valorObjetivo) {
            if(sol==-1){
                sol = card;
            }else {
                if(sol > card) {
                    sol = card;
                }
            }
        }
    }
    (sol == -1 || sol >= n)? sol = -1:1;
    return sol;
}

int buscarValorObjetivoPD_bottom_up( long long int &valorObjetivo, vector<int> &conj, int& sol,long long int &pasos) {
    sol = -1;
    if(conj.size() == 0){
        if(valorObjetivo == 0) {
            sol = 0;
            return 0;
        }else {
            sol = -1;
            return -1;
        }

    }else if(valorObjetivo == 0) {
        sol = 0;
        return 0;
    }
    vector<vector<int> > soluciones;
    soluciones.resize(conj.size()+1);
    for (int i = 0; i < conj.size()+1; ++i) {
        soluciones[i].resize(valorObjetivo+1);
    }
    vector<int> c = conj; //agrego un elemento a C para poder indexar de forma mas clara
    c.push_back(0);
    int aux = c[c.size()-1];
    c[c.size()-1] = c[0];
    c[0] = aux;
    for (int j = 0; j < valorObjetivo+1; ++j) {
        soluciones[0][j]=c.size();
    }
    soluciones[0][0]=0;
    for (int i = 1; i < c.size(); i++) {

        for (int j = 0; j < valorObjetivo+1; ++j) {
            pasos++;
            if(c[i] > j ) {//si el elemento actual supera el valor objetivo no puede ser solucion
                soluciones[i][j] = soluciones[i-1][j];
            }else {
                soluciones[i][j] = min((soluciones[i-1][j]),1+(soluciones[i-1][j-c[i]]));
            }
        }

    }
    sol = soluciones[c.size()-1][valorObjetivo];
    (sol > conj.size() || sol == 0)? sol=-1 : 1;
    return sol;
}

int buscarValorObjetivoPD_top_down(vector<int> conj,  long long int valorObjetivo, int& sol,long long int &pasos) {
    sol=-1;
    if(conj.size() == 0){
        if(valorObjetivo == 0) {
            sol=0;
            return 0;
        }else{
            sol=-1;
            return -1;
        }
    }else if(valorObjetivo == 0) {
        sol=0;
        return 0;
    }
    vector<vector<int> > soluciones(conj.size(),vector<int>((valorObjetivo+1),-1)); //inicializo matriz n x V con -1
    soluciones[0][0]=0;
    sol=aux_top_down(conj,valorObjetivo,conj.size(),soluciones, conj.size()-1,pasos);
    (sol > conj.size() || sol == 0)? sol=-1 : 1;
    return sol;
}

int aux_top_down(vector<int> conjActual,  long long int valorObjetivo, int tamInicial, vector<vector<int> > & memoria, int desde,long long int &pasos) {
    pasos++;
    if(desde < 0) {
        if(valorObjetivo==0) {
            return 0;
        }else {
            return tamInicial+1;
        }
    }
    if(memoria[desde][valorObjetivo] == -1) { // si no esta calculada
        if(conjActual[desde] == valorObjetivo){
            memoria[desde][valorObjetivo] = 1;
        }else if(conjActual[desde] > valorObjetivo) {  //si el valor actual es mayor que el valor objetivo no puede ser parte de la sol
            memoria[desde][valorObjetivo] = aux_top_down(conjActual,valorObjetivo,tamInicial,memoria,desde-1,pasos);
        }else{
            memoria[desde][valorObjetivo] = min(1 + aux_top_down(conjActual,valorObjetivo-conjActual[desde],tamInicial,memoria,desde-1,pasos), aux_top_down(conjActual,valorObjetivo,tamInicial,memoria,desde-1,pasos));
        }
    }
    return memoria[desde][valorObjetivo];
}
