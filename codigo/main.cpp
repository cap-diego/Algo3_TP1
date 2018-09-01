#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <chrono>
#include <tuple>
using namespace std;

void mostrarPartesBinario(vector<vector<int> > &partes);
void mostrarPartes(vector<vector<int> > &v);
void mostrarSumasDePartes(vector<tuple<int,int> > &);
void generarPartes(vector<vector<int> > &v);
void computarSumasDeSubConj(vector<vector<int> > &partes,vector<tuple<int,int> > &res, vector<int> &conj);
void inicializacionDeAuxiliares(vector<vector<int> > &partes,vector<tuple<int,int> > &res, vector<int> &conj);
int buscarValorObjetivo_BF(vector<vector<int> > &partes, vector<tuple<int,int> > &res, vector<int> &conj,int& valorObjetivo, int &cantidadPasos);
int buscarValorObjetivo_BT(vector<vector<int> > &partes,vector<tuple<int,int> > &res, int& valorObjetivo, vector<int> &conj, int &cantidadPasos);


int main() {
    vector<tuple<int,int> > res;//aca voy a guardar la informacion solo necesaria de cada subconjunto, es decir, la suma de sus elemen y la cardinalidad
    vector<vector<int> > partes;//aca voy a obtener mi pre-cjto de partes, es decir, la combinacion binaria de cada subconjunto de conj
    vector<int> conj;//este serian mis elementos;
    int valorASumar=1;//este seria el valor a sumar
    int cantElem=0;//esta es la cantidad de elementos que van a ingresar
    /*ifstream infile("../arch1.txt");
    infile >> cantElem >> valorASumar;
    cout <<"la cant de elem es: "<< cantElem << " y el vo: "<<valorASumar<<endl;
    int elem;
    while(infile >> elem) {
        cout << "elem: "<<elem<<endl;
        conj.push_back(elem);
    }*/

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    inicializacionDeAuxiliares(partes,res,conj);
    generarPartes(partes);

    //finalizo el clock
    end = std::chrono::system_clock::now();
    int elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>
            (end-start).count();
    time_t end_time = std::chrono::system_clock::to_time_t(end);
    cout <<"tardo:" <<elapsed_seconds<<endl;

    int pasosBF=0;
    cout << "CORRO ALGORITMO DE FUERZA BRUTA"<<endl;
    int cardResultado  = buscarValorObjetivo_BF(partes,res,conj,valorASumar,pasosBF);
    cout << ((cardResultado!=-1)? ("El conjunto que suma eso de menos elementos tiene: " ) : ("No hay conjunto:")) <<cardResultado<<endl;
    cout <<"realizo: " <<pasosBF<<" pasos"<<endl;



    cout << "CORRO ALGORITMO DE BACKTRACKING"<<endl;
    int pasosBT=0;
    int cardResultado2  = buscarValorObjetivo_BT(partes,res,valorASumar,conj,pasosBT);
    cout << ((cardResultado!=-1)? ("El conjunto que suma eso de menos elementos tiene: " ) : ("No hay conjunto:")) <<cardResultado2<<endl;
    cout <<"realizo: " <<pasosBT<<" pasos"<<endl;

    cout << "la proporcion del backtracking con respecto a fuerza bruta es de " << (pasosBF/pasosBT)<<" a 1" <<endl;
    return 0;
}


void generarPartes(vector<vector<int> > &partes) {
    int coefDeCargado=0;//O(1)
    int cargaActual=0, valor=0;
    for (int i = 0; i < partes[0].size(); ++i) { //O(n)
        coefDeCargado *=2; //O(1) * O(n)^= O(n)
        if(coefDeCargado==0) coefDeCargado=1;//O(1) * O(n)^= O(n)
        cargaActual = 0;//O(1) * O(n)= O(n)
        valor=0;//O(1) * O(n)= O(n)
        for (int j = 0; j < partes.size(); j++) {//O(2^n)
            if(cargaActual==coefDeCargado) {
                (valor==0)? valor=1:valor=0;//O(1) * O(2^n) * O(n)
                cargaActual=0;//O(1) * O(2^n) * O(n)
            }
            cargaActual++;//O(1) * O(2^n) * O(n)
            partes[j][i] = valor;//O(1) * O(2^n) * O(n)
        }
    }
}

void computarSumasDeSubConj(vector<vector<int> > &partes,vector<tuple<int,int> > &res, vector<int> &conj) {
    for (int i = 0; i < partes.size(); ++i) { //O(2^n)
        int sumaDeElementos=0, cardinalidad=0;
        for (int j = 0; j < partes[i].size(); j++) {//O(n)
            if(partes[i][j]==1) {
                sumaDeElementos+=conj[j];
                cardinalidad++;
            }
        }
        res[i] = make_tuple(sumaDeElementos,cardinalidad);
    }
}

void mostrarPartesBinario(vector<vector<int> > & partes) {

    for (int k = 0; k < partes.size(); ++k) {
        for (int i = 0; i < partes[i].size(); ++i) {
            cout << partes[k][i] << " ";
            //cout << "(i:"<<k<<"j:"<<i<<") "<< partes[k][i] << " ";
        }
        cout<<endl;
    }
}

void mostrarPartes(vector<vector<int> > &v){
    cout <<"{";
    for (int i = 0; i < v.size(); ++i) {
        cout <<"{";
        for (int j = 0; j < v[i].size(); ++j) {
            cout << v[i][j];
            if(j!= v[i].size()-1) cout <<",";
        }
        cout<<"}";

        if(i!=v.size()-1) cout <<",";
    }
    cout << "}";
}

void inicializacionDeAuxiliares(vector<vector<int> > &partes,vector<tuple<int,int> > &res, vector<int> &conj) {
    //inicializo conj de entrada
    if(conj.size()==0) {
        for (int k = 1; k <= 20; ++k) {
            conj.push_back(k);
        }
    }
    int n = conj.size(); //O(1)
    partes.resize(pow(2,n)); //O(n)
    res.resize(partes.size());//O(n)

    //
    for (int i = 0; i < partes.size(); ++i) {//O(2^n)
        partes[i].resize(n);//O(n)
    }
}

void mostrarSumasDePartes(vector<tuple<int,int> > &res) {
    cout << "{";
    for (int i = 0; i < res.size(); ++i) {
        cout << "(" << get<0>(res[i]) << "," << get<1>(res[i])<<"),";
    }
    cout <<"}"<<endl;
}

int buscarValorObjetivo_BF(vector<vector<int> > &partes, vector<tuple<int,int> > &res, vector<int> &conj,int& valorObjetivo, int &cantidadPasos){
    int found = -1;

    for (int i = 0; i < partes.size(); ++i) { //O(2^n)
        int sumaDeElementos=0, cardinalidad=0;
        for (int j = 0; j < partes[i].size(); j++) {//O(n)
            cantidadPasos++;
            if(partes[i][j]==1) {
                sumaDeElementos+=conj[j];
                cardinalidad++;
            }
        }
        res[i] = make_tuple(sumaDeElementos,cardinalidad);
    }

    for (int i = 0; i < res.size(); ++i) {
        cantidadPasos++;
        if(get<0>(res[i]) == valorObjetivo) {
            //2opc. la sol encontrada ahora es mejor opc que la sol actual (y si no hay actual es trivial q es la mejor sol entonces) o no es mejor que la actual
            if(found > get<1>(res[i]) || found==-1) found = get<1>(res[i]);
        }
    }
    return found;
}

int buscarValorObjetivo_BT(vector<vector<int> > &partes,vector<tuple<int,int> > &res, int& valorObjetivo, vector<int> &conj,int &cantidadPasos) {
    //poda por factibilidad: si voy sumando y encuentro que la suma actual sobrepsa al valor objetivo, dejo de computar.
    //poda por opt: 1_si encuentor un 0 en el conjunto dejo de computar. 2_si la cantidad de elemen actual supera a la de la solucion actual, dejo de sumar
    int solucion=-1;

    bool deboSeguir=true;
    for (int i = 0; i < partes.size(); ++i) { //O(2^n)
        int sumaDeElementos=0, cardinalidad=0;
        deboSeguir=true;
        for (int j = 0; j < partes[i].size() && deboSeguir; j++) {//O(n)
            cantidadPasos++;
            if(partes[i][j]==1) {
                sumaDeElementos+=conj[j];
                cardinalidad++;
                if(conj[j] == 0) { //poda por optimalidad
                    deboSeguir=false;
                    cardinalidad=-1; //uso el valor de -1 para reflejar que no es un dato que necesite
                    sumaDeElementos=-1; //idem arriba
                }
                //res[i].push_back(conj[j]); //cambio a res por un vector de tuplas, donde cada tupla tiene suma y cardinalida del subset
            }
            if(solucion != -1 && solucion < cardinalidad) { //poda por optimalidad
                deboSeguir=false;
                cardinalidad=-1; //uso el valor de -1 para reflejar que no es un dato que necesite
                sumaDeElementos=-1; //idem arriba
            }
            if(j==partes[i].size()-1 && sumaDeElementos==valorObjetivo) solucion=cardinalidad;
            if(sumaDeElementos > valorObjetivo) { //poda por factibilidad
                deboSeguir=false;
                cardinalidad=-1; //uso el valor de -1 para reflejar que no es un dato que necesite
                sumaDeElementos=-1; //idem arriba
            }
        }
        res[i] = make_tuple(sumaDeElementos,cardinalidad);
    }


    return solucion;
}