#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <chrono>
#include <tuple>
#include <bitset>
#include <limits>
#include <string>
#define VALOROBJETIVO 2

using namespace std;

typedef std::vector<vector<bool> > matriz;
void mostrarPartesBinario(matriz &partes);
void mostrarPartes(matriz &v);
void generatePrintBinary(int n);
void mostrarSumasDePartes(vector<tuple<int,int> > &);
void computarSumasDeSubConj(matriz &partes, vector<int> &conj);
void generarPartes(matriz &v, int &pasos);
void inicializacionDeAuxiliares(matriz &partes, vector<int> &conj);
int buscarValorObjetivo_BF( vector<int> &conj,int& valorObjetivo, long&,vector<pair<int,int> > &partesAct, matriz &partes);
int buscarValorObjetivo_BT(matriz &partes,int& valorObjetivo, vector<int> &conj, long &cantidadPasos);
int buscarValorObjetivoPD_bottom_up(int &valorObjetivo, vector<int> &conj, long &cantidadPasos);
int buscarValorObjetivoPD_top_down(vector<int>, int, long &, int);
void actualizarPartes(matriz &partes, vector<pair<int,int> > &,vector<int> &conj);
void cargarDatosEntrada(vector<int> &conj, int& valorObjetivo, int& cantidadElementos);
int aux_top_down(vector<int> conjActual, int valorObjetivo, long &cantPasos, int tamInicial, vector<vector<int> > & memoria, int desde);

void powerSet(string set, vector<string >& powerset);
void subsetSums(vector<int>, int n);

int bruteforce(vector<int>& c, int valorObjetivo,long &pasos);
int backtracking(vector<int> &c, int valorObjetivo, long &pasos);
void auxBacktracking(vector<int>& c, vector<int> partes_i, int desde, int valorObjetivo, int& sol, long &pasos);
int sumaConj(vector<int>& c, int valorObjetivo);
vector<int> agregar(vector<int> a, int b);

int main(int argc, char* argv[]) {
    using namespace chrono;
    int valorObjetivo = VALOROBJETIVO;
    vector<int> arr;
    int tam = 20;
    arr.resize(tam);
    for (int i = 0; i < tam; i++) {
        if(i!= 1) arr[i]=3+1;
    }
    arr[1]=5;
    int sol;
    long pasosBF=0;


    auto start = std::chrono::steady_clock::now();
    sol = bruteforce(arr,valorObjetivo,pasosBF);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> diff = duration_cast<duration<double, std::milli> >(end - start);
    cout <<"sol = " <<sol<<" y TARDO: (ms,pasos)"<< diff.count()<<","<<pasosBF<<endl;

    long pasosBT=0;
    start = std::chrono::steady_clock::now();
    sol = backtracking(arr,valorObjetivo,pasosBT);
    end = std::chrono::steady_clock::now();
    diff = duration_cast<duration<double, std::milli> >(end - start);
    cout <<"sol = " <<sol<<" y TARDO: (ms,pasos)"<< diff.count()<<","<<pasosBT<<endl;


    matriz partes;//aca voy a obtener mi pre-cjto de partes, es decir, la combinacion binaria de cada subconjunto de conj
    vector<int> conj;//este serian mis elementos;
    vector<pair<int,int> > partesAct;

    int valorASumar=VALOROBJETIVO;//este seria el valor a sumar
    int cantElem=0;//esta es la cantidad de elementos que van a ingresar
    float tiempoEjecucionTotal=0;


    conj = arr;

    //inicializacionDeAuxiliares(partes,conj);
    //generarPartes(partes, cPasos);
    //actualizarPartes(partes,partesAct,conj);
/*
    cout <<endl<< "-------------------------------------------------------------------------------"<<endl;
    pasosBF=0;
    cout << "CORRO ALGORITMO DE FUERZA BRUTA"<<endl;
    start = std::chrono::steady_clock::now();
    int cardResultado  = buscarValorObjetivo_BF(conj,valorASumar,pasosBF,partesAct,partes);
    end = std::chrono::steady_clock::now();
    diff = duration_cast<duration<double, std::milli> >(end - start);
    cout << ((cardResultado!=-1)? ("El conjunto que suma eso de menos elementos tiene: " ) : ("No hay conjunto:")) <<cardResultado<<endl;
    cout <<"realizo: " <<pasosBF<<" pasos"<<endl;
    cout << "EL ALGORITMO DE FUERZA BRUTA TARDO: " << diff.count() << " ms"<< endl;*/


    cout<<endl << "-------------------------------------------------------------------------------"<<endl;
    cout << "CORRO ALGORITMO DED PROGRAMACION DINAMICA, BOTTOM UP"<<endl;
    long pasosPD=0;
    start = std::chrono::steady_clock::now();
    int cardResultado3 = buscarValorObjetivoPD_bottom_up(valorASumar, conj, pasosPD);
    end = std::chrono::steady_clock::now();
    diff = duration_cast<duration<double, std::milli> >(end - start);
    cout << ((cardResultado3!=-1)? ("El conjunto que suma eso de menos elementos tiene: " ) : ("No hay conjunto:")) <<cardResultado3<<endl;
    cout << "realizo: "<< pasosPD << " pasos"<<endl;
    tiempoEjecucionTotal+=diff.count();
    cout << "EL ALGORITMO DE PD TARDO: " << diff.count() << " ms"<<endl;
    cout << "la proporcion de PD con respecto a fuerza bruta es de " << ceil(float(pasosBF)/float(pasosPD))<<" a 1" <<endl;
    cout << "la proporcion de PD con respecto a backtracking es de " << ceil(float(pasosBT)/float(pasosPD))<<" a 1" <<endl;

    cout<<endl << "-------------------------------------------------------------------------------"<<endl;
    cout << "CORRO ALGORITMO DED PROGRAMACION DINAMICA, TOP DOWN"<<endl;
    long pasosPD_down=0;
    start = std::chrono::steady_clock::now();
    int cardResultado3_down = buscarValorObjetivoPD_top_down(conj,valorASumar, pasosPD_down, conj.size());
    end = std::chrono::steady_clock::now();
    diff = duration_cast<duration<double, std::milli> >(end - start);
    cout << ((cardResultado3_down!=-1)? ("El conjunto que suma eso de menos elementos tiene: " ) : ("No hay conjunto:")) <<cardResultado3_down<<endl;
    cout << "realizo: "<< pasosPD_down << " pasos"<<endl;
    tiempoEjecucionTotal+=diff.count();
    cout << "EL ALGORITMO DE PD TARDO: " << diff.count() << " ms"<<endl;
    cout << "la proporcion de PD con respecto a fuerza bruta es de " << ceil(float(pasosBF)/float(pasosPD_down))<<" a 1" <<endl;
    cout << "la proporcion de PD 0con respecto a backtracking es de " << ceil(float(pasosBT)/float(pasosPD_down))<<" a 1" <<endl;
    cout << "la proporcion de PD_down con respecto a PD es de " << ceil(float(pasosPD)/float(pasosPD_down))<<" a 1" <<endl;

    /*ofstream outputFile;
    outputFile.open("../output.csv");
    if(outputFile.is_open()) {
        //escri
        outputFile << "PD\tBT\tBf"<<endl;
        outputFile<<pasosPD<<"\t"<<pasosBT<<"\t"<<pasosBF<<endl;
        outputFile.close();
    }else cout << "error"<<endl;
    */


    return 0;
}

void actualizarPartes(matriz &partes, vector<pair<int,int> > &partesAct,vector<int> &conj) {
    //partesAct.resize(partes.size());

    for (int i = 0; i < partes.size(); ++i) {
        int sum=0, card=0;
        for (int j = 0; j < partes[i].size(); ++j) {
            if(partes[i][j]==1) {
                sum+=conj[j];
                card++;
            }

        }
        //pair a = make_pair(sum,card);
        partesAct.push_back(make_pair(sum,card));
    }
}

void generatePrintBinary(int n)
{
    // Create an empty queue of strings
    vector<string> q;

    // Enqueue the first binary number
    q.push_back("1");

    // This loops is like BFS of a tree with 1 as root
    // 0 as left child and 1 as right child and so on
    while (n--)
    {
        // print the front of queue
        string s1 = q.front();
        q.pop_back();
        cout << s1 << "\n";

        string s2 = s1;  // Store s1 before changing it

        // Append "0" to s1 and enqueue it
        q.push_back(s1.append("0"));

        // Append "1" to s2 and enqueue it. Note that s2 contains
        // the previous front
        q.push_back(s2.append("1"));
    }
}

void generarPartes(matriz &partes, int &pasos) {

    int coefDeCargado=0;//O(1)
    int cargaActual=0, valor=0;
    for (int i = 0; i < partes[0].size(); ++i) { //O(n)
        coefDeCargado *=2; //O(1) * O(n)^= O(n)
        if(coefDeCargado==0) coefDeCargado=1;//O(1) * O(n)^= O(n)
        cargaActual = 0;//O(1) * O(n)= O(n)
        valor=0;//O(1) * O(n)= O(n)
        for (int j = 0; j < partes.size(); j++) {//O(2^n)
            pasos++;
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
    //cout <<"{";
    for (int i = 0; i < v.size(); ++i) {
        //cout <<"{";
        for (int j = 0; j < v[i].size(); ++j) {
            cout << v[i][j];
            if(j!= v[i].size()-1) cout <<" ";
        }
        //cout<<"}";
        cout <<endl;
    }
    //cout << "}";
}


void powerSet(string set, vector<string>& powerset) {

    string currentSet;
    int vo=7;
    int valorInicioAscii=48;
    int sol=-1;
    unsigned long power = (1 << set.size()) ; // pow(2, set.size())
    //powerset.resize(power);
    for(unsigned long i = 0; i < power; ++i) {

        unsigned long working = i;

        currentSet.clear();
        int suma=0,card=0;
        for(int j = 0; j < sizeof(unsigned long) * 8; ++j) {

            if (working & 1) {
                suma += int(set[j])-valorInicioAscii;
                card++;
                currentSet.push_back(set[j]);
            }

            working >>= 1;

        }
        powerset.push_back(currentSet);
        if(vo == suma) {
            sol = card;
        }
    }
    cout <<"SOL: "<<sol<<endl;

}


void inicializacionDeAuxiliares(matriz &partes, vector<int> &conj) {
    //inicializo conj de entrada
    if(conj.size()==0) {
        for (int k = 0; k < 20; ++k) {
            conj.push_back(k);
            //if(k % 2 !=0) conj.push_back(k);
        }
    }

    int n = conj.size(); //O(1)
    partes.resize(pow(2,n)); //O(n)


    for (int i = 0; i < partes.size(); ++i) {//O(2^n)
        //cout <<i<<" "<<endl;
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

int buscarValorObjetivo_BF(vector<int> &conj,int& valorObjetivo, long &cantidadPasos,vector<pair<int,int> > &partesAct, matriz & partes){
    int sol = -1;

    for (int i = 0; i < partes.size(); ++i) { //O(2^n)
        int sumaDeElementos=0, cardinalidad=0;
        for (int j = 0; j < partes[i].size(); j++) {//O(n)
            cantidadPasos++;
            if(partes[i][j]==1) {
                sumaDeElementos+=conj[j];
                cardinalidad++;
            }
            if(j == partes[i].size()-1 ){
                if(sol ==-1 && sumaDeElementos == valorObjetivo) sol = cardinalidad;
                if(sol > cardinalidad && sumaDeElementos==valorObjetivo) sol = cardinalidad;
            }
        }

    }

/*
    for (int i = 0; i < partesAct.size(); ++i) {//O(2^N)
        cantidadPasos++;
        if(partesAct[i].first == valorObjetivo) {
            //2opc. la sol encontrada ahora es mejor opc que la sol actual (y si no hay actual es trivial q es la mejor sol entonces) o no es mejor que la actual
            if(sol > partesAct[i].second || sol==-1) sol = partesAct[i].second;
        }
    }
    */
    return (sol > conj.size() || sol == 0)? -1:sol;
}

int buscarValorObjetivo_BT(matriz &partes, int& valorObjetivo, vector<int> &conj, int &cantidadPasos) {
    //poda por factibilidad: si voy sumando y encuentro que la suma actual sobrepsa al valor objetivo, dejo de computar.
    //poda por opt: 1_si encuentor un 0 en el conjunto dejo de computar. 2_si la cantidad de elemen actual supera a la de la solucion actual, dejo de sumar
    int solucion=-1;

    bool deboSeguir=true; //flag para las podas

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
                }

            }
            if(solucion != -1 && solucion < cardinalidad) { //poda por optimalidad
                deboSeguir=false;
            }
            if(j==partes[i].size()-1 && sumaDeElementos==valorObjetivo) {
                solucion=cardinalidad;
                /*cout<<"la sol con card: "<<cardinalidad<< " es la: ";
                for (int k = 0; k < partes[i].size(); ++k) {
                    if(partes[i][k]!=0)cout <<conj[k]<<" ";
                }
                cout <<endl;*/
            }
            if(sumaDeElementos > valorObjetivo) { //poda por factibilidad
                deboSeguir=false;
            }
        }
    }


    return (solucion== 0 || solucion > conj.size())? -1:solucion;
}


void subsetSums(int arr[], int n)
{
    int vo = VALOROBJETIVO;
    // There are totoal 2^n subsets
    long long total = 1<<n;

    // Consider all numbers from 0 to 2^n - 1
    for (long long i=0; i<total; i++)
    {
        //if(i % 1000)  cout<<"faltan: "<< total-i<<endl;
        long long sum = 0;

        // Consider binary reprsentation of
        // current i to decide which elements
        // to pick.
        for (int j=0; j<n; j++)
            if (i & (1<<j))
                sum += arr[j];

        // Print sum of picked elements.
        if(sum == vo) {
            cout << "hay sol"<<endl;
            return;
        }

        ///cout << sum << " ";
    }

    cout <<"sin sol"<<endl;
}

int sumaConj(vector<int>& conj){
    int s = 0;
    for (int i = 0; i < conj.size(); ++i) {
        s+=conj[i];
    }
    return s;
}

void auxBacktracking(vector<int>& c, vector<int> partes_i, int desde, int valorObjetivo, int& sol, long &pasos){
    /*cout <<"Llega el: "<<endl<<"[";
    for (int j = 0; j < partes_i.size(); ++j) {
        cout <<partes_i[j]<<",";
    }
    cout<<"]" <<endl;*/

    if(sol==-1 || (sol != -1 && partes_i.size() < sol)){
        if(partes_i.size() > 0) {
            //cout <<"temp tiene elementos"<<endl;
            if(sumaConj(partes_i) == valorObjetivo) {
                if(sol==-1){
                    sol = partes_i.size();
                    //cout << "Primera sol: "<<sol<<endl;

                }
                else {
                    if(sol > partes_i.size()) {
                        sol = partes_i.size();
                        //cout << " SOL MEJOR: " <<sol<<endl;
                    }

                }

            }
        }
    }
    pasos++;
    long pasosInt= 0;
    //if(partes_i.size() == 0) partes_i.push_back(c[0]); //linea agregada para el mejor caso
    for (int i = desde; i < c.size(); ++i) {
        pasosInt++;
        //if(sumaConj(partes_i) < valorObjetivo) {//cota por fact
        if((sumaConj(partes_i) < valorObjetivo) && (sol==-1 || sol!=-1 && sol < partes_i.size())) {//cota por fact y optº
            partes_i.push_back(c[i]);
            desde++;
            //cout << "meto a partes: "<<c[i]<<endl;
            auxBacktracking(c,partes_i,desde,valorObjetivo,sol,pasos);
            //cout << "saco el : "<<partes_i[partes_i.size()-1]<<endl;
            partes_i.pop_back();
        }

    }
    pasos+=pasosInt;

}

vector<int> agregar(vector<int> a, int b) {
    vector<int> c = a;
    c.push_back(b);
    return c;
}

int backtracking(vector<int>& c, int valorObjetivo,long &pasos) {
    vector<int> partes_i;
    int sol = -1;
    auxBacktracking(c,partes_i,0,valorObjetivo,sol,pasos);
    if(sol==-1) {
        return -1;
    }else{
        return sol;
    }

}

int bruteforce(vector<int>& c, int valorObjetivo,long& pasos){
    // There are totoal 2^n subsets
    int n = c.size();
    long long total = 1<<n;
    int sol=-1;
    // Consider all numbers from 0 to 2^n - 1
    for (long long i=0; i<total; i++)
    {
        int card=0;
        long long sum = 0;
        long long aux =1;
        // Consider binary reprsentation of
        // current i to decide which elements
        // to pick.
        pasos += n;
        for (int j=0; j<n; j++){
            if(j==0) aux = 1;
            if(j!= 0) {
                aux*=2;
            }
            if (i & aux) {
                card++;
                sum += c[j];
            }

        }

        // Print sum of picked elements.
        if(sum == valorObjetivo) {
            if(sol==-1){
                sol = card;
            }else {
                if(sol > card) {
                    sol = card;
                }
            }

        }

        ///cout << sum << " ";
    }
    return sol;
}

int buscarValorObjetivoPD_bottom_up(int &valorObjetivo, vector<int> &conj, long &cantidadPasos) {
   /* int sol = -1;
    cantidadPasos = 1;
    vector<vector<int> > soluciones;
    soluciones.resize(conj.size());

    for (int i = 0; i < conj.size(); ++i) {
        soluciones[i].resize(valorObjetivo);
    }

    vector<int> c = conj;*/
    /*c.push_back(0);
    int aux = c[c.size()-1];
    c[c.size()-1] = c[0];
    c[0] = aux;*/
    /*for (int j = 0; j < valorObjetivo; ++j) {
        soluciones[0][j]=c.size();
    }
    for (int i = 0; i < c.size(); i++) {

        for (int j = 0; j < valorObjetivo; ++j) {
            cantidadPasos++;
            if(c[i] > j+1 ) {//si el elemento actual supera el valor objetivo no puede ser solucion
                soluciones[i][j] = soluciones[i][j];
            }else {
                soluciones[i][j] = min((soluciones[i][j]),1+(soluciones[i][j+1-c[i]]));
            }
        }

    }
    sol = soluciones[c.size()-1][valorObjetivo-1];
    //mostrarPartes(soluciones);
    //cout << "SOL ES: "<< sol << " y conj size es: "<<conj.size()<<endl;
    cout <<  "sol up: "<<sol<<endl;
    return (sol < conj.size() || sol == 0)? sol:-1;*/



    //sol vieja
    int sol = -1;
    cantidadPasos = 1;
    vector<vector<int> > soluciones;
    soluciones.resize(conj.size()+1);
    for (int i = 0; i < conj.size()+1; ++i) {
        soluciones[i].resize(valorObjetivo+1);
    }

    vector<int> c = conj;
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
            cantidadPasos++;
            if(c[i] > j ) {//si el elemento actual supera el valor objetivo no puede ser solucion
                soluciones[i][j] = soluciones[i-1][j];
            }else {
                soluciones[i][j] = min((soluciones[i-1][j]),1+(soluciones[i-1][j-c[i]]));
            }
        }

    }
    sol = soluciones[c.size()-1][valorObjetivo];
    //mostrarPartes(soluciones);
    return (sol <= conj.size())? sol:-1;

}

int buscarValorObjetivoPD_top_down(vector<int> conj, int valorObjetivo, long &cantPasos, int tamInicial) {
    //cout << "entro a aux pd, conj actual tam: "<< conjActual.size()<< " y tam inicialm es: " << tamInicial << " y el valor objetivo actual es: " << valorObjetivo<<endl;
    cantPasos++;
    vector<vector<int> > soluciones;
    soluciones.resize(conj.size());
    for (int i = 0; i < conj.size(); ++i) {
        soluciones[i].resize(valorObjetivo+1);
        for (int j = 0; j < soluciones[i].size(); ++j) {
            soluciones[i][j]=-1;
            if(j==i && i == 0) soluciones[i][j]=0;
        }
    }
    int sol=aux_top_down(conj,valorObjetivo,cantPasos,conj.size(),soluciones, conj.size()-1);
    //cout << "SOL TOP DOWN: " <<sol<<endl;
    return (sol > conj.size() || sol == 0)? -1 : sol;

}

int aux_top_down(vector<int> conjActual, int valorObjetivo, long &cantPasos, int tamInicial, vector<vector<int> > & memoria, int desde) {
    cantPasos++;
    //cout << "MUESTRO MAATRIZ"<<endl;
    //mostrarPartes(memoria);
    if(desde < 0) {
        if(valorObjetivo==0) {
            return 0;
        }else {
            return tamInicial+1;
        }
    }
    if(memoria[desde][valorObjetivo] == -1) { // si no esta calculada
            if(conjActual[desde] > valorObjetivo) {  //si el valor actual es mayor que el valor objetivo no puede ser parte de la sol
                //cout <<"desde > 0 y conj[desde] > vo"<<endl;
                //cout <<"conj act; "<<conjActual[desde]<<endl;
                memoria[desde][valorObjetivo] = aux_top_down(conjActual,valorObjetivo,cantPasos,tamInicial,memoria,desde-1);
            }else{
                //cout <<"desde>0 y conj[desde] <=vo"<<endl;
                memoria[desde][valorObjetivo] = min(1 + aux_top_down(conjActual,valorObjetivo-conjActual[desde],cantPasos,tamInicial,memoria,desde-1), aux_top_down(conjActual,valorObjetivo,cantPasos,tamInicial,memoria,desde-1));
            }

    }
    //cout <<"dsp de aplica la funcion"<<endl;
    //mostrarPartes(memoria);
    return memoria[desde][valorObjetivo];
}

void cargarDatosEntrada(vector<int> &conj, int& valorObjetivo, int& cantidadElementos) {

    ifstream infile;
    infile >> cantidadElementos >> valorObjetivo;
    cout <<"la cant de elem es: "<< cantidadElementos << " y el vo: "<<valorObjetivo<<endl;
    int elem;
    while(infile >> elem) {
        cout << "elem: "<<elem<<endl;
        conj.push_back(elem);
    }
}




