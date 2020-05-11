/* 
 * File:   main.cpp
 * Author: carlos
 *
 * Created on 26 de noviembre de 2017, 0:05
 */

#include <vector>
#include <list>
#include <iostream>
#include "Build_cBiK.h"
#include "NodeSKQ.h"
//Para leer archivos
#include <fstream>
#include <string>

using namespace std;

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

int main(int argc, char **argv) {
    //Parametros generales
    string dataSet = "datasetTweets3M";
    //string dataSet = "datasetPOIs";
    //string dataSet = "dataset15Puntos";
    
    //Parametros de consulta
    //double query[2] = {4.1, 3.2};
    //double query2[2] = {4.1, 3.2};
    //int k = 1;
    int numKeys = 3;
    int nTest = 1000;
    //double alpha = 0.3;
    //string queryDataSet = dataSet+"_query"+to_string(numKeys);
    //string queryDataSet = "datasetTweets0M_query"+to_string(numKeys);
    string queryDataSet = dataSet+"_rQuery_20K_"+to_string(numKeys);
    
    /**************************************************************************/
    /***************************** BUILD cBiK *********************************/
    /**************************************************************************/
    /*memory_monitor::start();
    auto start = timer::now();
    
    Build_cBiK cBiK = Build_cBiK(dataSet);
    
    auto stop = timer::now();
    memory_monitor::stop();
    
    cout << "\n =================== ESTADISTICAS ===================" << endl;
    cout << "construction time in seconds: " << duration_cast<seconds>(stop-start).count() << endl;
    std::cout << "peak usage = " << memory_monitor::peak() / (1024*1024) << " MB" << std::endl;
    
    cBiK.printMapa();*/
    
    /**************************************************************************/
    /****************************** FIND KNN *********************************/
    /**************************************************************************/
    /*Build_cBiK cBiK = Build_cBiK();
    cBiK.load_cBiK(dataSet);
    cBiK.printMapa();
    
    vector<std::string> queryKey;
    int numToken;
    int ciclo=0;
    string palabras = "";
    
    vector<vector<std::string>> arregloPalabras(nTest);
    double arregloPuntos[nTest][2];
    
    ifstream lecturaQuery;
    lecturaQuery.open(queryDataSet,ios::in);
    for(string linea; getline(lecturaQuery, linea); ) {
        stringstream registro(linea);
        if(ciclo < nTest) {
            palabras = "";
            numToken = 0;
            queryKey.clear();
            for (string dato; getline(registro, dato, ' '); ) {
                if(numToken == 0) {
                    //LATITUD
                    arregloPuntos[ciclo][0] = stod(dato);
                }else if(numToken == 1) {
                    //LONGITUD
                    arregloPuntos[ciclo][1] = stod(dato);
                }else if(numToken > 1){
                    //PALABRAS CLAVES
                    palabras = palabras + dato + " - ";
                    queryKey.push_back(dato);
                    //cout << dato << ", ";
                }
                numToken++;
            }
            
            //Se cargan los datos
            arregloPalabras[ciclo] = queryKey;
        }
        ciclo++;
    }
    lecturaQuery.close();
    
    //Ejecuta las pruebas
    auto start = timer::now();
    for(int i=0; i<arregloPalabras.size(); i++) {
        //cout << (ciclo+1) << "/" << nTest << " >> " << palabras << "(" << query[0] << " , " << query[1] << ") ===>> ";
        //cBiK.findRKNN(arregloPuntos[i], arregloPalabras[i], k, alpha);
        cBiK.findBKNN(arregloPuntos[i], arregloPalabras[i], k);
    }
    auto stop = timer::now();
    
    cout << "******************************* TEST *******************************" << endl;
    cout << "Tiempo BKNN: " << (duration_cast<milliseconds>(stop-start).count())/1000.0 << " == Errores: " << cBiK.contFailB << endl;
    cout << "********************************************************************" << endl;
    */
    /*cout << "******************************* TEST *******************************" << endl;
    cout << "Tiempo RKNN: " << (duration_cast<milliseconds>(stop-start).count())/1000.0 << " == Errores: " << cBiK.contFailR << endl;
    cout << "********************************************************************" << endl;*/
    
    
    /**************************************************************************/
    /*************************** RANGE SEARCH *********************************/
    /**************************************************************************/
    Build_cBiK cBiK = Build_cBiK();
    cBiK.load_cBiK(dataSet);
    cBiK.printMapa();
    
    double timeRangeSearch = 0.0;
    vector<std::string> queryKey;
    int numToken;
    int ciclo=0;
    string palabras = "";
    
    vector<vector<std::string>> arregloPalabras(nTest);
    double arregloPuntos[nTest][4];
    
    ifstream lecturaQuery;
    lecturaQuery.open(queryDataSet,ios::in);
    for(string linea; getline(lecturaQuery, linea); ) {
        stringstream registro(linea);
        if(ciclo < nTest) {
            palabras = "";
            numToken = 0;
            queryKey.clear();
            for (string dato; getline(registro, dato, ' '); ) {
                if(numToken == 0) {
                    //LATITUD
                    //query[0] = stod(dato);
                    arregloPuntos[ciclo][0] = stod(dato);
                }else if(numToken == 1) {
                    //LONGITUD
                    //query[1] = stod(dato);
                    arregloPuntos[ciclo][1] = stod(dato);
                }else if(numToken == 2) {
                    //LONGITUD
                    //query2[0] = stod(dato);
                    arregloPuntos[ciclo][2] = stod(dato);
                }else if(numToken == 3) {
                    //LONGITUD
                    //query2[1] = stod(dato);
                    arregloPuntos[ciclo][3] = stod(dato);
                }else if(numToken > 3){
                    //PALABRAS CLAVES
                    palabras = palabras + dato + " - ";
                    queryKey.push_back(dato);
                    //cout << dato << ", ";
                }
                numToken++;
            }
            
            //Se cargan los datos
            arregloPalabras[ciclo] = queryKey;
        }
        ciclo++;
    }
    lecturaQuery.close();
    
    //Ejecuta las pruebas
    auto start = timer::now();
    for(int i=0; i<arregloPalabras.size(); i++) {
        cBiK.rangeQuery(arregloPuntos[i], arregloPalabras[i]);
    }
    auto stop = timer::now();
    
    cout << "******************************* TEST *******************************" << endl;
    cout << "Tiempo RSKQ: " << (duration_cast<milliseconds>(stop-start).count())/1000.0 << " == Errores: " << cBiK.contFailRS << endl;
    cout << "********************************************************************" << endl;
    
    
    //Se libera la memoria de las coordenadas
    /*for (int i = 0; i < cantidadElementos; i++) {
        delete [] coordinates[i];
    }
    delete [] coordinates;*/
    
    return 0;
}
