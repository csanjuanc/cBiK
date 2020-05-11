/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Build_cBiK.h
 * Author: carlos
 *
 * Created on January 16, 2018, 3:28 PM
 */

#ifndef BUILD_CBIK_H
#define BUILD_CBIK_H

#include <stdbool.h>
#include <stdlib.h>
#include <vector>
#include <list>
#include <iostream>
#include <queue>
#include "NodeSKQ.h"
#include "MaxHeap.h"
#include "MinHeap.h"
#include <set>
//Para los tiempos
#include <time.h>
#include <sys/time.h>
#include <sdsl/bit_vectors.hpp>
//Para leer archivos
#include <fstream>
//Para la funcion hash
#include <unordered_map>
#include <string>
#include <locale>

using namespace sdsl;
using namespace std;

class Build_cBiK {
public:
    int contFailB = 0;
    int contFailR = 0;
    int contFailRS = 0;
    Build_cBiK(string datasetName);
    Build_cBiK();
    
    void load_cBiK(string datasetName);
    void printMapa();
    double findBKNN(const double* query, const vector<std::string> &queryKey, int k);
    double findRKNN(const double* query, const vector<std::string> &queryKey, int k, double alpha);
    double rangeQuery(const double* query, const vector<std::string> &queryKey);
    void rangeSearchQuery(const double* query1, const double* query2, vector<int> queryKey, int start, int end, int profundidad, vector<int> &result);
    
private:
    //Atributos generales de construccion
    string datasetName;
    //En esta variable se guarda los puntos recibidos
    vector<NodeSKQ> nodosSKQ;
    //En esta variable se guardan las referencias finales del KD-Tree
    vector<NodeSKQ> nodosKdTree;
    long long int numKeys = 0;
    
    bit_vector mapa;
    sd_vector<> sdVectorKeywords;
    sd_vector<> sdVectorResume;
    //vector<long long int> vectorResume;
    set<long long int> vectorResumeFinal;
    set<long long int> vectorKeys;
    
    //variable global para la distancia maxima de un punto
    double flag[2];
    
    //para convertir las palabras
    locale loc;
    
    //Para el rank
    rank_support_v<0> mapRank;
    rank_support_v<1> mapRank_1;
    
    //Atributos del cBiK
    std::unordered_map<std::string,int> hashMapKeys;
    
    //Distancia maxima entre dos puntos
    double dmax = 162000;
    
    void loadDataset();
    void create_cBiK();
    void export_cBiK();
    
    string toLower(string &palabra);
    
    void initializeReference(vector<NodeSKQ>& coordinates, vector<NodeSKQ>& reference);
    double superKeyCompare(const double *a, const double *b, const long p, const long dim);
    void mergeSort(vector<NodeSKQ> &reference, vector<NodeSKQ>& temporary, const long low, const long high, const long p, const long dim);
    long removeDuplicates(vector<NodeSKQ>& reference, const long i, const long dim);
    void buildKdTree(vector< vector<NodeSKQ> >& references, vector<NodeSKQ>& temporary, const long start, const long end, const long dim, const long depth);
    set<int> generateResume(int start, int end);
    
    double searchEuclideanDistance(const double* p1, const double* p2);
    
    bool checkKeywords(vector<int> query, int pos);
    bool checkSummaryLeft(vector<int> query, int pos, int start, int end);
    bool checkSummaryRight(vector<int> query, int pos, int start, int end);
    
    double timeval_diff(struct timeval *a, struct timeval *b);
    
    //Helper Indices cBiK
    long long int getNewSummaryIndex(int i);
    
    void printVector();
    void printTree();
    void printKdTree(int start, int end, int depth);
    void printTuple(double* tuple);
    void printKeys(vector<bit_vector> arrBits);
    void printKey(bit_vector arrBits);
    
    //Funciones para las ponderaciones
    double getSpatialScore(const double* query, int pos);
    double getKeywordsTextualScore(vector<int> queryKey, int pos);
    double getTotalScore(double &spatialScore, double &textualScore, double &alpha);
    double* getSummariesTextualScores(vector<int> queryKey, int pos, int start, int end);
    
    MaxHeap searchBKNN(const double* query, vector<int> queryKey, int start, int end, const int profundidad, MaxHeap &heap, int totalK);
    MinHeap searchRKNN(const double* query, vector<int> queryKey, int start, int end, MinHeap &heap, int totalK, double alpha);
};

#endif /* BUILD_CBIK_H */

