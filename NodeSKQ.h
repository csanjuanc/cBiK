/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   NodeSKQ.h
 * Author: carlos
 *
 * Created on January 13, 2018, 3:57 AM
 */

#ifndef NODESKQ_H
#define NODESKQ_H

#include <vector>
#include <list>
#include <iostream>
using namespace std;

//#include <sdsl/bit_vectors.hpp>
//using namespace sdsl;

class NodeSKQ {
    
//private:
    
    
public:
    double* coordenadas;
    vector<int> palabras;
    
    NodeSKQ();
    NodeSKQ(double* c);
    NodeSKQ(double* c, vector<int> &p);

};

#endif /* NODESKQ_H */

