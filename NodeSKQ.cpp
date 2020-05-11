/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   NodeSKQ.cpp
 * Author: carlos
 * 
 * Created on January 13, 2018, 3:57 AM
 */

#include "NodeSKQ.h"
NodeSKQ::NodeSKQ() {
}

NodeSKQ::NodeSKQ(double* c) {
    this->coordenadas = c;
}

NodeSKQ::NodeSKQ(double* c, vector<int> &p) {
    this->coordenadas = c;
    this->palabras = p;
}

