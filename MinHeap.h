/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MinHeap.h
 * Author: carloco
 *
 * Created on December 20, 2017, 2:41 PM
 */

#ifndef MINHEAP_H
#define MINHEAP_H

class MinHeap {
    
//private:
    double** heap;
    int numElementos = 0;
    int maxElementos = 0;
    
    void Intercambia(int i, int j);
    void minHeapify(int i, int j);
    
public:
    //MinHeap(const MinHeap& orig);
    //virtual ~MinHeap();
    MinHeap(int maxsize);
    void insert(long posElement, float distanceToQuery);
    void buildHeap();
    int getMinPoint();
    double getMinScore() ;
    double** getArray();
    int getCantidad();
    void showArray();
    void imprimirDistancias();

};

#endif /* MINHEAP_H */

