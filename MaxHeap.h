/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MaxHeap.h
 * Author: carloco
 *
 * Created on December 20, 2017, 2:41 PM
 */

#ifndef MAXHEAP_H
#define MAXHEAP_H

class MaxHeap {
    
//private:
    double** heap;
    int numElementos = 0;
    int maxElementos = 0;
    
    void Intercambia(int i, int j);
    void maxHeapify(int i, int j);
    
public:
    //MaxHeap(const MaxHeap& orig);
    //virtual ~MaxHeap();
    MaxHeap(int maxsize);
    void insert(long posElement, float distanceToQuery);
    void buildHeap();
    double getDistanceMax() ;
    double** getArray();
    int getCantidad();
    void showArray();
    void imprimirDistancias();

};

#endif /* MAXHEAP_H */

