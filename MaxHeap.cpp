/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MaxHeap.cpp
 * Author: carloco
 * 
 * Created on December 20, 2017, 2:41 PM
 */

#include "MaxHeap.h"
#include <iostream>
#include <limits>

using namespace std;

/*MaxHeap::MaxHeap(const MaxHeap& orig) {
}

//Destructor
MaxHeap::~MaxHeap() {
}*/

//Constructor con el tamannio del Heap
MaxHeap::MaxHeap(int maxsize) {
    heap = new double *[maxsize];

    for (int i = 0; i < maxsize; i++) {
        heap[i] = new double[2];
        heap[i][0] = 0;//indice del punto mas cercano
        heap[i][1] = numeric_limits<double>::max();//distancia del punto a la consulta
    }

    maxElementos = maxsize;
}


void MaxHeap::Intercambia(int i, int j) {
    double temp[2];
    temp[0] = heap[i][0];
    temp[1] = heap[i][1];

    heap[i][0] = heap[j][0];
    heap[i][1] = heap[j][1];
    heap[j][0] = temp[0];
    heap[j][1] = temp[1];
}


void MaxHeap::maxHeapify(int i, int j) {
    if ((2 * i + 1) <= j) {
        int k;

        if ((2 * i + 2) <= j) {
            k = (heap[2 * i + 2][1] >= heap[2 * i + 1][1]) ? 2 * i + 2 : 2 * i + 1;
        } else {
            k = 2 * i + 1;
        }

        //Ordena en base a la distancia del punto
        if (heap[i][1] < heap[k][1]) {
            Intercambia(i, k);
            maxHeapify(k, j);
        }
    }
}


void MaxHeap::insert(long posElement, float distanceToQuery) {
    //verifica que llene el heap
    /*if(numElementos < maxElementos) {
        heap[numElementos][0] = posElement;
        heap[numElementos][1] = distanceToQuery;
        numElementos++;
        buildHeap();
    }else {
        //Se asume que el punto que viene es mayor al maxheap
        heap[0][0] = posElement;
        heap[0][1] = distanceToQuery;
        maxHeapify(0, maxElementos - 1);
    }*/
    if(numElementos < maxElementos) {
        numElementos++;
    }
    heap[0][0] = posElement;
    heap[0][1] = distanceToQuery;
    maxHeapify(0, maxElementos - 1);
}


void MaxHeap::buildHeap() {
    //for (int i = (sizeof (heap) / sizeof (heap[0])) / 2; i >= 0; i--) {
    //cout << "Max Elementos: " << maxElementos << endl;
    for (int i = (maxElementos * 0.5); i >= 0; i--) {
        //cout << "I: " << i << endl;
        maxHeapify(i, maxElementos - 1);
    }
}


double MaxHeap::getDistanceMax() {
    return heap[0][1];
}


double** MaxHeap::getArray() {
    return heap;
}


int MaxHeap::getCantidad() {
    return numElementos;
}

void MaxHeap::showArray() {
    for(int i=0; i<maxElementos; i++) {
        cout << heap[i][0] << "|";
    }
    cout << endl;
}

void MaxHeap::imprimirDistancias() {
    for(int i=0; i<maxElementos; i++) {
        cout << heap[i][1] << "|";
    }
    cout << endl;
}