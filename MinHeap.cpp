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

#include "MinHeap.h"
#include <iostream>
#include <limits>

using namespace std;

/*MinHeap::MinHeap(const MinHeap& orig) {
}

//Destructor
MinHeap::~MinHeap() {
}*/

//Constructor con el tamanio del Heap
MinHeap::MinHeap(int maxsize) {
    heap = new double *[maxsize];

    for (int i = 0; i < maxsize; i++) {
        heap[i] = new double[2];
        heap[i][0] = 0;//indice del punto mas cercano
        heap[i][1] = numeric_limits<double>::min();//ponderacion del punto
    }

    maxElementos = maxsize;
}


void MinHeap::Intercambia(int i, int j) {
    double temp[2];
    temp[0] = heap[i][0];
    temp[1] = heap[i][1];

    heap[i][0] = heap[j][0];
    heap[i][1] = heap[j][1];
    heap[j][0] = temp[0];
    heap[j][1] = temp[1];
}


void MinHeap::minHeapify(int i, int j) {
    if ((2 * i + 1) <= j) {
        int k;

        if ((2 * i + 2) <= j) {
            k = (heap[2 * i + 2][1] <= heap[2 * i + 1][1]) ? 2 * i + 2 : 2 * i + 1;
        } else {
            k = 2 * i + 1;
        }

        //Ordena en base a la distancia del punto
        if (heap[i][1] > heap[k][1]) {
            Intercambia(i, k);
            minHeapify(k, j);
        }
    }
}


void MinHeap::insert(long posElement, float distanceToQuery) {
    //verifica que llene el heap
    /*if(numElementos < maxElementos) {
        heap[numElementos][0] = posElement;
        heap[numElementos][1] = distanceToQuery;
        numElementos++;
        buildHeap();
    }else {
        //Se asume que el punto que viene es menor al MinHeap
        heap[0][0] = posElement;
        heap[0][1] = distanceToQuery;
        minHeapify(0, maxElementos - 1);
    }*/
    if(numElementos < maxElementos) {
        numElementos++;
    }
    heap[0][0] = posElement;
    heap[0][1] = distanceToQuery;
    minHeapify(0, maxElementos - 1);
}


void MinHeap::buildHeap() {
    //for (int i = (sizeof (heap) / sizeof (heap[0])) / 2; i >= 0; i--) {
    //cout << "Max Elementos: " << maxElementos << endl;
    for (int i = (maxElementos * 0.5); i >= 0; i--) {
        //cout << "I: " << i << endl;
        minHeapify(i, maxElementos - 1);
    }
}

int MinHeap::getMinPoint() {
    return heap[0][0];
}

double MinHeap::getMinScore() {
    return heap[0][1];
}


double** MinHeap::getArray() {
    return heap;
}


int MinHeap::getCantidad() {
    return numElementos;
}

void MinHeap::showArray() {
    for(int i=0; i<maxElementos; i++) {
        cout << heap[i][0] << "|";
    }
    cout << endl;
}

void MinHeap::imprimirDistancias() {
    for(int i=0; i<maxElementos; i++) {
        cout << heap[i][1] << "|";
    }
    cout << endl;
}