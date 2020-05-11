/* 
 * File:   Build_cBiK.cpp
 * Author: carlos
 * 
 * Created on January 16, 2018, 3:28 PM
 */

#include "Build_cBiK.h"

Build_cBiK::Build_cBiK(string datasetName) {
    this->datasetName = datasetName;
    loadDataset();
    create_cBiK();
    export_cBiK();
}

Build_cBiK::Build_cBiK() {
    
}

string Build_cBiK::toLower(string &palabra) {
    string str;
    for (string::size_type i=0; i<palabra.length(); i++){
        str = str + tolower(palabra[i],loc);
    }
    return str;
}

void Build_cBiK::load_cBiK(string datasetName) {
    this->datasetName = datasetName;
    
    /**************************************************************************/
    /********************** CARGAR HASHMAP(KEY,ID) **************************/
    /**************************************************************************/
    cout << "CARGANDO HASHMAP >> " << "cBiK_"+datasetName+"_hashmap" << endl;
    ifstream lecturaHashMap;
    lecturaHashMap.open("cBiK_"+datasetName+"_hashmap",ios::in);
    
    //HASHING PALABRA-ID
    int contTkn=0;
    string keyHash = "";
    int idHash = 0;
    //Se lee cada linea del archivo de texto
    for(string linea; getline(lecturaHashMap, linea); ) {
        stringstream registro(linea);
        
        contTkn = 0;
        keyHash = "";
        idHash = 0;
        //Lee cada elemento de la linea
        for (string dato; getline(registro, dato, ' '); ) {
            if(contTkn == 0) {
                //PALABRA
                keyHash = dato;
            }else if(contTkn == 1) {
                //ID
                idHash = stoi(dato);
            }
            contTkn++;
        }
        hashMapKeys[keyHash] = idHash;
    }
    lecturaHashMap.close();
    
    cout << "TOTAL PALABRAS:" << hashMapKeys.size() << endl;
    this->numKeys = hashMapKeys.size();
    
    /**************************************************************************/
    /************************** CARGAR PUNTOS *******************************/
    /**************************************************************************/
    cout << "CARGANDO PUNTOS >> " << "cBiK_"+datasetName+"_points" << endl;
    long long int totalObjects = 0;
    ifstream lecturaPalabras;
    lecturaPalabras.open("cBiK_"+datasetName+"_points",ios::in);
    
    //Se lee cada linea del archivo de texto
    for(string linea; getline(lecturaPalabras, linea); ) {
        totalObjects++;
    }
    lecturaPalabras.close();
    cout << "Total de objetos >> " << totalObjects << endl;
    
    lecturaPalabras.open("cBiK_"+datasetName+"_points",ios::in);
    
    //arreglo para guardar los puntos
    double **coordinates;
    coordinates = new double *[totalObjects];
    for(int i = 0; i<totalObjects; i++) {
        coordinates[i] = new double[2];//dos dimensiones
    }
    
    //HASHING PALABRA-ID
    int contToken=0;
    int i=0;
    int idKey = 0;
    
    //Se lee cada linea del archivo de texto
    for(string linea; getline(lecturaPalabras, linea); ) {
        stringstream registro(linea);
        
        contToken = 0;
        //Lee cada elemento de la linea
        for (string dato; getline(registro, dato, ' '); ) {
            if(contToken == 0) {
                //LATITUD
                coordinates[i][0] = stod(dato);
            }else if(contToken == 1) {
                //LONGITUD
                coordinates[i][1] = stod(dato);
            }
            contToken++;
        }
        i++;
    }
    lecturaPalabras.close();
    
    nodosKdTree.resize(totalObjects);
    for (int i = 0; i < nodosKdTree.size(); i++) {
        nodosKdTree[i] = NodeSKQ(&(coordinates[i][0]));
    }
    
    /**************************************************************************/
    /********************** CARGAR BITMAP MAPA ******************************/
    /**************************************************************************/
    cout << "CARGANDO BITMAP MAPA >> " << "cBiK_"+datasetName+"_map.sdsl" << endl;
    load_from_file(mapa, "cBiK_"+datasetName+"_map.sdsl");
    cout << "\t NUMERO ELEMENTOS MAPA >> " << mapa.size()*0.5 << endl;
    //Se crea el rank del bitmap del mapa
    mapRank_1 = rank_support_v<1>(&mapa);
    
    /**************************************************************************/
    /********************** CARGAR BITMAP KEYWORDS **************************/
    /**************************************************************************/
    cout << "CARGANDO BITMAP KEYWORDS >> " << "cBiK_"+datasetName+"_keywords.sdsl" << endl;
    load_from_file(sdVectorKeywords, "cBiK_"+datasetName+"_keywords.sdsl");
    cout << "\t NUMERO ELEMENTOS KEYWORDS >> " << sdVectorKeywords.size() << endl;
    
    /**************************************************************************/
    /********************** CARGAR BITMAP SUMMARY ***************************/
    /**************************************************************************/
    cout << "CARGANDO BITMAP SUMMARY >> " << "cBiK_"+datasetName+"_summary.sdsl" << endl;
    load_from_file(sdVectorResume, "cBiK_"+datasetName+"_summary.sdsl");
    cout << "\t NUMERO ELEMENTOS SUMMARY >> " << sdVectorResume.size() << endl;
}

void Build_cBiK::export_cBiK() {
    /**************************************************************************/
    /********************** EXPORTAR HASHMAP(KEY,ID) **************************/
    /**************************************************************************/
    cout << "EXPORTANDO HASHMAP >> " << "cBiK_"+datasetName+"_hashmap" << endl;
    ofstream archivoHM("cBiK_"+datasetName+"_hashmap");
    for (auto& x: hashMapKeys) {
        archivoHM << x.first << " " << x.second << endl;
    }
    archivoHM.close();
    
    /**************************************************************************/
    /************************** EXPORTAR PUNTOS *******************************/
    /**************************************************************************/
    cout << "EXPORTANDO PUNTOS >> " << "cBiK_"+datasetName+"_points" << endl;
    ofstream archivoP("cBiK_"+datasetName+"_points");
    for(int i=0; i<nodosKdTree.size(); i++) {
        archivoP << setprecision(16) << nodosKdTree[i].coordenadas[0] << " " << nodosKdTree[i].coordenadas[1] << endl;
    }
    archivoP.close();
    
    /**************************************************************************/
    /********************** EXPORTAR BITMAP MAPA ******************************/
    /**************************************************************************/
    cout << "EXPORTANDO BITMAP MAPA >> " << "cBiK_"+datasetName+"_map.sdsl" << endl;
    store_to_file(mapa, "cBiK_"+datasetName+"_map.sdsl");
    
    /**************************************************************************/
    /********************** EXPORTAR BITMAP KEYWORDS **************************/
    /**************************************************************************/
    cout << "EXPORTANDO BITMAP KEYWORDS >> " << "cBiK_"+datasetName+"_keywords.sdsl" << endl;
    store_to_file(sdVectorKeywords, "cBiK_"+datasetName+"_keywords.sdsl");
    
    /**************************************************************************/
    /********************** EXPORTAR BITMAP SUMMARY ***************************/
    /**************************************************************************/
    cout << "EXPORTANDO BITMAP SUMMARY >> " << "cBiK_"+datasetName+"_summary.sdsl" << endl;
    store_to_file(sdVectorResume, "cBiK_"+datasetName+"_summary.sdsl");
}

void Build_cBiK::loadDataset() {
    long long int totalObjects = 0;
    
    cout << "LEYENDO DATASET >> " << datasetName << endl;
    
    /**************************************************************************/
    /*************************** CUENTA OBJETOS *******************************/
    /**************************************************************************/
    ifstream lecturaPalabras;
    lecturaPalabras.open(datasetName,ios::in);
    
    //Se lee cada linea del archivo de texto
    for(string linea; getline(lecturaPalabras, linea); ) {
        totalObjects++;
    }
    lecturaPalabras.close();
    cout << "Total de objetos >> " << totalObjects << endl;
    
    /**************************************************************************/
    /*************************** CARGAR OBJETOS *******************************/
    /**************************************************************************/
    lecturaPalabras.open(datasetName,ios::in);
    
    //arreglo para guardar los puntos
    double **coordinates;
    coordinates = new double *[totalObjects];
    for(int i = 0; i<totalObjects; i++) {
        coordinates[i] = new double[2];//dos dimensiones
    }
    
    //arreglo para guardar las palabras claves
    vector< vector<int> > indexKeys;
    
    //HASHING PALABRA-ID
    int contToken=0;
    int i=0;
    int idKey = 0;
    vector<int> arregloClaves;
    std::unordered_map<std::string,int>::const_iterator got;
    
    //Se lee cada linea del archivo de texto
    for(string linea; getline(lecturaPalabras, linea); ) {
        stringstream registro(linea);
        contToken = 0;
        arregloClaves.clear();
        
        //Lee cada elemento de la linea
        for (string dato; getline(registro, dato, ' '); ) {
            if(contToken == 2) {
                //LATITUD
                coordinates[i][0] = stod(dato);
            }else if(contToken == 3) {
                //LONGITUD
                coordinates[i][1] = stod(dato);
            }else if(contToken > 3){
                //PALABRAS CLAVES
                string laPalabra = toLower(dato);
                got = hashMapKeys.find(laPalabra);
                if(got == hashMapKeys.end()) {
                    //No existe la palabra
                    hashMapKeys[laPalabra] = idKey;
                    arregloClaves.push_back(idKey);
                    idKey++;
                }else {
                    //Ya existe, se agrega al arreglo
                    arregloClaves.push_back(got->second);
                }
            }
            contToken++;
        }
        indexKeys.push_back(arregloClaves);
        i++;
    }
    lecturaPalabras.close();
    
    /**************************************************************************/
    /*************************** VECTOR nodosSKQ*******************************/
    /**************************************************************************/
    cout << "TOTAL PALABRAS:" << hashMapKeys.size() << endl;
    this->numKeys = hashMapKeys.size();
    
    //Se guardan las coordenadas y las palabras en un objeto
    this->nodosSKQ.clear();
    this->nodosSKQ.resize(totalObjects);
    for (int i = 0; i < this->nodosSKQ.size(); i++) {
        this->nodosSKQ[i] = NodeSKQ(&(coordinates[i][0]), (indexKeys[i]));
    }
}

void Build_cBiK::create_cBiK(){
    cout << "\nCREANDO INDICE cBiK >> " << endl;
    
    struct timeval iniTime2, endTime2;
    double secs;
    
    //Marco el tiempo de inicio
    gettimeofday(&iniTime2, NULL);

    vector<NodeSKQ>& coordinates = this->nodosSKQ;
    const long numDimensions = 2;//numero de dimensiones
    
    //Inicializa y ordena la referencia del arreglo
    vector< vector<NodeSKQ> > references(numDimensions, vector<NodeSKQ>(coordinates.size()));
    vector<NodeSKQ> temporary(coordinates.size());

    for (long i = 0; i < references.size(); i++) {
        initializeReference(coordinates, references.at(i));
        mergeSort(references[i], temporary, 0, references[i].size() - 1, i, numDimensions);
    }
    
    //Elimine las referencias a coordenadas duplicadas mediante una pasada a través de cada matriz de referencia.
    vector<long> end(references.size());
    for (long i = 0; i < end.size(); i++) {
        end[i] = removeDuplicates(references[i], i, numDimensions);
    }

    //Verifique que se elimino la misma cantidad de referencias de cada matriz de referencia.
    for (long i = 0; i < end.size() - 1; i++) {
        for (long j = i + 1; j < end.size(); j++) {
            if (end.at(i) != end.at(j)) {
                cout << "reference removal error" << endl;
                exit(1);
            }
        }
    }
    
    //Se crea el arreglo final
    long long int totalSize = (end.at(0) + 1);
    nodosKdTree.resize(totalSize);
    
    //Se crea el mapa del doble de la cantidad de elementos
    //La mitad es para el mapeo de nodos internos y hojas, la otra mitad es para indicar los nodos con resumenes
    mapa = bit_vector(totalSize*2,0);
    
    cout << "\t...CONSTRUYENDO iKD-TREE >> " << endl;
    //Crea el KD-Tree por Medianas y marca el mapa
    buildKdTree(references, temporary, 0, end.at(0), numDimensions, 0);
    
    //Crea el constructor del arreglo compacto de palabras claves
    sd_vector_builder localCompactKeywords((totalSize*numKeys), vectorKeys.size());
    
    //Se copian las palabras claves al arreglo compacto
    for (set<long long int>::iterator it=vectorKeys.begin(); it!=vectorKeys.end(); ++it) {
        localCompactKeywords.set(*it);
    }
    
    //Se limpia el contenido de vectorKeys para ahorrar memoria
    vectorKeys.clear();
    
    //Crea el arreglo compacto de palabras claves
    sd_vector<> sdKeywords(localCompactKeywords);
    sdVectorKeywords = sdKeywords;
    
    //Se crea el rank del bitmap del mapa
    mapRank_1 = rank_support_v<1> (&mapa);
    long long int nodosInterno = mapRank_1(mapa.size()*0.5);
    
    cout << "\t...CONSTRUYENDO SUMMARY >> " << endl;
    //Se cuentan los bits resumenes encendido para crear el arreglo
    set<int> totalBitsResumen = generateResume(0, totalSize-1);
    totalBitsResumen.clear();

    cout << "\t...CONSTRUYENDO COMPACT SUMMARY >> " << endl;
    
    //cout << "TOTAL DE UNOS:" << vectorResumeFinal.size() << endl;
    long long int nodosInternosResumen = mapRank_1(mapa.size());
    
    //Crea el arreglo compacto de resumenes
    sd_vector_builder localCompactResume(((nodosInternosResumen-nodosInterno)*2*numKeys), vectorResumeFinal.size());
    
    //Copia los bits del arreglo compacto
    int i = 0;
    for (set<long long int>::iterator it=vectorResumeFinal.begin(); it!=vectorResumeFinal.end(); ++it) {
        localCompactResume.set(*it);  
    }
    
    //Se limpia el contenido de vectorResumeFinal para ahorrar memoria
    vectorResumeFinal.clear();
    
    sd_vector<> sdResume(localCompactResume);
    sdVectorResume = sdResume;

    //marco el tiempo final
    gettimeofday(&endTime2, NULL);
    secs = timeval_diff(&endTime2, &iniTime2);
    
    printf("Tiempo: %.16g microseg -- %.16g miliseg -- %.16g seg\n", secs * 1000000.0, secs * 1000.0, secs);
    cout << endl;
    
    // Verify the k-d tree and report the number of KdNodes.
    /*long numberOfNodes = root->verifyKdTree(numDimensions, 0);
    cout << endl << "Number of nodes = " << numberOfNodes << endl;*/
}

/*
 * Initialize a reference array by creating references into the coordinates array.
 *
 * calling parameters:
 *
 * coordinates - a vector<NodeSKQ> of pointers to each of the (x, y, z, w...) tuples
 * reference - a vector<NodeSKQ> that represents one of the reference arrays
 */
void Build_cBiK::initializeReference(vector<NodeSKQ>& coordinates, vector<NodeSKQ>& reference) {
    for (long j = 0; j < coordinates.size(); j++) {
        reference.at(j) = coordinates.at(j);
    }
}

/*
 * The superKeyCompare method compares two double arrays in all k dimensions,
 * and uses the sorting or partition coordinate as the most significant dimension.
 *
 * calling parameters:
 *
 * a - a double*
 * b - a double*
 * p - the most significant dimension
 * dim - the number of dimensions
 *
 * returns: a double result of comparing two double arrays
 */
double Build_cBiK::superKeyCompare(const double *a, const double *b, const long p, const long dim) {
    double diff = 0;
    for (long i = 0; i < dim; i++) {
        long r = i + p;
        // A fast alternative to the modulus operator for (i + p) < 2 * dim.
        r = (r < dim) ? r : r - dim;
        diff = a[r] - b[r];
        if (diff != 0) {
            break;
        }
    }
    return diff;
}

/*
 * The mergeSort function recursively subdivides the array to be sorted
 * then merges the elements. Adapted from Robert Sedgewick's "Algorithms
 * in C++" p. 166. Addison-Wesley, Reading, MA, 1992.
 *
 * calling parameters:
 *
 * reference - a vector<NodeSKQ> that represents the reference array to sort
 * temporary - a temporary array into which to copy intermediate results;
 *             this array must be as large as the reference array
 * low - the start index of the region of the reference array to sort
 * high - the high index of the region of the reference array to sort
 * p - the sorting partition (x, y, z, w...)
 * dim - the number of dimensions
 * depth - the depth of subdivision
 */
void Build_cBiK::mergeSort(vector<NodeSKQ> &reference, vector<NodeSKQ>& temporary, const long low, const long high, const long p, const long dim) {
    long i, j, k;

    if (high > low) {
        // Evite el desbordamiento al calcular la mediana.
        const long mid = low + ((high - low) >> 1);

        // Recursivamente subdivide las mitades inferior y superior de la matriz.
        mergeSort(reference, temporary, low, mid, p, dim);
        mergeSort(reference, temporary, mid + 1, high, p, dim);

        // Combina los resultados para este nivel de subdivisión.
        for (i = mid + 1; i > low; i--) {
            temporary[i - 1] = reference[i - 1];
        }
        for (j = mid; j < high; j++) {
            temporary[mid + (high - j)] = reference[j + 1]; // Evite el desbordamiento de direcciones.
        }
        for (k = low; k <= high; k++) {
            if(superKeyCompare(temporary[i].coordenadas, temporary[j].coordenadas, p, dim) < 0) {
                reference[k] = temporary[i++];
            }else {
                reference[k] = temporary[j--];
            }
        }
    }
}

/*
 * Check the validity of the merge sort and remove duplicates from a reference array.
 *
 * calling parameters:
 *
 * reference - a vector<NodeSKQ> that represents one of the reference arrays
 * i - the leading dimension for the super key
 * dim - the number of dimensions
 *
 * returns: the end index of the reference array following removal of duplicate elements
 */
long Build_cBiK::removeDuplicates(vector<NodeSKQ>& reference, const long i, const long dim) {
    long end = 0;
    for (long j = 1; j < reference.size(); j++) {
        double compare = superKeyCompare(reference[j].coordenadas, reference[j - 1].coordenadas, i, dim);
        if (compare < 0) {
            cout << "merge sort failure: superKeyCompare(ref[" << j << "], ref["
                    << j - 1 << "], (" << i << ") = " << compare << endl;
            exit(1);
        } else if (compare > 0) {
            reference[++end] = reference[j];
        }
    }
    return end;
}

/*
 * This function builds a k-d tree by recursively partitioning the
 * reference arrays and adding kdNodes to the tree.  These arrays
 * are permuted cyclically for successive levels of the tree in
 * order that sorting occur on x, y, z, w...
 *
 * calling parameters:
 *
 * references - a vector< vector<NodeSKQ> > of pointers to each of the (x, y, z, w...) tuples
 * temporary - a vector<NodeSKQ> that is used as a temporary array
 * start - start element of the reference arrays
 * end - end element of the reference arrays
 * dim - the number of dimensions
 * depth - the depth in the tree
 *
 * returns: a KdNode pointer to the root of the k-d tree
 */
void Build_cBiK::buildKdTree(vector< vector<NodeSKQ> >& references, vector<NodeSKQ>& temporary, const long start, const long end, const long dim, const long depth) {
    // The axis permutes as x, y, z, w... and addresses the referenced data.
    long axis = depth % dim;
    
    if (end == start) {
        //Se agrega la key al punto, siempre es una hoja
        nodosKdTree[end] = references[0][end];
        
        //Se copian las palabras claves segun la posicion del futuro arreglo de bits
        for(int i=0; i<nodosKdTree[end].palabras.size(); i++) {
            vectorKeys.insert((end*numKeys)+nodosKdTree[end].palabras[i]);
        }
        
        //cout << "Start: " << start << endl;
    } else if (end == start + 1) {
        // Two references were passed to this function in sorted order, so store the start
        // element at this level of the tree and store the end element as the > child.
        nodosKdTree[start] = references[0][start];
        nodosKdTree[end] = references[0][end];
        
        //Marco los nodos internos en el mapa
        mapa[start] = 1;
        
        //Se copian las palabras claves segun la posicion del futuro arreglo de bits
        for(int i=0; i<nodosKdTree[start].palabras.size(); i++) {
            vectorKeys.insert((start*numKeys)+nodosKdTree[start].palabras[i]);
        }
        for(int i=0; i<nodosKdTree[end].palabras.size(); i++) {
            vectorKeys.insert((end*numKeys)+nodosKdTree[end].palabras[i]);
        }
        
        //cout << "Start: " << start << " End: " << end << endl;
    } else if (end == start + 2) {
        // Three references were passed to this function in sorted order, so
        // store the median element at this level of the tree, store the start
        // element as the < child and store the end element as the > child.
        nodosKdTree[start + 1] = references[0][start + 1];
        nodosKdTree[start] = references[0][start];
        nodosKdTree[end] = references[0][end];
        
        //Marco los nodos internos en el mapa
        mapa[start + 1] = 1;
        
        //Se copian las palabras claves segun la posicion del futuro arreglo de bits
        for(int i=0; i<nodosKdTree[(start + 1)].palabras.size(); i++) {
            vectorKeys.insert(((start + 1)*numKeys)+nodosKdTree[(start + 1)].palabras[i]);
        }
        for(int i=0; i<nodosKdTree[start].palabras.size(); i++) {
            vectorKeys.insert((start*numKeys)+nodosKdTree[start].palabras[i]);
        }
        for(int i=0; i<nodosKdTree[end].palabras.size(); i++) {
            vectorKeys.insert((end*numKeys)+nodosKdTree[end].palabras[i]);
        }
        
        //cout << "Start: " << start << " Start + 1: " << (start+1) << " End: " << end << endl;
    } else if (end > start + 2) {
        // More than three references were passed to this function, so
        // the median element of references[0] is chosen as the tuple about
        // which the other reference arrays will be partitioned.  Avoid
        // overflow when computing the median.
        const long median = start + ((end - start) / 2);
        
        //Se guardan los puntos de la mediana del arreglo
        nodosKdTree[median] = references[0][median];
        
        //Marco los nodos internos en el mapa
        mapa[median] = 1;
        //Marco los nodos internos que tienen resumenes explicitos sumando la cantidad de elementos
        mapa[median+nodosKdTree.size()] = 1;
        
        //Se copian las palabras claves segun la posicion del futuro arreglo de bits
        for(int i=0; i<nodosKdTree[median].palabras.size(); i++) {
            vectorKeys.insert((median*numKeys)+nodosKdTree[median].palabras[i]);
        }
        
        //cout << "Median: " << median << endl;
        
        // Copy references[0] to the temporary array before partitioning.
        for (long i = start; i <= end; i++) {
            temporary[i] = references[0][i];
        }

        // Process each of the other reference arrays in a priori sorted order
        // and partition it by comparing super keys.  Store the result from
        // references[i] in references[i-1], thus permuting the reference
        // arrays.  Skip the element of references[i] that that references
        // a point that equals the point that is stored in the new k-d node.
        long lower, upper, lowerSave, upperSave;
        for (long i = 1; i < dim; i++) {

            // Process one reference array.  Compare once only.
            lower = start - 1;
            upper = median;
            for (long j = start; j <= end; j++) {
                double compare = superKeyCompare(references[i][j].coordenadas, nodosKdTree[median].coordenadas, axis, dim);
                if (compare < 0) {
                    references[i - 1][++lower] = references[i][j];
                } else if (compare > 0) {
                    references[i - 1][++upper] = references[i][j];
                }
            }

            // Check the new indices for the reference array.
            if (lower < start || lower >= median) {
                cout << "incorrect range for lower at depth = " << depth << " : start = "
                        << start << "  lower = " << lower << "  median = " << median << endl;
                exit(1);
            }

            if (upper <= median || upper > end) {
                cout << "incorrect range for upper at depth = " << depth << " : median = "
                        << median << "  upper = " << upper << "  end = " << end << endl;
                exit(1);
            }

            if (i > 1 && lower != lowerSave) {
                cout << " lower = " << lower << "  !=  lowerSave = " << lowerSave << endl;
                exit(1);
            }

            if (i > 1 && upper != upperSave) {
                cout << " upper = " << upper << "  !=  upperSave = " << upperSave << endl;
                exit(1);
            }

            lowerSave = lower;
            upperSave = upper;
        }

        // Copy the temporary array to references[dim-1] to finish permutation.
        for (long i = start; i <= end; i++) {
            references[dim - 1][i] = temporary[i];
        }

        // Recursively build the < branch of the tree.
        buildKdTree(references, temporary, start, lower, dim, depth + 1);
        
        // Recursively build the > branch of the tree.
        buildKdTree(references, temporary, median + 1, upper, dim, depth + 1);
    } else if (end < start) {

        // This is an illegal condition that should never occur, so test for it last.
        cout << "error has occurred at depth = " << depth << " : end = " << end
                << "  <  start = " << start << endl;
        exit(1);
    }
}

long long int Build_cBiK::getNewSummaryIndex(int i) {
    long long int nodosInternos = mapRank_1(mapa.size()*0.5);
    long long int posResumen = mapRank_1(i+nodosKdTree.size());
    
    return 2*numKeys*(posResumen-nodosInternos);
}

set<int> Build_cBiK::generateResume(int start, int end) {
    set<int> vectorRetorno;
    
    if (end == start) {
        //Se agrega la key al punto, siempre es una hoja
        for(int i=0; i<nodosKdTree[end].palabras.size(); i++) {
            vectorRetorno.insert(nodosKdTree[end].palabras[i]);
        }
    } else if (end == start + 1) {
        //Padre start
        //HD end
        //variable pos tiene la posicion del bit de inicio del correspondiente resumen
        for(int i=0; i<nodosKdTree[end].palabras.size(); i++) {
            vectorRetorno.insert(nodosKdTree[end].palabras[i]);
        }
    } else if (end == start + 2) {
        //Padre start + 1
        //HI start
        //HD end
        for(int i=0; i<nodosKdTree[start].palabras.size(); i++) {
            vectorRetorno.insert(nodosKdTree[start].palabras[i]);
        }
        for(int i=0; i<nodosKdTree[end].palabras.size(); i++) {
            vectorRetorno.insert(nodosKdTree[end].palabras[i]);
        }
    } else if (end > start + 2) {
        // Padre median
        const long median = start + ((end - start) / 2);
        long long int pos = getNewSummaryIndex(median);

        // Recursively build the < branch of the tree.
        set<int> resumenI = generateResume(start, median-1);

        //Copio el hijo izquierdo
        int hijoI = (start+median-1) * 0.5;
        for(int i=0; i<nodosKdTree[hijoI].palabras.size(); i++) {
            resumenI.insert(nodosKdTree[hijoI].palabras[i]);
        }
        
        //Agrego al vector de resumenes las posiciones del resumen izquierdo
        for (set<int>::iterator it=resumenI.begin(); it!=resumenI.end(); ++it){
            vectorRetorno.insert(*it);
            vectorResumeFinal.insert(pos+*it);
        }
        resumenI.clear();
        
        // Recursively build the > branch of the tree.
        set<int> resumenD = generateResume(median + 1, end);
        
        //Copio el hijo derecho
        int hijoD = (median + 1 + end) * 0.5;
        for(int i=0; i<nodosKdTree[hijoD].palabras.size(); i++) {
            resumenD.insert(nodosKdTree[hijoD].palabras[i]);
        }
        
        //Agrego al vector de resumenes las posiciones del resumen derecho
        for (set<int>::iterator it=resumenD.begin(); it!=resumenD.end(); ++it){
            vectorRetorno.insert(*it);
            vectorResumeFinal.insert(pos+numKeys+*it);
        }
        resumenD.clear();
    } else if (end < start) {
        // This is an illegal condition that should never occur, so test for it last.
        cout << "error has occurred at : end = " << end << "  <  start = " << start << endl;
        exit(1);
    }
    
    return vectorRetorno;
}

double Build_cBiK::searchEuclideanDistance(const double* p1, const double* p2) {
    //calcula la distancia para dos dimensiones
    return (p2[0] - p1[0])*(p2[0] - p1[0]) + (p2[1] - p1[1])*(p2[1] - p1[1]);
}

bool Build_cBiK::checkKeywords(vector<int> query, int pos) {
    bool resp = true;
    int i=0;
    long long int posInicial = pos*numKeys;
    
    while(i<query.size() && resp){
        if(sdVectorKeywords[(posInicial+query[i])] == 0) {
            resp = false;
        }
        i++;
    }
        
    return resp;
}

bool Build_cBiK::checkSummaryLeft(vector<int> query, int pos, int start, int end) {
    bool resp = true;
    if(mapa[pos] == 0) {
        //Es una hoja
        resp = false;
    }else {
        //Es un nodo interno, se revisa si es un un nodo interno final
        if (end == start) {
            //Siempre es una hoja, no tiene resumen
            resp = false;
        } else if (end == start + 1) {
            //NI start
            //HD end
            //No tiene resumen izquierdo
            resp = false;
        } else if (end == start + 2) {
            //NI start + 1
            //HI start
            //HD end
            //El resumen izquierdo corresponde a la keyword del HI
            int i=0;
            long long int posInicial = start*numKeys;

            while(i<query.size() && resp){
                if(sdVectorKeywords[(posInicial+query[i])] == 0) {
                    resp = false;
                }
                i++;
            }
        } else if (end > start + 2) {
            //NI pos
            //El resumen izquierdo corresponde al summary de pos
            int i=0;
            long long int posInicial = getNewSummaryIndex(pos);
        
            while(i<query.size() && resp){
                if(sdVectorResume[(posInicial+query[i])] == 0) {
                    resp = false;
                }
                i++;
            }
        } else if (end < start) {
            // This is an illegal condition that should never occur, so test for it last.
            cout << "error has occurred at : end = " << end << "  <  start = " << start << endl;
            exit(1);
        }
        
    }
    
    return resp;
}

bool Build_cBiK::checkSummaryRight(vector<int> query, int pos, int start, int end) {
    bool resp = true;
    
    if(mapa[pos] == 0) {
        //Es una hoja
        resp = false;
    }else {
        //Es un nodo interno, se revisa si es un un nodo interno final
        if (end == start) {
            //Siempre es una hoja, no tiene resumen
            resp = false;
        } else if (end == start + 1) {
            //NI start
            //HD end
            //El resumen derecho corresponde a la keyword del HD
            int i=0;
            long long int posInicial = end*numKeys;

            while(i<query.size() && resp){
                if(sdVectorKeywords[(posInicial+query[i])] == 0) {
                    resp = false;
                }
                i++;
            }
        } else if (end == start + 2) {
            //NI start + 1
            //HI start
            //HD end
            //El resumen derecho corresponde a la keyword del HD
            int i=0;
            long long int posInicial = end*numKeys;

            while(i<query.size() && resp){
                if(sdVectorKeywords[(posInicial+query[i])] == 0) {
                    resp = false;
                }
                i++;
            }
        } else if (end > start + 2) {
            //NI pos
            //El resumen derecho corresponde al summary de pos
            int i=0;
            long long int posInicial = getNewSummaryIndex(pos);
        
            while(i<query.size() && resp){
                if(sdVectorResume[(posInicial+numKeys+query[i])] == 0) {
                    resp = false;
                }
                i++;
            }
        } else if (end < start) {
            // This is an illegal condition that should never occur, so test for it last.
            cout << "error has occurred at : end = " << end << "  <  start = " << start << endl;
            exit(1);
        }
        
    }
    
    return resp;
}

double Build_cBiK::timeval_diff(struct timeval *a, struct timeval *b){
  return
    (double)(a->tv_sec + (double)a->tv_usec/1000000) -
    (double)(b->tv_sec + (double)b->tv_usec/1000000);
}

void Build_cBiK::printMapa() {
    //printTree();
    cout << endl;
    cout << "******************************* MAPA *******************************" << endl;
    cout << "Largo mapa:" << mapa.size()*0.5 << endl;
    cout<< "Mapa size original: " << size_in_bytes(mapa)<< endl;
    //cout << mapa << endl;
    
    cout << "***************************** KEYWORDS *****************************" << endl;
    cout << "Largo:" << sdVectorKeywords.size() << endl;
    cout<< "[keywords size]: " << size_in_bytes(sdVectorKeywords)<< endl;
    /*for(int i=0; i<mapa.size()*0.5; i++) {
        for(int j=0; j<numKeys; j++) {
            cout << sdVectorKeywords[(i*numKeys)+j];
        }
        cout << endl;
    }*/
    //cout << sdVectorKeywords << endl;
    
    cout << endl;
    
    cout << "***************************** RESUMENES ****************************" << endl;
    cout << "Largo:" << sdVectorResume.size() << endl;
    cout<< "[Resume size]: " << size_in_bytes(sdVectorResume)<< endl;
    
    /*long long int nodosInterno = mapRank_1(mapa.size()*0.5);
    long long int nodosInternosResumen = mapRank_1(mapa.size());
    
    for(int i=0; i<(nodosInternosResumen-nodosInterno); i++) {
        for(int j=0; j<numKeys; j++) {
            cout << sdVectorResume[(i*2*numKeys)+j];
        }
        cout << " - ";
        for(int j=0; j<numKeys; j++) {
            cout << sdVectorResume[(i*2*numKeys)+j+numKeys];
        }
        cout << endl;
    }*/
    //cout << sdVectorResume << endl;
    cout << endl;
}

void Build_cBiK::printVector() {
    /*cout << endl;
    cout << "****************************** VECTOR ******************************" << endl;
    for (int i = 0; i < nodosKdTree.size(); i++) {
        cout << i << ".- ";
        printTuple(nodosKdTree[i].coordenadas);
        cout << " - ";
        //Imprime las keys propias del punto
        cout << *nodosKdTree[i].palabras << endl;
    }
    cout << "********************************************************************" << endl;
    cout << endl;*/
}

void Build_cBiK::printTree() {
    cout << endl;
    cout << "************************** ARBOL KD-TREE ***************************" << endl;
    printKdTree(0, (nodosKdTree.size() - 1), 0);
    cout << "********************************************************************" << endl;
    cout << endl;
}

void Build_cBiK::printKdTree(int start, int end, int depth) {
    int mid = (start + end) / 2;

    if (end > mid) {
        printKdTree(mid + 1, end, depth + 1);
    }

    for (int i = 0; i < depth; i++) {
        cout << "         ";
    }
    cout << start << "," << end << " Pos:" << mid << " == ";
    printTuple(nodosKdTree[mid].coordenadas);
    cout << " == ";
    
    int pos = numKeys*mid;
    for(int i=0; i<numKeys; i++) {
        if(sdVectorKeywords[i+pos] == 1) {
            cout << "1";
        }else {
            cout << "0";
        }
    }
    
    /*for(int i=0; i<numKeys ; i++) {
        bool existe = false;
        for(int j=0; j<nodosKdTree[mid].palabras.size(); j++) {
            //cout << "Revisa: " << nodosKdTree[mid].palabras[j] << " - " << i << endl;
            if(nodosKdTree[mid].palabras[j] == i) {
                existe = true;
            }
        }
        if(existe) {
            cout << "1";
        }else {
            cout << "0";
        }
    }*/

    cout << endl;
    if (start < mid) {
        printKdTree(start, mid - 1, depth + 1);
    }
}

void Build_cBiK::printTuple(double* tuple) {
    printf("%.16g, %.16g", tuple[0], tuple[1]); 
    //cout << "(" << tuple[0] << "," << tuple[1] << ")";
}

void Build_cBiK::printKeys(vector<bit_vector> arrBits) {
    for(int i=0; i<arrBits.size(); i++) {
        cout << "(" << arrBits[i] << ")";
    }
}

void Build_cBiK::printKey(bit_vector arrBits) {
    cout << "(" << arrBits << ")";
}

double Build_cBiK::rangeQuery(const double* query, const vector<std::string> &queryKey) {
    //struct timeval iniTime2, endTime2;
    //double secs;
    
    //Atributos basicos para la busqueda
    int start = 0;
    int end = nodosKdTree.size() - 1;
    int depth = 0;
    //Indice de los keys
    vector<int> indexQueryKey;
    string str;
    //Variable donde se guardan los indices a los puntos que estan dentro de la region
    vector<int> result;
    
    //gettimeofday(&iniTime2, NULL);
    //Se buscan los indices de las palabras claves recibidas
    for(int i=0; i<queryKey.size(); i++) {
        //printTuple(query);
        str = queryKey[i];
        str = toLower(str);
        indexQueryKey.push_back(hashMapKeys[str]);
        //cout << str << " >> " << hashMapKeys[str] << endl;
    }
    
    //Se ordenan los puntos del rectangulo
    double minPoint[2];
    double maxPoint[2];
    
    //query[0] = latitud p1
    //query[1] = longitud p1
    //query[2] = latitud p2
    //query[3] = longitud p2
    minPoint[0] = (query[0] < query[2]) ? query[0] : query[2];
    maxPoint[0] = (query[0] > query[2]) ? query[0] : query[2];
    minPoint[1] = (query[1] < query[3]) ? query[1] : query[3];
    maxPoint[1] = (query[1] > query[3]) ? query[1] : query[3];

    //Se realiza la busqueda por rango
    rangeSearchQuery(minPoint, maxPoint, indexQueryKey, start, end, depth, result);
    //gettimeofday(&endTime2, NULL);
    
    /*if(result.empty()) {
        cout << "\nNo hay puntos dentro del rango dado." << endl;
        contFailRS++;
    }else {
        for(int i=0; i<result.size(); i++) {
            printTuple(nodosKdTree[result[i]].coordenadas);
            cout << endl;
        }
    }*/
    
    //cout << endl;
    
    //secs = timeval_diff(&endTime2, &iniTime2);
    
    return 0;//secs * 1000.0;//
}

void Build_cBiK::rangeSearchQuery(const double* query1, const double* query2, vector<int> queryKey, int start, int end, int profundidad, vector<int> &result) {
    //Determina el eje x - y dada la profundidad
    const int p = profundidad % 2;//Debido a las dos dimensiones
    int mid = (start + end) * 0.5;
    
    //Revisa si el punto esta contenido en el rectangulo
    if( (query1[0] <= nodosKdTree[mid].coordenadas[0] && nodosKdTree[mid].coordenadas[0] <= query2[0]) &&
        (query1[1] <= nodosKdTree[mid].coordenadas[1] && nodosKdTree[mid].coordenadas[1] <= query2[1])){
        //Revisa si posee las palabras claves
        if(checkKeywords(queryKey, mid)) {
            result.push_back(mid);
        }
    }
    
    //cout << "Punto visitado: ";
    //printTuple(nodosKdTree[mid].coordenadas);
    //cout << endl;
    
    //Si es distinto es porque es un nodo interno
    if(start != end) {
        //Me aseguro de obtener los valores menor y mayor correctos
        double menor;
        double mayor;
        if(query1[p] > query2[p]) {
            mayor = query1[p];
            menor = query2[p];
        }else {
            mayor = query2[p];
            menor = query1[p];
        }
        
        if(mayor <= nodosKdTree[mid].coordenadas[p]) {
            //Revisa si existe la palabra clave por la izquierda
            if(checkSummaryLeft(queryKey, mid, start, end)) {
                rangeSearchQuery(query1, query2, queryKey, start, mid-1, profundidad + 1, result);
            }
        }else if(menor > nodosKdTree[mid].coordenadas[p]) {
            //Revisa si existe la palabra clave por la derecha
            if(checkSummaryRight(queryKey, mid, start, end)) {
                rangeSearchQuery(query1, query2, queryKey, mid+1, end, profundidad + 1, result);
            }
        }else {
            //Se revisan ambos lados
            if(checkSummaryLeft(queryKey, mid, start, end)) {
                rangeSearchQuery(query1, query2, queryKey, start, mid-1, profundidad + 1, result);
            }
            if(checkSummaryRight(queryKey, mid, start, end)) {
                rangeSearchQuery(query1, query2, queryKey, mid+1, end, profundidad + 1, result);
            }
        }
    }
}

double Build_cBiK::findBKNN(const double* query, const vector<std::string> &queryKey, int k) {
    //cout << "**************************** KNN KDTREE ****************************" << endl;
    struct timeval iniTime2, endTime2;
    double secs;
    
    //Atributos basicos para la busqueda
    int start = 0;
    int end = nodosKdTree.size() - 1;
    int depth = 0;
    string str;
    
    int totalPuntos = nodosKdTree.size();
    
    /*long long int posInicial = 187361*numKeys;
    for(int i=0; i<numKeys; i++) {
        if(sdVectorKeywords[posInicial+i] == 1){
            cout << "indice: " << i << endl;
        }
            
    }*/
    //cout << "ESPECIAL: " << mapa[187361] << endl;

    //Indice de los keys
    vector<int> indexQueryKey;

    if (k < 1) {
        cout << "k debe ser mayor o igual a 1" << endl;
    }else if (k > totalPuntos) {
        cout << "k debe ser menor o igual al total de puntos (" << totalPuntos << ")" << endl;
    }else{
        MaxHeap heap = MaxHeap(k);
        gettimeofday(&iniTime2, NULL);
        for(int i=0; i<queryKey.size(); i++) {
            //printTuple(query);
            str = queryKey[i];
            str = toLower(str);
            indexQueryKey.push_back(hashMapKeys[str]);
            //cout << str << " >> " << hashMapKeys[str] << endl;
        }
        
        heap = searchBKNN(query, indexQueryKey, start, end, depth, heap, k);
        /*gettimeofday(&endTime2, NULL);
        
        secs = timeval_diff(&endTime2, &iniTime2);
        //Muestra el tiempo por cada consulta
        cout << secs * 1000.0 << endl;

        //Se pide el arreglo (HEAP) para mostrar los puntos
        double** elArreglo;
        elArreglo = heap.getArray();

        for (int i = 0; i < k; i++) {
            if(elArreglo[i][1] == numeric_limits<double>::max()) {
                cout << "\tK=" << (i+1) << " => No existe." << endl;
                //contFailB++;
            }else {
                cout << "\tK=" << (i+1);
                cout << " => Punto: ";
                printTuple(nodosKdTree[(int)elArreglo[i][0]].coordenadas);
                cout << " -- ";
                //printKey(*nodosKdTree[(int)elArreglo[i][0]].palabras);
                cout << " Distancia: " << elArreglo[i][1];
                
                cout << endl;
            }
        }*/
    }
    //cout << "********************************************************************" << endl;
    //cout << endl;
    
    return secs * 1000.0;//retorna las milesimas
}

MaxHeap Build_cBiK::searchBKNN(const double* query, vector<int> queryKey, int start, int end, const int profundidad, MaxHeap &heap, int totalK) {
    //Determina el eje x - y dada la profundidad
    const int p = profundidad % 2;//Debido a las dos dimensiones
    int mid = (start + end) * 0.5;
    
    //cout << "Punto visitado: ";
    //printTuple(nodosKdTree[mid].coordenadas);
    //cout << endl;
    
    //Es un nodo interno
    if(start != end) {
        if(query[p] <= nodosKdTree[mid].coordenadas[p]) {
            //Busca por la izquierda si existen las palabras claves
            //cout << "START: " << start << " == MID: " << mid << " == END: " << end << " == CHECK LEFT: " << checkSummaryLeft(queryKey, mid, start, end) << endl;
            if(checkSummaryLeft(queryKey, mid, start, end)) {
                //Actualiza el Heap
                searchBKNN(query, queryKey, start, mid-1, profundidad + 1, heap, totalK);

                //Busca si intersecta un subespacio mediante la circunferencia
                if ((abs(nodosKdTree[mid].coordenadas[p] - query[p]) < heap.getDistanceMax()) || (heap.getCantidad() < totalK)) {
                    //Intersecta asi que busca coincidencia de las keys hacia la derecha
                    if(checkSummaryRight(queryKey, mid, start, end)) {
                        searchBKNN(query, queryKey, mid+1, end, profundidad + 1, heap, totalK);
                    }
                }
            }else {
                //cout << "START: " << start << " == MID: " << mid << " == END: " << end << " == CHECK RIGHT: " << checkSummaryRight(queryKey, mid, start, end) << endl;
                if(checkSummaryRight(queryKey, mid, start, end)) {
                    //Actualiza el Heap
                    searchBKNN(query, queryKey, mid+1, end, profundidad + 1, heap, totalK);
                }
            }
        }
        
        if(query[p] >= nodosKdTree[mid].coordenadas[p]) {
            //Busca por la derecha si existen las palabras claves
            //cout << "START: " << start << " == MID: " << mid << " == END: " << end << " == CHECK RIGHT: " << checkSummaryRight(queryKey, mid, start, end) << endl;
            if(checkSummaryRight(queryKey, mid, start, end)) {
                //Actualiza el Heap
                searchBKNN(query, queryKey, mid+1, end, profundidad + 1, heap, totalK);

                //Busca si intersecta un subespacio mediante la circunferencia
                if ((abs(nodosKdTree[mid].coordenadas[p] - query[p]) < heap.getDistanceMax()) || (heap.getCantidad() < totalK)) {
                    //Intersecta asi que busca coincidencia de las keys hacia la izquierda
                    if(checkSummaryLeft(queryKey, mid, start, end)) {
                        //Actualiza el Heap
                        searchBKNN(query, queryKey, start, mid-1, profundidad + 1, heap, totalK);
                    }
                }
            }else {
                //cout << "START: " << start << " == MID: " << mid << " == END: " << end << " == CHECK LEFT: " << checkSummaryLeft(queryKey, mid, start, end) << endl;
                if(checkSummaryLeft(queryKey, mid, start, end)) {
                    //Actualiza el Heap
                    searchBKNN(query, queryKey, start, mid-1, profundidad + 1, heap, totalK);
                }
            }
        }
    }
    
    //Revisa el punto
    if ((searchEuclideanDistance(query, nodosKdTree[mid].coordenadas) < heap.getDistanceMax()) || (heap.getCantidad() < totalK)) {
        if(checkKeywords(queryKey, mid)) {
            heap.insert(mid, searchEuclideanDistance(query, nodosKdTree[mid].coordenadas));
        }
    }

    return heap;
}

double Build_cBiK::findRKNN(const double* query, const vector<std::string> &queryKey, int k, double alpha) {
    //cout << "**************************** KNN KDTREE ****************************" << endl;
    //struct timeval iniTime2, endTime2;
    double secs;
    
    //Atributos basicos para la busqueda
    int start = 0;
    int end = nodosKdTree.size() - 1;
    string str;
    
    int totalPuntos = nodosKdTree.size();

    //Indice de los keys
    vector<int> indexQueryKey;

    //Se valida k
    if (k < 1) {
        cout << "k debe ser mayor o igual a 1" << endl;
    }else if (k > totalPuntos) {
        cout << "k debe ser menor o igual al total de puntos (" << totalPuntos << ")" << endl;
    }else{
        //Se valida alfa
        if(alpha < 0) {
            cout << "El valor de alfa no puede ser negativo" << endl;
        }else if(alpha > 1) {
            cout << "El valor de alfa no puede ser mayor a 1" << endl;
        }else {
            MinHeap heap = MinHeap(k);
            //gettimeofday(&iniTime2, NULL);
            for(int i=0; i<queryKey.size(); i++) {
                //printTuple(query);
                str = queryKey[i];
                str = toLower(str);
                indexQueryKey.push_back(hashMapKeys[str]);
                //cout << str << " >> " << hashMapKeys[str] << endl;
            }

            heap = searchRKNN(query, indexQueryKey, start, end, heap, k, alpha);
            //gettimeofday(&endTime2, NULL);
            
            //secs = timeval_diff(&endTime2, &iniTime2);
            //Muestra el tiempo por cada consulta
            /*cout << secs * 1000.0 << endl;

            //Se pide el arreglo (HEAP) para mostrar los puntos
            double** elArreglo;
            elArreglo = heap.getArray();

            for (int i = 0; i < k; i++) {
                if(elArreglo[i][1] == numeric_limits<double>::min()) {
                    cout << "\tK=" << (i+1) << " => No existe." << endl;
                    //contFailR++;
                }else {
                    cout << "\tK=" << (i+1);
                    cout << " => Punto: ";
                    printTuple(nodosKdTree[(int)elArreglo[i][0]].coordenadas);
                    cout << " -- ";
                    //printKey(*nodosKdTree[(int)elArreglo[i][0]].palabras);
                    cout << " Score: " << elArreglo[i][1];

                    cout << endl;
                }
            }*/
        }
    }
    
    //cout << endl;
    
    return secs * 1000.0;//retorna las milesimas
}

MinHeap Build_cBiK::searchRKNN(const double* query, vector<int> queryKey, int start, int end, MinHeap &heap, int totalK, double alpha) {
    int mid = (start + end) * 0.5;
    
    //Revisa el puntaje de clasificacion del punto
    double spatialScore = getSpatialScore(query, mid);
    double textualScore = getKeywordsTextualScore(queryKey, mid);
    
    double totalScore = getTotalScore(spatialScore, textualScore, alpha);
    
    //cout << "MID: " << mid << " ++ Punto visitado: ";
    //printTuple(nodosKdTree[mid].coordenadas);
    //cout << " >> SPATIALSCORE: " << spatialScore << " == TEXTUAL SCORE: " << textualScore << " == TOTAL SCORE: " << totalScore;
    //cout << endl;
    
    //Es un nodo interno
    if(start != end) {
        double *summaryTextualScores = getSummariesTextualScores(queryKey, mid, start, end);
        double scoreLeft = getTotalScore(spatialScore, summaryTextualScores[0], alpha);
        double scoreRight = getTotalScore(spatialScore, summaryTextualScores[1], alpha);
        
        //cout << "SCORE LEFT: " << scoreLeft << " == SCORE RIGHT: " << scoreRight << endl;
        
        if(scoreLeft > scoreRight) {//Cuando el score es mayor por la izquierda
            //Ya que es mayor, me aseguro de que existen palabras claves y no es necesario hacer scoreLeft > 0
            //Busco por la izquierda
            searchRKNN(query, queryKey, start, mid-1, heap, totalK, alpha);
            
            //Si no se tienen los k elementos, se busca por la derecha
            //cout << "Cantidad >> " << heap.getCantidad() << endl;
            if((scoreRight > heap.getMinScore()) || (scoreRight > 0 && (heap.getCantidad() < totalK))) {
                searchRKNN(query, queryKey, mid+1, end, heap, totalK, alpha);
            }
        }else if(scoreRight > scoreLeft){//Cuando el score es mayor por la derecha o igual
            //Busco por la derecha
            searchRKNN(query, queryKey, mid+1, end, heap, totalK, alpha);

            //Si no se tienen los k elementos, se busca por la izquierda
            //cout << "Cantidad >> " << heap.getCantidad() << endl;
            if((scoreLeft > heap.getMinScore()) || (scoreLeft > 0 && (heap.getCantidad() < totalK))) {
                searchRKNN(query, queryKey, start, mid-1, heap, totalK, alpha);
            }
        }else {
            //Son iguales, da lo mismo que lado elegir, por defecto derecha
            //Ya que izquierda y derecha son iguales solo me aseguro de que sean distinto a cero
            //Si son cero no es necesario buscar por ningun lado
            if(scoreRight != 0) {
                //Busco por la derecha
                searchRKNN(query, queryKey, mid+1, end, heap, totalK, alpha);
                
                //Como son iguales, necesito buscar por el otro lado ya que existe el mismo score
                searchRKNN(query, queryKey, start, mid-1, heap, totalK, alpha);
            }
        }
    }

    if ((totalScore > heap.getMinScore()) || (heap.getCantidad() < totalK)) {
        //Si es distinto a cero significa q existe al menos un termino buscado
        if(totalScore != 0) {
            heap.insert(mid, totalScore);
        }
    }

    return heap;
}

double Build_cBiK::getTotalScore(double &spatialScore, double &textualScore, double &alpha) {
    //Si no existen ninguna palabra clave se da puntaje cero, se desprecia el puntaje espacial
    return textualScore == 0 ? 0 : ((alpha*spatialScore)+(1-alpha)*(textualScore));
}

double Build_cBiK::getSpatialScore(const double* query, int pos) {
    //cout << setprecision(16) << "SPATIAL SCOOOOOOOOOOOOOOOOOOORE: " << (1 - (searchEuclideanDistance(query, nodosKdTree[pos].coordenadas)/dmax)) << endl;
    return (1 - (searchEuclideanDistance(query, nodosKdTree[pos].coordenadas)/dmax));
}

double Build_cBiK::getKeywordsTextualScore(vector<int> queryKey, int pos) {
    long long int posInicial = pos*numKeys;
    //Relevancia textual
    double scoreText = 0;
    int i=0;
    
    while(i<queryKey.size()){
        if(sdVectorKeywords[(posInicial+queryKey[i])] == 1) {
            scoreText += (((double) 1)/((double) queryKey.size()));
        }
        i++;
    }
    
    return scoreText;
}

double* Build_cBiK::getSummariesTextualScores(vector<int> queryKey, int pos, int start, int end) {
    //orden del arreglo en puntaje izquierdo y derecho respectivamente
    double* textualScores = new double[2];
    textualScores[0] = 0;
    textualScores[1] = 0;
    
    //Es un nodo interno
    if(mapa[pos] == 1) {
        if (end == start + 1) {//Es un nodo interno implicito, se revisa en Keywords
            //NI start
            //HD end
            //El resumen derecho corresponde a la keyword del HD, no tiene resumen izquierdo
            long long int posInicial = end*numKeys;

            int i=0;
            while(i<queryKey.size()){
                //Resumen derecho
                if(sdVectorKeywords[(posInicial+queryKey[i])] == 1) {
                    textualScores[1] += (((double) 1)/((double) queryKey.size()));
                }
                i++;
            }
        } else if (end == start + 2) {//Es un nodo interno implicito, se revisa en Keywords
            //NI start + 1
            //HI start
            //HD end
            //El resumen izquierdo corresponde a la keyword del HI
            long long int posInicialHI = start*numKeys;
            //El resumen derecho corresponde a la keyword del HD
            long long int posInicialHD = end*numKeys;

            int i=0;
            while(i<queryKey.size()){
                //Resumen izquierdo
                if(sdVectorKeywords[(posInicialHI+queryKey[i])] == 1) {
                    textualScores[0] += (((double) 1)/((double) queryKey.size()));
                }
                //Resumen derecho
                if(sdVectorKeywords[(posInicialHD+queryKey[i])] == 1) {
                    textualScores[1] += (((double) 1)/((double) queryKey.size()));
                }
                i++;
            }
        }else if(end > start + 2) {//Es un nodo interno explicito
            long long int posInicial = getNewSummaryIndex(pos);
            
            int i=0;
            while(i<queryKey.size()){
                //Resumen izquierdo
                if(sdVectorResume[(posInicial+queryKey[i])] == 1) {
                    textualScores[0] += (((double) 1)/((double) queryKey.size()));
                }
                //Resumen derecho
                if(sdVectorResume[(posInicial+numKeys+queryKey[i])] == 1) {
                    textualScores[1] += (((double) 1)/((double) queryKey.size()));
                }
                i++;
            }
        }
    }
        
    return textualScores;
}