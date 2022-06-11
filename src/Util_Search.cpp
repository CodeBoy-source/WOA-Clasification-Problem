/**
 * @file Util_Search.cpp
 * @version 1.0
 * @date 11/06/2022
 * @author Brian Sena Simons 3ºA-A2
 **/
#include "../tools/Search.h"
#include "../inc/eigen-3.4.0/Eigen/Dense"
#include "../inc/random.hpp"
#include <fstream>
#include <math.h>
#include <iostream>

using namespace std;
using namespace Eigen;
using namespace std::chrono;
using Random = effolkronium::random_static;

/*
 * @brief Búsqueda local utilizada en la práctica anterior. La idea es ir
 * modificando el vector de soluciones Weights con una distribución normal
 * de media 0 y varianza 0.3 aleatoriamente y aceptando los cambios cuando
 * se produce una mejora. El máximo a evaluar es pasado por parámetro y el
 * máximo de vecinos a explorar también (20*n_columnas).
 */
RowVectorXd LocalSearch(MatrixXd allData,vector<char> label, RowVectorXd Weights,
unsigned int& eval_num, unsigned int max_eval, unsigned int maxTilBetter, vector<float>& fitness, float alpha,long int seed){
    Random::seed(seed);
    eval_num = 0;
    vector<int> indexGrid, Evaluations;
    fillRange(indexGrid,allData.cols());
    std::normal_distribution<double> distribution(0.0, sqrt(0.3));
    unsigned int tilBetter = 0, ran_num = 0, reduct, right;
    RowVectorXd Weights_before(allData.cols());
    float function_after=0, function_before=0, reduction_rate, right_rate;
    if(fitness.size()<=1){
        fitness.clear();
        fitness.push_back(0);
        fitness.push_back(0);
    }else{
        function_before = fitness[0] + fitness[1];
    }
    while(eval_num < max_eval and tilBetter<maxTilBetter ){
        /// Modificamos el indice barajado correspondiente para modificar los pesos.
        Weights[indexGrid[ran_num++]] += Random::get(distribution);
        /// Si hemos modificados volvemos a barajar.
        if(ran_num>=indexGrid.size()){
            Random::shuffle(indexGrid);
            ran_num = 0;
        }

        getReductRight(allData,label,Weights,right,reduct);

        /// Evaluamos la función.
        reduction_rate = (1-alpha)*(float(reduct)/float(Weights.cols()));
        right_rate = alpha*(float(right)/float(allData.rows()));
        function_after = right_rate + reduction_rate;
        eval_num++;

        /// Verificamos si hemos mejorado, actuamos acorde.
        if(function_after <= function_before){
            tilBetter++;
            Weights = Weights_before;
        }else{
            fitness.at(0) = right_rate;
            fitness.at(1) = reduction_rate;
            tilBetter = 0;
            Weights_before = Weights;
            function_before = function_after;
            Random::shuffle(indexGrid);
            ran_num = 0;
        }

    } // END WHILE
    return Weights;

}

RowVectorXd LocalSearchStrong(MatrixXd allData,vector<char> label, RowVectorXd Weights,
unsigned int& eval_num, unsigned int max_eval, unsigned int maxTilBetter, vector<float>& fitness, float alpha,long int seed){
    Random::seed(seed);
    eval_num = 0;
    vector<int> Evaluations;
    std::normal_distribution<double> distribution(0.0, sqrt(0.4));
    unsigned int tilBetter = 0, reduct, right;
    RowVectorXd Weights_before(allData.cols());
    float function_after=0, function_before=0, reduction_rate, right_rate;
    if(fitness.size()<=1){
        fitness.clear();
        fitness.push_back(0);
        fitness.push_back(0);
    }else{
        function_before = fitness[0] + fitness[1];
    }
    unsigned int amount = 0;
    while(eval_num < max_eval and tilBetter<maxTilBetter ){
        amount = Random::get<unsigned>(1,Weights.cols()-1);
        for(unsigned int i=0;i<amount;i++){
            Weights[Random::get<unsigned>(0,Weights.cols()-1)] += Random::get(distribution);
        }
        getReductRight(allData,label,Weights,right,reduct);

        /// Evaluamos la función.
        reduction_rate = (1-alpha)*(float(reduct)/float(Weights.cols()));
        right_rate = alpha*(float(right)/float(allData.rows()));
        function_after = right_rate + reduction_rate;
        eval_num++;

        /// Verificamos si hemos mejorado, actuamos acorde.
        if(function_after <= function_before){
            tilBetter++;
            Weights = Weights_before;
        }else{
            fitness.at(0) = right_rate;
            fitness.at(1) = reduction_rate;
            tilBetter = 0;
            Weights_before = Weights;
            function_before = function_after;
        }

    } // END WHILE
    return Weights;

}

/*
 * @brief La idea es mutar un número de veces determinado por la esperanza
 * matemática de la mutación. Y para aumentar la diversidad permito que se
 * mute más de una componente de una fila dada. Además guardo las posiciones que
 * se mutaron por si tengo que volver a calcular su fitness.
 */
void Mutate(MatrixXd* NP2, vector<int>&indexGrid,unsigned int Mutacion){
    std::normal_distribution<double> distribution(0.0, sqrt(0.3));
    int Rest = Mutacion; unsigned int Diversidad=0, i=0;
    indexGrid.clear();
    while(Rest > 0){
        i = Random::get<unsigned>(0,NP2->rows()-1);
        Diversidad = Random::get<unsigned>(1,NP2->cols()-1);
        for(unsigned int j=0;j<Diversidad;++j){
            (*NP2)(i,Random::get<unsigned>(0,NP2->cols()-1)) += Random::get(distribution);
        }

        indexGrid.push_back(i);
        Rest = Rest - Diversidad;
    }
}

/*
 * @brief Función que cruzan la matriz de forma aleatoria, cruzando los N primeros
 * después de barajar.
 */
int randomOnly(MatrixXd data, vector<char> Tlabel, MatrixXd* P1,MatrixXd* NP2,
        MatrixXd& GenData,int CrossType, unsigned int Cruzes,unsigned int Mutacion){

    MatrixXd NewGenData(GenData.rows(),GenData.cols());
    unsigned int pair = 0,nuevos=0,i=0,size=0,cont=0,parent=0;
    vector<int> indexGrid;
    shuffleFit(*P1,GenData,-1);
    RowVectorXd Cross1, Cross2;
    for(cont=2*Cruzes;cont<P1->rows();cont++){
        NP2->row(cont) = P1->row(cont);
        NewGenData.row(cont) = GenData.row(cont);
    }
    for(i=0,size=Cruzes;i<size;++i){
        if(CrossType == 0)
            BLXCross(P1->row(pair),P1->row(pair+1),Cross1,Cross2);
        else
            ArithmeticCross(P1->row(pair),P1->row(pair+1),Cross1,Cross2);

        NP2->row(nuevos) = Cross1;
        NP2->row(nuevos+1) = Cross2;

        pair+=2;
        nuevos+=2;
    }

    for(parent=2*Cruzes;parent<P1->rows();parent++){
        GenData.row(parent) = NewGenData.row(parent);
    }

    Mutate(NP2, indexGrid,Mutacion);
    unsigned int eval_cont = 0;
    nuevos = 0;
    for(i=0,size=Cruzes;i<size;++i){
        Cross1 = NP2->row(nuevos);
        GenData.row(nuevos) = get1Fit(data,Tlabel,Cross1,0.5);
        NP2->row(nuevos) = Cross1;
        Cross2 = NP2->row(nuevos+1);
        GenData.row(nuevos+1) = get1Fit(data,Tlabel,Cross2,0.5);
        NP2->row(nuevos+1) = Cross2;
        eval_cont += 2;
        nuevos+=2;
    }
    for(i=0,size=indexGrid.size();i<size;++i){
        if(indexGrid[i]>=int(2*Cruzes)){
            Cross1 = NP2->row(indexGrid[i]);
            GenData.row(indexGrid[i]) = get1Fit(data,Tlabel,Cross1,0.5);
            NP2->row(indexGrid[i]) = Cross1;
            eval_cont++;
        }
    }
    return eval_cont;
}

/*
 * @brief función que cruza los N primeros mejores padres de la matriz de soluciones
 */
int onlyBestCrossing(MatrixXd data, vector<char> Tlabel, MatrixXd* P1,MatrixXd* NP2,
        MatrixXd& GenData,int CrossType, unsigned int Cruzes,unsigned int Mutacion){

    MatrixXd NewGenData(GenData.rows(),GenData.cols());
    unsigned int pair = 0,nuevos=0,i=0,size=0,cont=0,parent=0;
    vector<int> indexGrid;
    getBest(GenData,indexGrid,P1->rows());
    RowVectorXd Cross1, Cross2;
    for(cont=2*Cruzes;cont<P1->rows();cont++){
        NP2->row(cont) = P1->row(indexGrid[cont]);
        NewGenData.row(cont) = GenData.row(indexGrid[cont]);
    }
    for(i=0,size=Cruzes;i<size;++i){
        if(CrossType == 0)
            BLXCross(P1->row(indexGrid[pair]),P1->row(indexGrid[pair+1]),Cross1,Cross2);
        else
            ArithmeticCross(P1->row(indexGrid[pair]),P1->row(indexGrid[pair+1]),Cross1,Cross2);

        NP2->row(nuevos) = Cross1;
        NP2->row(nuevos+1) = Cross2;

        pair+=2;
        nuevos+=2;
    }

    for(parent=2*Cruzes;parent<P1->rows();parent++){
        GenData.row(parent) = NewGenData.row(parent);
    }

    Mutate(NP2, indexGrid,Mutacion);
    unsigned int eval_cont = 0;
    nuevos = 0;
    for(i=0,size=Cruzes;i<size;++i){
        Cross1 = NP2->row(nuevos);
        GenData.row(nuevos) = get1Fit(data,Tlabel,Cross1,0.5);
        NP2->row(nuevos) = Cross1;
        Cross2 = NP2->row(nuevos+1);
        GenData.row(nuevos+1) = get1Fit(data,Tlabel,Cross2,0.5);
        NP2->row(nuevos+1) = Cross2;
        eval_cont += 2;
        nuevos+=2;
    }
    for(i=0,size=indexGrid.size();i<size;++i){
        if(indexGrid[i]>=int(2*Cruzes)){
            Cross1 = NP2->row(indexGrid[i]);
            GenData.row(indexGrid[i]) = get1Fit(data,Tlabel,Cross1,0.5);
            NP2->row(indexGrid[i]) = Cross1;
            eval_cont++;
        }
    }
    return eval_cont;
}

/*
 * @brief función que cruza de forma aleatorioa los N primeros padres y  además
 * guarda el padre de la pareja con mejor puntuación aumentando la presión social.
 */
int randomCrossKeepBest(MatrixXd data, vector<char> Tlabel, MatrixXd* P1,MatrixXd* NP2,
        MatrixXd& GenData,int CrossType, unsigned int Cruzes,unsigned int Mutacion){

    unsigned int pair = 0,nuevos=0,i=0,size=0,parent=0;
    vector<int> indexGrid, randomIndex;
    shuffleFit(*P1,GenData,-1);
    MatrixXd NewGenData = GenData;
    RowVectorXd Cross1, Cross2;
    bool full = false;
    for(i=0,size=Cruzes;i<size&&!full;++i){
        if(CrossType == 0)
            BLXCross(P1->row(pair),P1->row(pair+1),Cross1,Cross2);
        else
            ArithmeticCross(P1->row(pair),P1->row(pair+1),Cross1,Cross2);

        NP2->row(nuevos) = Cross1;
        NP2->row(nuevos+1) = Cross2;

        NP2->row(nuevos+2) = (GenData.row(pair).sum()>GenData.row(pair+1).sum())? P1->row(pair): P1->row(pair+1);
        NewGenData.row(nuevos+2) = GenData.row(parent);
        pair+=2;
        if(nuevos+5<NP2->rows()) nuevos+=3;
        else full=true;
    }
    if(full){
        for(;nuevos<NP2->rows();nuevos++){
            NP2->row(nuevos) = P1->row(nuevos);
            NewGenData.row(nuevos) = GenData.row(nuevos);
        }
    }

    Mutate(NP2,indexGrid,Mutacion);

    full = false;
    unsigned int eval_cont = 0;
    nuevos = 0;
    for(i=0,size=Cruzes;i<size&&!full;++i){
        Cross1 = NP2->row(nuevos);
        NewGenData.row(nuevos) = get1Fit(data,Tlabel,Cross1,0.5);
        NP2->row(nuevos) = Cross1;
        Cross2 = NP2->row(nuevos+1);
        NewGenData.row(nuevos+1) = get1Fit(data,Tlabel,Cross2,0.5);
        NP2->row(nuevos+1) = Cross2;

        eval_cont += 2;
        if(nuevos+5>NP2->rows()) full=true;
        else nuevos+=3;
    }
    for(i=0,size=indexGrid.size();i<size;++i){
        if(indexGrid[i]%3==2){
            Cross1 = NP2->row(indexGrid[i]);
            NewGenData.row(indexGrid[i]) = get1Fit(data,Tlabel,Cross1,0.5);
            NP2->row(indexGrid[i]) = Cross1;
            ++eval_cont;
        }
    }
    GenData = NewGenData;
    return eval_cont;
}


void SimulatedAnnealing(MatrixXd data, vector<char>Tlabel, MatrixXd::RowXpr NewSolution,
RowVectorXd& BestSolution, RowVectorXd& best_score, float max_evaluations,
float max_vecinos,float max_exitos,float mu,float phi){
    SimulatedAnnealing(data,Tlabel,NewSolution,BestSolution,best_score,max_evaluations,max_vecinos,max_exitos,mu,phi);
}

void SimulatedAnnealing(MatrixXd data, vector<char>Tlabel, RowVectorXd& NewSolution,
RowVectorXd& BestSolution, RowVectorXd& best_score, float max_evaluations,
float max_vecinos,float max_exitos,float mu,float phi){
    unsigned int cols = NewSolution.cols();
    if(max_vecinos==0)
        max_vecinos = 10*cols;
    if(max_exitos==0)
        max_exitos = 0.1*max_vecinos;


    std::normal_distribution<double> distribution(0.0, sqrt(0.4));

    float T_0=0.1, T_f=0.001;
    RowVectorXd score(2), old_score(2), sabest(2),
                Solution(cols),SABest(cols);

    old_score = get1Fit(data,Tlabel,NewSolution);
    T_0 = ( mu * old_score.sum()) / (-1*log(phi) );
    if(T_0 == 0)
        T_0 = 0.1;
    if(T_0 <= T_f)
        T_f = T_0 * pow(10,-3);

    float M = max_evaluations/max_vecinos;
    float beta = (T_0 - T_f) / (M*T_0*T_f),diff;
    unsigned int vecinos = 0, exitos = 0, evaluations=0, index;
    SABest = Solution = NewSolution;
    sabest = old_score;
    evaluations = 1;
    while(evaluations<=max_evaluations){
        vecinos = exitos = 0;
        while(vecinos<=max_vecinos && exitos<=max_exitos && evaluations<=max_evaluations){
            index = Random::get<unsigned>(0,cols-1);
            // Apply move
            NewSolution[index] += Random::get(distribution);
            score = get1Fit(data,Tlabel,NewSolution,0.5);
            vecinos++;evaluations++;
            diff = old_score.sum()-score.sum();
            if((diff<0.0) || (Random::get(0.0,1.0) <= 1.0/exp(diff/T_0))){
                Solution[index] = NewSolution[index];
                old_score = score;
                vecinos = 0;
                exitos++;
                if(sabest.sum() < score.sum()){
                    SABest = Solution;
                    sabest = score;
                }
            }
            NewSolution[index] = Solution[index];
        }// END WHILE2
        T_0 = (T_0) / (1+beta*T_0);
        if(exitos==0) break;
    }
    if(best_score.sum()<sabest.sum()){
        BestSolution = SABest;
        best_score = sabest;
    }
}

/*
 * @brief Esta función nos calcula la reducción modificando el vector de pesos y además
 * calcula la tasa de clasificación y nos devuelve el valor por los parámetros por
 * pasados por referencia reduct y right.
 */
void getReductRight(MatrixXd data, vector<char> Tlabel, RowVectorXd& Weights, unsigned int &right, unsigned int &reduct){
        unsigned int i,size, ManualNeighbour;
        reduct = 0;
        for(i=0,size=Weights.cols();i<size;i++){
            if(*(Weights.data() + i) < 0.1){
                *(Weights.data() + i) = 0;
                reduct++;
            }
            else if(*(Weights.data()+i)>1)
                *(Weights.data()+i) = 1;
        }
        /// Verificamos el porcentaje de acierto obtenido con la modificación.
        right = 0;
        for(i=0,size=data.rows();i<size;i++){
            ManualEuclideanDistance(Weights,data.row(i),data,i, ManualNeighbour);
            if(Tlabel[i] == Tlabel[ManualNeighbour])
                right++;
        }
}

/*
 * @brief Esta función nos devuelve apenas un vector con las funciones agregadas
 * pero como al final quería generar las tablas y estudiar más a hondo lso resultados
 * he cambiado por una matriz con la tasa de clasificación y la reducción.
 */
RowVectorXd getOnlyFit(MatrixXd data, vector<char> Tlabel, MatrixXd& Solutions,float alpha ){
    unsigned int right=0,  reduct = 0;
    unsigned int j;
    RowVectorXd Weights;
    RowVectorXd results(Solutions.rows());

    for(j=0;j<Solutions.rows();j++){
        /// Truncamos los pesos a 0 y 1; Contamos las reducciones a 0.
        Weights = Solutions.row(j);
        getReductRight(data,Tlabel,Weights,right,reduct);
        Solutions.row(j) = Weights;
        results(j) = alpha*(float(right)/float(data.rows())) + (1-alpha)*(float(reduct)/float(Weights.cols()));
    }
    return results;
}

/*
 * @brief esta es la actualización de getOnlyFit() explicado anteriormente. La
 * idea es la misma, obtener el fitness. Sin embargo, aquí lo almacenamos en una
 * matriz GenData de 2 columnas, donde la primera es la tasa de clasificación
 * multiplicada por alpha y la segunda es la de reducción multiplicada por 1-alpha.
 */
RowVectorXd getFit(MatrixXd data, vector<char> Tlabel, MatrixXd& Solutions, MatrixXd& GenData, float alpha ){
    unsigned int right=0,  reduct = 0;
    unsigned int j;
    RowVectorXd Weights;
    RowVectorXd results(Solutions.rows());
    if(GenData.rows()!=Solutions.rows() || GenData.cols()!=2)
        GenData.resize(Solutions.rows(),2);

    for(j=0;j<Solutions.rows();j++){
        /// Truncamos los pesos a 0 y 1; Contamos las reducciones a 0.
        Weights = Solutions.row(j);
        getReductRight(data,Tlabel,Weights,right,reduct);
        Solutions.row(j) = Weights;
        GenData(j,0) = alpha*(float(right)/float(data.rows()));
        GenData(j,1) = (1-alpha)*(float(reduct)/float(Weights.cols()));
        results(j) = GenData(j,0) + GenData(j,1);
    }
    return results;
}

/*
 * @brief No siempre quiero calcular el fitness de toda la matriz, por ello
 * he creado esta función que calcula apenas el fitness de 1 solución
 */
RowVectorXd get1Fit(MatrixXd data, vector<char> Tlabel, RowVectorXd& Weights, float alpha ){
    unsigned int right=0,  reduct = 0;
    RowVectorXd results(2);
    getReductRight(data,Tlabel,Weights,right,reduct);
    results(0) = alpha*(float(right)/float(data.rows()));
    results(1) = (1-alpha)*(float(reduct)/float(Weights.cols()));
    return results;
}

