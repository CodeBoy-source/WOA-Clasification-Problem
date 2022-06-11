/**
 * @file Search.h
 * @version 1.0
 * @date 11/06/2022
 * @author Brian Sena Simons 3ºA-A2
 **/
#ifndef SEARCH_H
#define SEARCH_H

#include "../inc/eigen-3.4.0/Eigen/Dense"
#include "../inc/random.hpp"
#include "Genetics.h"
#include <string.h>
#include <vector>

using namespace std;
using namespace Eigen;
using namespace std::chrono;
using Random = effolkronium::random_static;


void SimulatedAnnealing(MatrixXd data, vector<char>Tlabel, MatrixXd::RowXpr NewSolution,
RowVectorXd& BestSolution, RowVectorXd& best_score, float max_evaluations,
float max_vecinos,float max_exitos,float mu,float phi);

void SimulatedAnnealing(MatrixXd data, vector<char>Tlabel, RowVectorXd& NewSolution,
RowVectorXd& BestSolution, RowVectorXd& best_score, float max_evaluations,
float max_vecinos,float max_exitos,float mu,float phi);

/*
 *@brief Aplicamos búqueda local desde 0 hasta max_eval, con un máximo de vecinos
 visitados igual a maxTilBetter;
    Los valores que devuelven son: El nuevo peso por return, el valor de fitness
    de ese peso y el número de evaluaciones utilizados;
 *@param allData Matriz con los datos.
 *@param label Etiquetas del vector.
 *@param Weight Peso a mejorar.
 *@param eval_num parte de 0 y devuelve el número de evaluaciones obtenidos;
 *@param max_eval número máximo de evaluaciones
 *@param fitness puntuación obtenida por los pesos.
 *@param alpha ponderación de la función.
 */
RowVectorXd LocalSearch(MatrixXd allData,vector<char> label, RowVectorXd Weights,
unsigned int& eval_num, unsigned int max_eval, unsigned int maxTilBetter, vector<float>& fitness, float alpha=0.5, long int seed=1);

RowVectorXd LocalSearchStrong(MatrixXd allData,vector<char> label, RowVectorXd Weights,
unsigned int& eval_num, unsigned int max_eval, unsigned int maxTilBetter, vector<float>& fitness, float alpha=0.5, long int seed=1);

/*
 * @brief Data una matriz de datos con sus etiquetas y una matriz de pesos,
 * para cada fila de la matriz de pesos computamos el valor resultante del 1NN
 * ponderado con los parámetros de reducción y tasa de acierto multiplicados por
 * el valor alpha.
 * @param data matriz con los datos,
 * @param Tlabel vector con las etiquetas,
 * @param Solutions matriz con los pesos
 * @param alpha Ponderación entre reducción y Tasa de Aciertos.
 * @return Devolvemos un vector con la puntuación de cada fila en sus columnas.
 */
RowVectorXd getOnlyFit(MatrixXd data, vector<char> Tlabel, MatrixXd& Solutions,float alpha=0.5);

RowVectorXd getFit(MatrixXd data, vector<char> Tlabel, MatrixXd& Solutions, MatrixXd& GenData, float alpha=0.5);

RowVectorXd get1Fit(MatrixXd data, vector<char> Tlabel, RowVectorXd& Weights, float alpha=0.5 );

#endif
