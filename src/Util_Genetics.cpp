/**
 * @file Util_Genetics.cpp
 * @version 2.3
 * @date 09/05/2022
 * @author Brian Sena Simons 3ºA-A2
 * @brief Funciones diferentes que pueda necesitar para el desarrollo de los
 * algoritmos geneticos
 * @code
 * [...]
    for(i=0,size=Cruzes;i<size;++i){
        BLXCross(Solutions.row(pair),Solutions.row(pair + 1),Cross1,Cross2);
        NewPoblation.row(nuevos) = Cross1;
        NewPoblation.row(nuevos+1) = Cross2;
        NewPoblation.row(nuevos+2) = (Fitness(pair)>Fitness(pair+1))? Solutions.row(pair) : Solutions.row(pair+1);
        pair+=2;
        nuevos+=3;
    }
 *[...]
 * @endcode
 **/
#include "../tools/Genetics.h"
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
 * @brief Gracias a la librería Eigen, podemos simplesmente sumar los padres y
 * multiplicar por un scalar generado aleatoriamente entre 0 y 1; Y es como
 * si estuvieramos haciendo un bucle for y sumando componente a componente, multiplicando
 * por ese alpha fijo.
 */
void ArithmeticCross(RowVectorXd parent1, RowVectorXd parent2, RowVectorXd& res1, RowVectorXd& res2, long int seed){
    if(seed!=-1)
        Random::seed(seed);

    if(res1.cols() != parent1.cols())
        res1.resize(parent1.cols());
    if(res2.cols() != parent1.cols())
        res2.resize(parent2.cols());

    res1 = (parent1+parent2)*Random::get(0,1);
    res2 = (parent1+parent2)*Random::get(0,1);
}

/*
 * @brief Aquí no hay escapatoria, tenemos que hacer un bucle for para ir
 * mirando el máximo y el mínimo componente a componente para poder calcular el
 * intervalo a partir de la distancia entre esos y el valor alpha. Una vez tengamos
 * el valor mínimo LOW y el máximo HIGH entonces podemos llamar al generador
 * aleatorio para que nos genere un valor dentro del intervalo.
 */
void BLXCross(RowVectorXd parent1, RowVectorXd parent2,RowVectorXd& res1, RowVectorXd& res2, float alpha, long int seed){
    if(seed!=-1)
        Random::seed(seed);

    float Cmin, Cmax, Distance, LO, HI;

    if(res1.cols() != parent1.cols())
        res1.resize(parent1.cols());
    if(res2.cols() != parent1.cols())
        res2.resize(parent2.cols());

    for(unsigned int i=0;i < parent1.cols(); i++){
        Cmin = min(parent1(i),parent2(i));
        Cmax = max(parent1(i),parent2(i));
        Distance = Cmax - Cmin;
        LO = Cmin - Distance * alpha;
        HI = Cmax + Distance * alpha;

        res1(i) = Random::get(LO,HI);
        res2(i) = Random::get(LO,HI);
    }
}
