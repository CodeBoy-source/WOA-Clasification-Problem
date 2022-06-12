#include "../inc/eigen-3.4.0/Eigen/Dense"
// From https://github.com/effolkronium/random
#include "../inc/random.hpp"

#include "../tools/mytools.h"
#include "../tools/Euclidean.h"
#include "../tools/ReadData.h"
#include "../tools/Genetics.h"
#include "../tools/Search.h"
#include <vector>
#include <string>
#include <iostream>
#include <chrono>
#include <thread>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <fstream>
#include <unistd.h>


using namespace std;
using namespace Eigen;
using namespace std::chrono;
using Random = effolkronium::random_static;


int main(int argc, char** argv){
    if(argc<=7) {
        cerr << "[ERROR]: Couldn't resolve file name;" << endl;
        cerr << "[EXECUTION]: ./main (filename) (label1) (label2) (0-print,1-write) (seed)[int] (0-Normal, 1-ShuffleData,2-BalanceData) (PopSize)[+int]" << endl;
        cerr << "[EXAMPLE]: ./main parkinsons.arff 1 2 1 150421 0" << endl;
        exit(-1);
    }
    bool debuggin = false;
    string filename = argv[1];
    char type1 = *argv[2];
    char type2 = *argv[3];
    int streambus = atoi(argv[4]);
    long int seed = atol(argv[5]);
    int shuffle = atoi(argv[6]);
    unsigned int popSize = atoi(argv[7]);
    srand(seed);
    Random::seed(seed);
    bool printing = (streambus>=1)?false:true;
    ofstream plot,myfile;
    string writefile = "", plot_path = "", output;
    string path = get_selfpath();
    path = path.substr(0,path.find_last_of("/\\") + 1);
    bool plotting = false;
    if(streambus>=1) {
        //https://www.codegrepper.com/code-examples/cpp/c%2B%2B+get+filename+from+path
        // get filename
        std::string base_filename = filename.substr(filename.find_last_of("/\\") + 1);
        // remove extension from filename
        std::string::size_type const p(base_filename.find_last_of('.'));
        std::string file_without_extension = base_filename.substr(0, p);

        string datafilename = "WOA_"+to_string(popSize)+file_without_extension+to_string(seed)+"_"+to_string(shuffle);
        writefile = path+"../results/"+datafilename;
        writefile += ".txt";
        myfile.open(writefile,ios::out|ios::trunc);
        if(!myfile.is_open()){
            cerr << "[ERROR]: Couldn't open file, printing enabled" << endl;
            printing = true;
        }
        if(streambus>1){
            plotting = true;
            if(plotting){
                plot_path = path+"../results/plots/"+datafilename+".txt";
                plot.open(plot_path,ios::out|ios::trunc);
                if(!plot.is_open()){
                    cerr << "[ERROR]: Couldn't open plot file, action dismissed\n";
                    plotting = false;
                }
            }
        }
    }
    filename = path+"../datos/"+filename;

    vector<char> label;
    MatrixXd allData = readValues(filename,label);

    std::normal_distribution<double> distribution(0.0, sqrt(0.3));

    high_resolution_clock::time_point momentoInicio, momentoFin;
    milliseconds tiempo;

    MatrixXd data, test, group1, group2;
    vector<char> Tlabel, Ttlabel, label_group1, label_group2;
    /// Inicializamos todas las variables que vamos a necesitar para almacenar información
    if(shuffle==1){
        cout << "[WARNING]: Data has been shuffled; " << endl;
        if(streambus>=1)
            myfile << "[WARNING]: Data has been shuffled; " << endl;
        shuffleData(allData,label,seed);
    }
    if(shuffle==2){
        cout << "[WARNING]: Data has been balanced and shuffled (inevitably); " << endl;
        if(streambus>=1)
            myfile << "[WARNING]: Data has been balanced and shuffled (inevitably); " << endl;
        group1 = getClassLabelled(allData,label, label_group1, type1);
        group2 = getClassLabelled(allData,label, label_group2, type2);
    }

    if(myfile.is_open()){
        myfile << " ### Algoritmo WOA ###\n ";
        myfile << "F\tclasific\treducir \tfitness \ttime\n";
    }

    MatrixXd WhalePop(popSize,allData.cols()), WhaleFit(popSize,2);
    MatrixXd::Index maxIndex,minIndex;
    RowVectorXd bestWhale(allData.cols()),Direction(allData.cols()),bestFit(2);
    int maxcols = allData.cols(), evaluations, maxEvaluations = 15000, iteration, maxIterations;
    unsigned int i;
    float alpha = 0.5, a, a2, A,C,b,l,p;
    vector<double> randomA,randomC,randomL;
    for(int x=0,folds=5;x<folds;x++){
        momentoInicio = high_resolution_clock::now();
        if(shuffle!=2)
            getFold(allData,label,data,Tlabel,test,Ttlabel,x);
        else
            getBalancedFold(group1,label_group1,group2,label_group2,data,Tlabel, test, Ttlabel,x,seed);
        // Generate random population
        WhalePop = (MatrixXd::Random(popSize,maxcols) + MatrixXd::Constant(popSize,maxcols,1))/2.0;
        getFit(data,Tlabel,WhalePop,WhaleFit,alpha);
        // Initialize Best Whale as worst
        WhaleFit.rowwise().sum().maxCoeff(&maxIndex);
        bestWhale = WhalePop.row(maxIndex);
        bestFit = WhaleFit.row(maxIndex);
        iteration = 0;
        evaluations = popSize + 1;
        maxIterations = (maxEvaluations-evaluations)/popSize;
        // ############################################################
        while(evaluations < maxEvaluations){
            a = 2.0 * (1.0-float(iteration)/float(maxIterations));
            a2 = -1.0+float(iteration)*((-1.0)/float(maxIterations));
            randomA = Random::get<std::vector>(0.0,1.0,popSize);
            randomC = Random::get<std::vector>(0.0,1.0,popSize);
            randomL = Random::get<std::vector>(0.0,1.0,popSize);
            for(i=0;i<popSize;i++){
                A = 2 * a * randomA[i] - a;
                C = 2 * randomC[i];
                b = 1;
                l = (a2-1)*randomL[i]+1;
                p = Random::get();
                if(p<0.5){
                    if(abs(A)>1){
                       Direction = (C * bestWhale - WhalePop.row(i)).array().abs();
                       WhalePop.row(i) = bestWhale - A * Direction;
                    }else{
                        do{
                        p = Random::get<unsigned>(0,popSize-1);
                        }while(p==i);

                        Direction = (C * WhalePop.row(p) - WhalePop.row(i)).array().abs();
                        WhalePop.row(i) = WhalePop.row(p) - A * Direction;
                    }
                }else{
                    Direction = (bestWhale - WhalePop.row(i)).array().abs();
                    WhalePop.row(i) = Direction * exp(b*l) * cos(2*M_PI*l) + bestWhale;
                }
            }
            // This takes care of clipping the values from the MatrixSolutions as well
            getFit(data,Tlabel,WhalePop,WhaleFit,alpha);
            evaluations += popSize;
            if(bestFit.sum() < WhaleFit.rowwise().sum().maxCoeff(&maxIndex)){
                bestWhale = WhalePop.row(maxIndex);
                bestFit = WhaleFit.row(maxIndex);
            }

            ++iteration;
            if(debuggin){
                WhaleFit.rowwise().sum().maxCoeff(&maxIndex);
                WhaleFit.rowwise().sum().minCoeff(&minIndex);
                cout << "###################################\n" ;
                cout << "[ITERATION NUMBER]: " << iteration  << " - " << maxIterations << "\n";
                cout << "[BEST FITNESS]: " << WhaleFit(maxIndex,0)/alpha << "\t" << WhaleFit(maxIndex,1)/(1-alpha) << "\t" << WhaleFit.row(maxIndex).sum() << "\n";
                cout << "[WORST FITNESS]: " << WhaleFit(minIndex,0)/alpha << "\t" << WhaleFit(minIndex,1)/(1-alpha) << "\t" << WhaleFit.row(minIndex).sum() << "\n";
                cout << "[BEST SAVED]: " << bestFit(0)/alpha << "\t" << bestFit(1)/(1-alpha) << "\t" << bestFit.sum() << "\n";
                //cout << WhalePop << endl;
            }else{
                progress_bar(float(x*maxEvaluations + evaluations)/float(maxEvaluations*folds));
            }
            if(plotting){
                WhaleFit.rowwise().sum().maxCoeff(&maxIndex);
                WhaleFit.rowwise().sum().minCoeff(&minIndex);
                output = to_string(iteration) + "\t" + to_string(WhaleFit(maxIndex,0)/alpha) +
                    "\t" + to_string(WhaleFit(maxIndex,1)/(1-alpha)) + "\t" + to_string(WhaleFit.row(maxIndex).sum());
                plot << std::setw(31) << output;
                output = "\t" + to_string(WhaleFit(minIndex,0)/alpha) + "\t" + to_string(WhaleFit(minIndex,1)/(1-alpha)) +
                    "\t" + to_string(WhaleFit.row(minIndex).sum()) + "\n";
                plot << std::setw(30) << output;
            }
        }
        // ############################################################
        momentoFin = high_resolution_clock::now();
        tiempo = duration_cast<milliseconds>(momentoFin - momentoInicio);
        bestFit = get1Fit(test,Ttlabel, bestWhale,alpha);
        if(printing){
            cout << "###################################\n" ;
            cout << "[BEST FITNESS]: " << bestFit(0)/alpha << "\t" << bestFit(1)/(1-alpha) << "\t" << bestFit.sum() << "\t";
            cout << tiempo.count() << endl;
        }else{
            output = to_string(x) + "\t" + to_string(bestFit(0)/alpha)
                    + "\t" +to_string(bestFit(1)/(1-alpha)) + "\t" +
                    to_string(bestFit.sum()) + "\t" + to_string(tiempo.count()) + "\n";
            myfile << std::setw(30) << output;
        }
    }
    cout << endl;
    return 0;
}
