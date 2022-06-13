#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "../tools/mytools.h"
#include <filesystem>
#include <bits/stdc++.h>

using namespace std;
namespace fs = filesystem;

struct parameters{
    float average;
    float lambda;
    float delta;
    float popSize;
} ;

void printMax(vector<parameters> average, string filename,unsigned int n_top){
    string output;
    cout << "############# " + filename + " ##############" << endl;
    for(unsigned int i=0;i<n_top;i++){
        auto index = max_element(average.begin(),
                average.end(),
                [](const auto& lhs, const auto& rhs) { return lhs.average < rhs.average; });
        output = "Top-" +to_string(i)+":" + "\tAvg:" +  to_string(index->average) +
            "\tLambda:" + to_string(index->lambda) +"\tDelta:" + to_string(index->delta) +
            "\tpopSize:" + to_string(index->popSize) + "\n";
        cout << std::setw(30) << output;
        average.erase(index);
        if(average.empty())
            break;
    }
}

void getBest(parameters& best, vector<vector<parameters>> total_average){
    best.average = 0;
    float avg_fit = 0;
    for(unsigned int j=0;j<total_average[0].size();j++){
        avg_fit += total_average[0][j].average;
        avg_fit += total_average[1][j].average;
        avg_fit += total_average[2][j].average;
        avg_fit/=3.0;
        if(best.average < avg_fit){
            best.average = avg_fit;
            best.lambda = total_average[0][j].lambda;
            best.delta = total_average[0][j].delta;
            best.popSize = total_average[0][j].popSize;
        }
        avg_fit = 0;
    }
}

int main(int argc, char** argv){
    if(argc<=3){
        cerr << "[ERROR]: Couldn't resolve file name" << endl;
        cerr << "[SOLVE]: You must specify the amount of folds and, optinually, maximum values to find and popSize" << endl;
        cerr << "[EXEC]: ./main 5 3 500" << endl;
    }
    unsigned int KFolds = atoi(argv[1]);
    unsigned int n_top = 2;
    if(argc>2){
        n_top = atoi(argv[2]);
    }
    unsigned int popSize = atoi(argv[3]);
    bool debuggin = false;
    if(argc>4){
        debuggin = (atoi(argv[4])==1)?true:false;
    }
    vector<parameters> average;
    if(argc<=5){
        string filename = "";
        string path = get_selfpath();
        path = path.substr(0,path.find_last_of("/\\") + 1);
        path = path +"../results/HSearch/";
        //https://www.codegrepper.com/code-examples/cpp/c%2B%2B+get+filename+from+path
        // get filename
        std::string base_filename = filename.substr(filename.find_last_of("/\\") + 1);
        filename = path + base_filename;
        string read = "", output;
        ifstream myfile;
        vector<vector<parameters>> total_average;
        for( const auto& entry: fs::directory_iterator(path)){
            read = entry.path();
            filename = read.substr(read.find_last_of("/\\") + 1 );
            ofstream testfile;
            if(filename.find("WOAH_"+to_string(popSize))!=string::npos and filename.find("HS")!=string::npos){
                myfile.open(path+filename);
                filename = filename.substr(0,filename.find("."));
                if(!myfile.is_open()){
                    cerr << "[ERROR]: Unable to open file at " + filename << endl;
                }
                float fold, fitness, lambda,delta, avg_fit,popSize;
                while(!myfile.eof()){
                    avg_fit = 0;
                    for(unsigned int i=0;i<KFolds;i++){
                        myfile >> fold >> fitness >> lambda >> delta >> popSize;
                        if(fold!=i)
                            break;
                        else
                            avg_fit += fitness;
                    }
                    avg_fit /= float(KFolds);
                    if(avg_fit!=0)
                        average.push_back({avg_fit,lambda,delta, popSize});
                }

                if(debuggin)
                    printMax(average,filename,n_top);
                total_average.push_back(average);
                average.clear();
                myfile.close();
            }
        }
        parameters best;
        if(debuggin)
            cout << "Total Lines computed: " << total_average[0].size()*5 << endl;
        getBest(best,total_average);
        cout << "[BEST]: " << best.average << " - " << best.lambda << " - " << best.delta << " - " << best.popSize << endl;
    }else{
        string filename = argv[5];
        string output;
        string path = get_selfpath();
        path = path.substr(0,path.find_last_of("/\\") + 1);
        path = path +"../results/";
        //https://www.codegrepper.com/code-examples/cpp/c%2B%2B+get+filename+from+path
        std::string base_filename = filename.substr(filename.find_last_of("/\\") + 1);
        filename = path + base_filename;
        ifstream myfile(filename);
        if(!myfile.is_open()){
            cerr << "[ERROR]: Unable to open file at " + filename << endl;
        }
        if(filename.find("HS")==string::npos){
            // READ HEADERS
            getline(myfile,output);
            getline(myfile,output);

            float fold, fitness, useless, avg_fit;
            while(!myfile.eof()){
                avg_fit = 0;
                for(unsigned int i=0;i<KFolds;i++){
                    myfile >> fold >> useless >> useless >> fitness >> useless;
                    if(fold!=i)
                        break;
                    else
                        avg_fit += fitness;
                }
                avg_fit /= float(KFolds);
                if(avg_fit !=0){
                    average.push_back({avg_fit,0,0,0});
                }
            }
            printMax(average,filename.substr(filename.find_last_of("/")+1),n_top);
        }else{
            cout << "Is it using the old_version?[0:NO,1:YES] ";
            string old;
            cin >> old;
            if(old=="0"){
                float fold, fitness, lambda,delta, avg_fit,popSize;
                while(!myfile.eof()){
                    avg_fit = 0;
                    for(unsigned int i=0;i<KFolds;i++){
                        myfile >> fold >> fitness >> lambda >> delta >> popSize;
                        if(fold!=i)
                            break;
                        else
                            avg_fit += fitness;
                    }
                    avg_fit /= float(KFolds);
                    if(avg_fit!=0)
                        average.push_back({avg_fit,lambda,delta,popSize});
                }
                printMax(average,filename.substr(filename.find_last_of("/")+1),n_top);
            }else{
                float fold, fitness, lambda, avg_fit,popSize;
                while(!myfile.eof()){
                    avg_fit = 0;
                    for(unsigned int i=0;i<KFolds;i++){
                        myfile >> fold >> fitness >> lambda >> popSize;
                        if(fold!=i)
                            break;
                        else
                            avg_fit += fitness;
                    }
                    avg_fit /= float(KFolds);
                    if(avg_fit!=0)
                        average.push_back({avg_fit,lambda,0,popSize});
                }
                printMax(average,filename.substr(filename.find_last_of("/")+1),n_top);
            }
        }
        myfile.close();
    }
    return 0;
}
