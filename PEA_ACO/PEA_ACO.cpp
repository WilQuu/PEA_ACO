#include <iostream>
#include <fstream>
#include<string>
#include<cstdlib>
#include <sstream>
#include <cstring>
#include<time.h>
#include<math.h>
#include<ctime>
#include<vector>
#include <algorithm>
#include<map>

using namespace std;

vector<unsigned> random_path(int citiesNum) {
    vector<unsigned> random_path;
    vector<int> tempNums;
    for (int i = 1; i < citiesNum; i++)
        tempNums.push_back(i);

    random_path.push_back(0);
    int currentVertex;
    for (int i = 1; i < citiesNum; i++) {
        currentVertex = rand() % tempNums.size();
        vector<int>::iterator it;
        it = tempNums.begin() + currentVertex;
        random_path.push_back(tempNums.at(currentVertex));
        tempNums.erase(it);

    }
    random_path.push_back(0);
    return random_path;
}
void print_path(vector<unsigned> path) {

    for (int i = 0; i < path.size(); i++) {
        if (i == path.size() - 1) {
            cout << path.at(i) << endl;
        }
        else
            cout << path.at(i) << " -> ";

    }
}

int count_cost(vector<unsigned> path, int** cityMatrix) {
    int cost = 0;
    for (int i = 0; i < path.size() - 1; i++) {
        cost += cityMatrix[path.at(i)][path.at(i + 1)];
    }

    cost += cityMatrix[path.at(path.size() - 1)][path.at(0)];
    return cost;
}

bool contain(vector<unsigned> vec, int value, int start, int end) {
    for (int i = start; i <= end; i++) {
        if (vec.at(i) == value) {
            return true;
        }
    }
    return false;
}



bool compare_paths(vector<unsigned> vec_a, vector<unsigned> vec_b) {
    for (int i = 0; i < vec_a.size(); i++) {
        if (vec_a.at(i) != vec_b.at(i))
            return false;
    }

    return true;
}

int get_nn_value(int citiesNum,int** cityMatrix) {
    bool* visitedTab = new bool[citiesNum];

    for (int i = 0; i < citiesNum; i++) {
        visitedTab[i] = false;
    }
    int hamiltonGreedyVersion = 0;
    int hamiltonMostOptFound = INT_MAX;
    int currentX, nextX; // indeksy do sprawdzania pozycji rozwazanego wierzcholka
    int startIndex;
    for (int y = 0; y < citiesNum; y++) {
        startIndex = y;
        visitedTab[startIndex] = true;
        currentX = startIndex;
        bool isThereMoreCities = true;
        int howManyCities = 0;

        while (isThereMoreCities) {
            int minEdge = INT_MAX;

            for (int i = 0; i < citiesNum; i++) {
                if (cityMatrix[currentX][i] < minEdge
                    && !visitedTab[i]
                    && cityMatrix[currentX][i] != 0) {

                    minEdge = cityMatrix[currentX][i];
                    nextX = i;
                    visitedTab[i] = true;
                }
            }

            hamiltonGreedyVersion += minEdge;
            currentX = nextX;
            isThereMoreCities = false;
            for (int i = 0; i < citiesNum; i++) {
                if (visitedTab[i] == false) {
                    isThereMoreCities = true;
                }
            }
            howManyCities++;
            hamiltonGreedyVersion += cityMatrix[currentX][startIndex];
        }
       
        howManyCities = 0;
        nextX = NULL;
        for (int z = 0; z < citiesNum; z++) {
            visitedTab[z] = false;
        }
        if (hamiltonMostOptFound > hamiltonGreedyVersion) {
            hamiltonMostOptFound = hamiltonGreedyVersion;
        }
        hamiltonGreedyVersion = 0;
    }
    delete[]visitedTab;
    return hamiltonMostOptFound;
}


bool time_stop(clock_t start_time, clock_t now_time, double max_time) {
    double  elapsed = (double)(now_time - start_time) / CLOCKS_PER_SEC;
    if (max_time > elapsed)
        return true;
    else
        return false;
}

struct Ant {
    vector<bool>visited;
    int starting_vertex;
    vector<unsigned>current_path;
};
bool visited_all(Ant* ant) {
    for (int i = 0; i < ant->visited.size(); i++) {
        if (!ant->visited[i]) {
            return false;
        }
    }
    return true;
}
void count_probability(Ant* ant, map<double, int>& probability_map, int& last_vertex, double** pheromones, double alpha, int citiesNum, double beta, int** cityMatrix) {
    double sum = 0;
    int current_vertex = ant->current_path[ant->current_path.size() - 1];
    for(int i=0;i<citiesNum;i++){
        if (ant->visited[i] == 0) {
            double licznik = pow(pheromones[current_vertex][i], alpha) / pow((cityMatrix[current_vertex][i]), beta);
            for (int j = 0; j < citiesNum; j++) {
                if (ant->visited[j] == false && i != j) {
                    sum += pow(pheromones[current_vertex][j], alpha) / pow((cityMatrix[current_vertex][j]), beta);
                }
            } 
            if (sum != 0) {
                double value = licznik / sum;
                probability_map[value] = i;
            }
            else {
                last_vertex = i;
                return;
            }
        }
    }
}



void update_pheromones(double** pheromones, int citiesNum,double rho) {
    for (int i = 0; i < citiesNum; i++) {
        for (int j = 0; j < citiesNum; j++) {
            pheromones[i][j] = pheromones[i][j] * rho;
        }
    }
}

int next_city(map<double, int> probability_map) {
    double rand_value = ((double)rand() / (RAND_MAX));
    double probability_value = 0;
    double omega = 0;
    double percentage;
    map<double,int >::iterator it;
    for (it = probability_map.begin(); it != probability_map.end(); it++) {
        omega += it->first;
    }
    for (it = probability_map.begin(); it != probability_map.end(); it++) {
        probability_value += it->first;
        percentage = probability_value / omega;
        if (percentage > rand_value) {
            return it->second;
        }
    }
}
int main() {
    srand(time(NULL));
    fstream file;
    string tempLine;
    int hamiltonOpt;
    int citiesNum;
    string fileName;
    cout << "Name of instance : ";
    cin >> fileName;
    file.open(fileName, std::ios::in | std::ios::out);
    if (file.is_open()) {
        cout << "Max time value [s] :";
        double max_time;
        cin >> max_time;
        cout << "starting algorithm ..." << endl;
        getline(file, tempLine);
        getline(file, tempLine);
        citiesNum = atoi(tempLine.c_str());
        int** cityMatrix = new int* [citiesNum];
        for (int i = 0; i < citiesNum; i++)
            cityMatrix[i] = new int[citiesNum];


        for (int i = 0; i < citiesNum; i++) {
            getline(file, tempLine);
            stringstream ss;
            ss << tempLine;
            int found;
            string temp;
            int y = 0;
            while (!ss.eof()) {
                ss >> temp;
                if (stringstream(temp) >> found) {
                    cityMatrix[i][y] = found;
                    y++;
                }
            }

            temp = "";
        }
        getline(file, tempLine);
        hamiltonOpt = atoi(tempLine.c_str());
            int phermone_rate;
            /* deklaracja zmiennych */
            double** pheromone = new double* [citiesNum];

            for (int i = 0; i < citiesNum; i++)
                pheromone[i] = new double[citiesNum];
            int cost = 0;
            vector<unsigned> temp_path;
            temp_path = random_path(citiesNum);
            double start_pheromone_value = (double)citiesNum/(count_cost(temp_path,cityMatrix)/5);
            cout << start_pheromone_value << endl;
            for (int i = 0; i < citiesNum; i++) {
                for (int j = 0; j < citiesNum; j++) {
                    pheromone[i][j] = start_pheromone_value;
                }
            }

            int alpha = 1;
            int beta = 2;
            double rho = 0.5;
            vector<unsigned> best_path;
            vector<unsigned> opt_road;
            best_path = random_path(citiesNum);

            clock_t start = clock();
            clock_t now_time;
            vector<unsigned> temp_best_path = random_path(citiesNum);
            double PRDs = (double)count_cost(temp_best_path, cityMatrix) / hamiltonOpt * 100;
            cout << " ant - " << PRDs << " %" << endl;
            int temp_best_value = count_cost(temp_best_path, cityMatrix);
            /*      poczatek algorytmu ACO    */
            int i = 0;
            do
            {
                /* petla mrowek */
                for (int j = 0; j < citiesNum; j++) {
                    Ant* ant = new Ant();
                    ant->visited.resize(citiesNum, false);
                    ant->current_path.push_back(rand() % citiesNum);
                    ant->visited[ant->current_path.at(0)] = true;
                    //int last_vertex = -1;
              
                    while (ant->current_path.size() <= citiesNum) {
                        //choose next city
                        int last_vertex = -1;
                        map<double, int> probability_map;
                        count_probability(ant, probability_map, last_vertex, pheromone, alpha, citiesNum, beta, cityMatrix);
                 
                        if (last_vertex != -1) {
                            ant->current_path.push_back(last_vertex);
                            ant->visited[last_vertex] = true;
                            //dodanie wierzchołka końcowego (pierwszego) do trasy
                            ant->current_path.push_back(ant->current_path.at(0));

                            pheromone[ant->current_path[ant->current_path.size() - 3]][last_vertex] += (double)citiesNum;
                            pheromone[last_vertex][ant->current_path.at(0)] += (double)citiesNum;



                        }
                        else {
                            int next_vertex = next_city(probability_map);
                            // rozpylenie feromonów
                            pheromone[ant->current_path[ant->current_path.size() - 1]][next_vertex] += (double)citiesNum;
                            // dodanie wierzchołka do drogi mrówki
                            ant->current_path.push_back(next_vertex);
                            // ustawienie go jako odwiedzonego 
                            ant->visited[next_vertex] = true;
                        }

                    }
                    ant->current_path.at(ant->current_path.size() - 1) = ant->current_path.at(0);
                    if (temp_best_value > count_cost(ant->current_path, cityMatrix) ) {
                        temp_best_path.clear();
                        temp_best_path = ant->current_path;
                        double PRD = (double)count_cost(temp_best_path, cityMatrix) / hamiltonOpt * 100;
                        cout << i << " wave by " << j << " ant - " << PRD << " %" << endl;
                        temp_best_value = count_cost(temp_best_path, cityMatrix);
                    }
                    delete ant;
                    now_time = clock(); 
                    if (!time_stop(start, now_time, max_time)) {
                        double  elapsed = (double)(now_time - start) / CLOCKS_PER_SEC;
                        cout << elapsed << " [s]" << endl;
                        double PRD = (double)count_cost(temp_best_path, cityMatrix) / hamiltonOpt * 100;
                        cout << PRD << " % - prd" << endl;
                        print_path(temp_best_path);
                        return 0;
                    }
                
                }
                update_pheromones(pheromone, citiesNum, rho);
                i++;

                
                now_time = clock();
            } while (time_stop(start, now_time, max_time)&&i<200);

             /*      koniec algorytmu          */
            clock_t stop = clock();

            double  elapsed = (double)(stop - start) / CLOCKS_PER_SEC;
            cout << elapsed << " [s]" << endl;
            double PRD = (double)count_cost(temp_best_path, cityMatrix) / hamiltonOpt * 100;
            cout << PRD << " % - prd" << endl;
            print_path(temp_best_path);

            /*        zwalnianie pamieci               */
            for (int i = 0; i < citiesNum; i++)
                delete[] cityMatrix[i];
            delete[]cityMatrix;

            for (int i = 0; i < citiesNum; i++)
                delete[] pheromone[i];
            delete[]pheromone;
        
        cout << "end" << endl;
    }
    else cerr << "Dostep do pliku zostal zabroniony!" << endl;



    file.close();
    return 0;
}
