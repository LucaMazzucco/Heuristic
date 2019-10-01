#include <iostream>
#include <random>
#include <algorithm>
#include "heuristic.h"
#include <cstring>
#include <time.h>

using namespace std;
Costs *temp;
Heuristic::Heuristic(string path){
    this->hasSolution = false;
    string line;
    string word;

    ifstream iffN(path.c_str());

    if (!iffN.is_open()) {
        cout << "Impossible to open" << path << endl;
        cin.get();
        exit(1);
    }

    getline(iffN, line);
    std::replace(line.begin(), line.end(), ';', ' ');
    istringstream iss(line);
    iss >> word;
    this->nCells = atoi(word.c_str());
    iss >> word;
    this->nTimeSteps = atoi(word.c_str());
    iss >> word;
    this->nCustomerTypes = atoi(word.c_str());

    // Memory allocation
    solution = new int***[nCells];
    problem.costs = new double***[nCells];
    for (int i = 0; i < this->nCells; i++) {
        problem.costs[i] = new double**[nCells];
        solution[i] = new int**[nCells];
        for (int j = 0; j < this->nCells; j++) {
            problem.costs[i][j] = new double*[nCustomerTypes];
            solution[i][j] = new int*[nCustomerTypes];
            for (int m = 0; m < this->nCustomerTypes; m++) {
                problem.costs[i][j][m] = new double[nTimeSteps];
                solution[i][j][m] = new int[nTimeSteps];
            }
        }
    }
    problem.n = new int[nCustomerTypes];
    problem.activities = new int[nCells];
    problem.usersCell = new int**[nCells];
    for (int i = 0; i < this->nCells; i++) {
        problem.usersCell[i] = new int*[nCustomerTypes];
        for (int m = 0; m < this->nCustomerTypes; m++) {
            problem.usersCell[i][m] = new int[nTimeSteps];
        }
    }

    getline(iffN, line);
    getline(iffN, line);
    std::replace(line.begin(), line.end(), ';', ' ');
    istringstream issN(line);
    for (int m = 0; m < nCustomerTypes; m++) {
        issN >> word;
        problem.n[m] = atoi(word.c_str());
    }

    getline(iffN, line);
    for (int m = 0; m < nCustomerTypes; m++) {
        for (int t = 0; t < nTimeSteps; t++) {
            getline(iffN, line);// linea con m e t
            for (int i = 0; i < nCells; i++) {
                getline(iffN, line);// linea della matrice c_{ij} per t ed m fissati
                istringstream issC(line);
                for (int j = 0; j < nCells; j++) {
                    issC >> word;
                    problem.costs[i][j][m][t] = atoi(word.c_str());
                }
            }
        }
    }

    getline(iffN, line);
    getline(iffN, line);
    std::replace(line.begin(), line.end(), ';', ' ');
    istringstream issA(line);
    for (int i = 0; i < nCells; i++) {
        issA >> word;
        problem.activities[i] = atoi(word.c_str());
    }

    getline(iffN, line);
    for (int m = 0; m < nCustomerTypes; m++) {
        for (int t = 0; t < nTimeSteps; t++) {
            getline(iffN, line);
            getline(iffN, line);
            std::replace(line.begin(), line.end(), ';', ' ');
            istringstream issU(line);
            for (int i = 0; i < nCells; i++) {
                issU >> word;
                problem.usersCell[i][m][t] = atoi(word.c_str());
            }
        }
    }
    for (int i = 0; i < nCells; i++)                     //inizializza la matrice a zero
            for (int j = 0; j < nCells; j++)                 //solution è la matrice degli user che sono richiesti dalla cella i alla cella j di tipo m al tempo t
                for (int m = 0; m < nCustomerTypes; m++)
                    for (int t = 0; t < nTimeSteps; t++)
                        solution[i][j][m][t] = 0;

}

void Heuristic::quicksort(Costs **A, int len)
{
  if (len < 2) return;

  double pivot = A[len / 2]->costo;

  int ii, jj;
  for (ii = 0, jj = len - 1; ; ii++, jj--)
  {
    while (A[ii]->costo < pivot) ii++;
    while (A[jj]->costo > pivot) jj--;

    if (ii >= jj) break;

    temp = A[ii];
    A[ii]     = A[jj];
    A[jj]     = temp;
  }

  quicksort(A, ii);
  quicksort(A + ii, len - ii);
}
void Heuristic::solveFast(vector<double>& stat, int timeLimit) {

    clock_t tStart = clock();
    double objFun=99999, objFunTemp=0, comb=0;
    int i, j, m, t, ordering[nCells], bestorder[nCells];
    int l, mark[nCells], ct=0, ff, rest = 0, n=0, lastcomb = 0,sol,flag = 0, betterflag, rand1, rand2;

    //FILE PER STAMPARE IL DETTAGLIO DELLE SCELTE
    //FILE *f;
    //f = fopen("output.txt", "w");
    int precision = 23;
    int N = precision*2*nCustomerTypes*nTimeSteps*4;

    COST *costs;

    struct modify{
        int val;
        int j;
        int m;
        int t;
        struct modify *next;
    }*mod, *modt;

    costs = (COST*)malloc(N*sizeof(COST));          //vettori di puntatori a struttura Costs


    for(i=0; i<N; i++)
        costs[i] = (COST)malloc(sizeof(struct Costs));

    for(comb = 0; lastcomb == 0; comb++){
                if((clock()-tStart)>4900){
                    lastcomb =1;
                }
                ct =0;
                ff = 0;

                mod = (struct modify*)malloc(sizeof(struct modify));
                mod-> next = NULL;
                modt = mod;
                if(lastcomb==1){
                    for(i=0; i<nCells; i++)
                        ordering[i] = bestorder[i];
                        objFun=0;
                }
                else{
                    if(flag==0){
                    for(l = 0; l<nCells; l++) mark[l] = 0;
                    while(ff<nCells&&ct<nCells*2){
                        l = rand()%nCells;
                        ct++;
                        if(mark[l]==0){
                            mark[l]=1;
                            ordering[ff] = l;
                            ff++;
                        }
                    }
                    for(l=0; l<nCells; l++){
                        if(mark[l]==0){
                            mark[l]=1;
                            ordering[ff] = l;
                            ff++;
                        }
                    }
               }

                else{

                    if(betterflag==0&&rand()>0.5){

                        l = ordering[rand1];
                        ordering[rand1] = ordering[rand2];
                        ordering[rand2] = l;

                    }
                    rand1 = rand()%nCells;
                    rand2 = rand()%nCells;

                    //printf("\n%d %d", rand1, rand2);
                    l = ordering[rand1];
                    ordering[rand1] = ordering[rand2];
                    ordering[rand2] = l;

                    betterflag = 0;
                }

                }

        for (n = 0; n <nCells; n++) {
            i = ordering[n];
            rest = 0;
            int demand = problem.activities[i], satisfaction=0;
            //fprintf(f, "\n n:%d, demand:%d;", n, demand);
            if(demand!=0){          //se la domanda è uguale a 0 salta tutto
                bool notSatisfied = true;
                int q = 0;


              for (m = 0; m < nCustomerTypes; m++) {                 //caricamento in vettore costs dei dati
                    for (t = 0; t < nTimeSteps; t++) {
                        for(int k = 1; k<precision; k++){
                            if((k+i)<nCells){
                            if(problem.usersCell[k+i][m][t]!=0){
                                costs[q]->costo = (double)problem.costs[k+i][i][m][t]/(double)problem.n[m];
                                costs[q]->j = k+i;
                                costs[q]->i = i;
                                costs[q]->m = m;
                                costs[q]->t = t;
                                q++;
                                satisfaction+=problem.n[m]*problem.usersCell[k+i][m][t];
                            }}}
                        for(int k = 1; k<precision; k++){
                            if((-k+i)>=0){
                            if(problem.usersCell[i-k][m][t]!=0){
                                costs[q]->costo = (double)problem.costs[i-k][i][m][t]/(double)problem.n[m];
                                costs[q]->j = i-k;
                                costs[q]->i = i;
                                costs[q]->m = m;
                                costs[q]->t = t;
                                q++;
                                satisfaction+=problem.n[m]*problem.usersCell[i-k][m][t];
                            }}
                            }
                    }
              }


                for(int k = precision; k<nCells&&satisfaction<demand; k++){

                     for (m = 0; m < nCustomerTypes; m++) {                 //caricamento in vettore costs dei dati
                            for (t = 0; t < nTimeSteps; t++) {
                                if((k+i)<nCells){
                                if(problem.usersCell[k+i][m][t]!=0){
                                    costs[q]->costo = (double)problem.costs[k+i][i][m][t]/(double)problem.n[m];
                                    costs[q]->j = k+i;
                                    costs[q]->i = i;
                                    costs[q]->m = m;
                                    costs[q]->t = t;
                                    q++;
                                    satisfaction+=problem.n[m]*problem.usersCell[k+i][m][t];
                                }}
                                if((-k+i)>=0){
                                if(problem.usersCell[i-k][m][t]!=0){
                                    costs[q]->costo = (double)problem.costs[i-k][i][m][t]/(double)problem.n[m];
                                    costs[q]->j = i-k;
                                    costs[q]->i = i;
                                    costs[q]->m = m;
                                    costs[q]->t = t;
                                    q++;
                                    satisfaction+=problem.n[m]*problem.usersCell[i-k][m][t];
                                }}
                        }
                    }
                }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                quicksort(costs,q);
                //quicksort funzione a riga circa 112 di questo file
                //ho dovuto spostare la tua struttura per renderla leggibile alla funzione
                //ricorsiva nel file h. io non so se è corretto, ma funziona. controlla :)
                //essendo stata spostata al di fuori, magari la struttura non necessita di essere
                //dichiarata in tutte le funzioni.
               /* for(int w=0; w<q-1; w++){                   //ordinamento vettore costs in base al costo proporzionato al numero di tasks
                    for(int e=w+1; e<q; e++){
                        if(costs[w]->costo >= costs[e]->costo){

                            if(costs[w]->costo == costs[e]->costo){

                               if(rand()>0.5){
                                    temp=costs[e];
                                    costs[e] = costs[w];
                                    costs[w] = temp;
                               }
                            }
                            else{

                                temp=costs[e];
                                costs[e] = costs[w];
                                costs[w] = temp;
                             }
                        }
                    }
                }*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                for(int w = 0; w < q && notSatisfied;  w++){        //scelta effettiva della soluzione prendendo gli elementi in ordine da costs
                                j = costs[w]->j;
                                m = costs[w]->m;
                                t = costs[w]->t;
                                sol = 0;
                                if (demand > problem.n[m] * problem.usersCell[j][m][t]) {   //se demand è maggiore di numero di users *numero di task che possono soddisfare gli user
                                                                                  //carica in solution tutti gli utenti disponibili
                                    sol = problem.usersCell[j][m][t];
                                    problem.usersCell[j][m][t] =0;     //elimina tutti gli utenti da quel percorso

                                }
                                else {

                                    sol = floor((double)demand / (double)problem.n[m]);

                                    problem.usersCell[j][m][t] -= sol;

                                }

                                if (sol != 0){
                                    if(lastcomb==1){
                                        solution[i][j][m][t] = sol;
                                        objFun += sol * problem.costs[j][i][m][t];

                                    }
                                    objFunTemp += sol * problem.costs[j][i][m][t]; //aggiunge costo corrente
                                    demand -= problem.n[m]*sol;

                                    modt->val = sol;
                                    modt->j = j;
                                    modt->m = m;
                                    modt->t = t;
                                    modt->next = (struct modify*)malloc(sizeof(struct modify));
                                    modt = modt->next;
                                    modt->next = NULL;
                                    if(objFunTemp>objFun) goto nextloop;
                                }
                                if(demand<=0) notSatisfied = false;
                            }

                if(notSatisfied||demand>0) {goto nextloop;}//SE NON SODDISFA LA DOMANDA SALTA ALLA PROSSIMA COMBINAZIONE
            }
        }

        if(objFunTemp<objFun&&n==nCells){
            //if(flag==1){LIBERO BEST PRECEDENTE}
            for(i=0; i<nCells; i++)
                bestorder[i] = ordering[i];

            objFun = objFunTemp;
            flag = 1;
            betterflag = 1;
        }
nextloop:
        modt = mod;
        while(modt->next != NULL){

            problem.usersCell[modt->j][modt->m][modt->t] += modt->val;
            mod = modt;
            modt = modt->next;
            free(mod);
        }


        objFunTemp = 0;

    }

     printf("\nbest sol found:%.f", objFun);


    printf("\nCombination number: %.f", comb);

    stat.push_back(objFun);
    stat.push_back((double)(clock() - tStart) / CLOCKS_PER_SEC);
    //fclose(f);
    hasSolution=true;
}



void Heuristic::writeKPI(string path, string instanceName, vector<double> stat){
    if (!hasSolution)
        return;

    ofstream fileO(path, ios::app);
    if(!fileO.is_open())
        return;

    fileO << instanceName << ";" << stat[0] << ";" << stat[1];
    for(int i=2; i<stat.size(); i++)
        fileO <<  ";" << stat[i];
    fileO << endl;

    fileO.close();

}

void Heuristic::writeSolution(string path) {
    if (!hasSolution)
        return;

    ofstream fileO(path);
    if(!fileO.is_open())
        return;

    fileO << this->nCells << "; " << this->nTimeSteps << "; " << this->nCustomerTypes << endl;
    for (int m = 0; m < this->nCustomerTypes; m++)
        for (int t = 0; t < this->nTimeSteps; t++)
            for (int i = 0; i < this->nCells; i++)
                for (int j = 0; j < this->nCells; j++)
                    if (solution[i][j][m][t] > 0)
                        fileO << i << ";" << j << ";" << m << ";" << t << ";" << solution[i][j][m][t] << endl;

    fileO.close();
}

eFeasibleState Heuristic::isFeasible(string path) {

    string line;
    string word;
    int nCellsN;
    int nTimeStepsN;
    int nCustomerTypesN;
    int i, j, m, t;


    ifstream iffN(path.c_str());

    if (!iffN.is_open()) {
        cout << "Impossible to open" << path << endl;
        exit(1);
    }

    getline(iffN, line);
    std::replace(line.begin(), line.end(), ';', ' ');
    istringstream iss(line);
    iss >> word; // nCells
    nCellsN = atoi(word.c_str());
    iss >> word; // nTimeSteps
    nTimeStepsN = atoi(word.c_str());
    iss >> word; // nCustomerTypes
    nCustomerTypesN = atoi(word.c_str());

    int**** solutionN = new int***[nCells];
    for (i = 0; i < nCellsN; i++) {
        solutionN[i] = new int**[nCells];
        for (j = 0; j < nCellsN; j++) {
            solutionN[i][j] = new int*[nCustomerTypes];
            for (m = 0; m < nCustomerTypesN; m++) {
                solutionN[i][j][m] = new int[nTimeSteps];
                for ( t = 0; t < nTimeStepsN; t++) {
                    solutionN[i][j][m][t] = 0;
                }
            }
        }
    }

    while (getline(iffN, line)) {
        std::replace(line.begin(), line.end(), ';', ' ');
        istringstream iss(line);
        iss >> word; // i
        i = atoi(word.c_str());
        iss >> word; // j
        j = atoi(word.c_str());
        iss >> word; // m
        m = atoi(word.c_str());
        iss >> word; // t
        t = atoi(word.c_str());
        iss >> word; // value
        solutionN[i][j][m][t] = atoi(word.c_str());
    }

    // Demand
    bool feasible = true;
    int expr;
    for (int i = 0; i < nCells; i++) {
        expr = 0;
        for (int j = 0; j < nCells; j++)
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++)
                    expr += problem.n[m] * solutionN[j][i][m][t];
        if (expr < problem.activities[i])
            feasible = false;
    }

    if (!feasible)
        return NOT_FEASIBLE_DEMAND;

    // Max Number of users
    for (int i = 0; i < nCells; i++)
        for (int m = 0; m < nCustomerTypes; m++)
            for (int t = 0; t < nTimeSteps; t++) {
                expr = 0;
                for (int j = 0; j < nCells; j++)
                    expr += solutionN[i][j][m][t];
                if (expr > problem.usersCell[i][m][t])
                    feasible = false;
            }

    if(!feasible)
        return NOT_FEASIBLE_USERS;

    return FEASIBLE;
}

void Heuristic::getStatSolution(vector<double>& stat) {
    if (!hasSolution)
        return;

    int* tipi = new int[nCustomerTypes];
    for (int m = 0; m < nCustomerTypes; m++)
        tipi[m] = 0;

    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++)
            for (int t = 0; t < nTimeSteps; t++)
                for (int m = 0; m < nCustomerTypes; m++)
                    if (solution[i][j][m][t] > 0)
                        tipi[m] += solution[i][j][m][t];
    for (int m = 0; m < nCustomerTypes; m++)
        stat.push_back(tipi[m]);

}
