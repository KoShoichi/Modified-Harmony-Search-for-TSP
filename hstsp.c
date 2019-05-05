//
//  hstsp.c
//  Modified Harmony Search For TSP
//
//  Created by 顾肖一 on 2017/05/17.
//  Copyright © 2017年 Bill Gu. All rights reserved.
//
//  hstsp.c: Defines the entry point for the console application.
/***Harmony Search (HS) is a novel intelligent algorithm firstly designed by Geem, Kim and Loganathan (2001). This is a solver based on HS and it is modified by introducing greedy algorithm for early initialization, which can be proved to obtain a much faster convergence speed. In order to ensure the intensification and diversification of Harmony Memory, and to refrain from being tendentious towards similar solutions during early period, it makes the starting-point of greedy paths distinctly. Decent combination of local search has also be utilized to avoid the uncertainty caused by randomness as far as possible.
 ***/
//  For academic use, please contact: xiaoyig@student.unimelb.edu.au

#include "StdAfx.h"
#include "hstsp.h"
#include "cpu_time.h"
#define MAXCITY 30000
#define MAX_STR 1024
#define HMS     10              // Harmony memory Size, HMS
#define HMCR    0.95              // Harmony memory considering rate, HMCR
#define PAR     0.2		   // Pitch adjusting rate, PAR
#define BW      0.01             // Default BandWidth

//Limits and buffer sizes:
//#define LineMAX 128
#define WORD_MAX 64

//Control values for anneal/hybrid algorithm:
#define SATISFIED 10000
#define START_TEMP (avg_distance/40.0)
#define END_TEMP (avg_distance/10.0)

//Control values for only anneal:
#define DELTA_TEMP (0.9999)

//Distance Calculation HS
#define dist(i,j)   ((int)( sqrt(  (gHarmonyMemory[i][j].X-gHarmonyMemory[i][j-1].X)*(gHarmonyMemory[i][j].X-gHarmonyMemory[i][j-1].X) + (gHarmonyMemory[i][j].Y-gHarmonyMemory[i][j-1].Y)* (gHarmonyMemory[i][j].Y-gHarmonyMemory[i][j-1].Y) ) + 0.5))
#define distNew(j)  ((int)( sqrt(  (newHarmonyCity[j].X-newHarmonyCity[j-1].X)*(newHarmonyCity[j].X-newHarmonyCity[j-1].X) + (newHarmonyCity[j].Y-newHarmonyCity[j-1].Y)* (newHarmonyCity[j].Y-newHarmonyCity[j-1].Y) ) + 0.5))


typedef struct{
    int row;
    int col;
} TSPINDEX;

//Define Golbal Variables
int CITYS = 0;
double superHMCR = HMCR;
double superPAR  = PAR;
double superBW   = BW;

TSPCITY gHarmonyMemory[HMS+1][MAXCITY]; //Define all the Harmony Memory,the last column is for the result;
TSPCITY newHarmonyCity[MAXCITY];      //Generator a new Solution;
TSPCITY tspBest = {0,0,0,0};
TSPCITY tspWorst= {0,0,0,0};

//FUNCTION PROTOTYPES for Harmony Search
void UpdateHMMemory(int i,int j,TSPCITY cityInfo);
void GeneratorNewHarmonyMemory(TSPdata *tspdata);
TSPCITY RandomCityGeneration(TSPdata *tspdata, int i);
void HMParameterSetting(double ConsideringRate, double AdjustingRate, double bandWidth);

//FUNCTION PROTOTYPES for Adjustment Algorithm:

void get_options(int argc, char **argv);
void set_best(int distance, int *path);
void swap(int i, int j, int *array);
void copy_array(int *to, int *from, int len);
int get_list_of_cities(int *list);
void print_distances();
void print_distance(int i, int j);
void print_cities();
void print_city(city *c);
void print_solution();
void calc_distances(int max_id);
int calc_distance(city *a, city *b);
int get_distance(int i, int j);
int calc_path_dist(int *path, int len);
void free_distances();
void nearest_neighbor(int *path, int len);
int swap_closest(int *remaining, int num_remaining);
void two_opt(int *path, int len);
void hybrid(int timelimit, int *path, int len, Vdata *vdata);
void anneal(int *path, int len);
int anneal_accept(int new_dst, int old_dst, double temp);
double change_temp(double old_temp);
void two_opt_swap(int i, int j, int *path);
int two_opt_dist(int old_dist, int i, int j, int *path, int len);
void sig_handler(int sig);
void install_sig_handlers();
double get_max(double a, double b);
int DataPass(int linenumber,Vdata *vdata,TSPdata *tspdata);
void executionGreedy(TSPdata *tspdata, Vdata *vdata, int fuko);
void Greedy(int *path, int len, TSPdata *tspdata);
void ReturnBestSolution(Vdata *vdata,TSPdata *tspdata);
void DataPassNew(Vdata *vdata, TSPdata *tspdata);
void GreedyGeneration(TSPdata *tspdata, Vdata *vdata, int timelimit);
void InitReversion(TSPdata *tspdata, Vdata *vdata);

//Static Variables, this will not be used crossing files

//The list of cities and their coordinates:
static city  *cities[MAXCITY];
static int num_cities;

//The matrix of distances between cities
//and the average distance between cities:
static int ** distances;
static int avg_distance;

//The optimal distance/path found thus far
//(printed on a SIGTERM or SIGINT):
static int best_distance;
static int * best_path;

//The command line options chosen:
static int use_anneal = 0;
static int use_nearest_neighbor = 0;
static int use_two_opt = 0;
static int verbose = 0;
static int debug = 0;

//The input and output options/filenames:
static int in_file = 0;
static char in_filename[WORD_MAX];
static int out_file = 0;
static char out_filename[WORD_MAX];


/**@Function:getMinNode
 * @Params: TSPdata *tspdata
 * @comments:get the Min_Node_No from tsp file
 * @author:  Bill Gu
 */
void getMinNode(TSPdata *tspdata)
{
    CITYS = tspdata->min_node_num;
}

/**@Function:tspWorstSetting
 * @comments:sets the line number of the worst solution in HM
 * @author:  Bill Gu
 */
void tspWorstSetting()
{
    tspWorst.CityIndex = 0;
    tspWorst.Fitness   = gHarmonyMemory[0][CITYS].Fitness;
    for(int i=1;i<HMS;i++)
    {
        if(gHarmonyMemory[i][CITYS].Fitness>tspWorst.Fitness)
        {
            tspWorst.CityIndex = i;
            tspWorst.Fitness   = gHarmonyMemory[i][CITYS].Fitness;
        }
    }
}

/**@Function:tspBestSetting
 * @comments:sets the line number of the best solution in HM after all iteration
 * @author:  Bill Gu
 */
void tspBestSetting()
{
    tspBest.CityIndex = 0;
    tspBest.Fitness   = gHarmonyMemory[0][CITYS].Fitness;
    for(int i=1;i<HMS;i++)
    {
        if(gHarmonyMemory[i][CITYS].Fitness<tspBest.Fitness)
        {
            tspBest.CityIndex = i;
            tspBest.Fitness   = gHarmonyMemory[i][CITYS].Fitness;
        }
    }
}
/**@Function:RandomCityOfIndex
 * @Params:  int Range, define Range;
 * @return:  int get index of city;
 * @comments: gets Random, generates a vaule between 0 to Range -1
 * @author:  Bill Gu
 */
int RandomCityOfIndex(int Range)
{
    return (rand()%Range);
}

/**@Function:RandomDoubleNext
 * @Params:  void
 * @return:  float between 0.00-1.00;
 * @comments: gets Random, generates a float between 0.00 to 1.00
 * @author:  Bill Gu
 */
double RandomDoubleNext()
{
    return (double)(rand()/(double)RAND_MAX);
}

/**@Function:CheckCityIndexExist
 * @Params:  int solIndex  choose one solution
 * @Params:  int index  choose one CityIndex
 * @Params:  int currentPosition, current plot has been filled up with cities in the array
 * @return:  True (not 0), the City exists, False(0), the City does not exist
 * @comments: check if the Current City exists or not in Current Solution.
 * @author:  Bill Gu
 */
int CheckCityIndexExist(int solIndex,int index,int currentPosition)
{
    int i;
    if (solIndex == -1)  //new Harmony Check when solIndex == -1, it is the new Harmony
    {
        for ( i=0;i<currentPosition;i++)
        {
            if (index == newHarmonyCity[i].CityIndex)
            {
                return 1; //this index is existing/selected
            }
        }
        return 0;
    }
    else
    {
        for (i =0;i<currentPosition;i++)
        {
            if (index == gHarmonyMemory[solIndex][i].CityIndex)
            {
                return 1; //this index is existing/selected
            }
        }
        return 0;
    }
}

/**@Function:memoryConsideration
 * @Params:  int index, get the Memory From HM record
 * @return:  TSPINDEX, get the new index of City column and row
 * @comments: considers picking up a new city from the HM record
 * @author:  Bill Gu
 */
TSPINDEX  memoryConsideration(int index)
{
    TSPINDEX tspIndex = {-1,-1};
    int HmIndex = RandomCityOfIndex(HMS);
    tspIndex.row = HmIndex;
    tspIndex.col = index;
    return tspIndex;
}

/**@Function:pitchAdjustment
 * @Params:  double bw, bandwidth to constrain the near cities chosen (how far from current city)
 * @Params:  int index: Get the Memory From HM record ;
 * @Params:  TSPdata *tspdata, get nearest city form tspdata according to BW parameter
 * @return:  TSPCITY, get the new index of City column and row
 * @comments: BW define the BandWidth of current city, to adjustment form near city range;
 * @author:  Bill Gu
 */
TSPCITY pitchAdjustment(double bw, int index,TSPdata *tspdata)
{
    //to do some adjustment
    int i;
    int NearCitys =0;
    int RandomNear;
    double distance;
    
    TSPCITY* adjustHarmonycity=(TSPCITY*)malloc(tspdata->n*sizeof(TSPCITY));
    TSPCITY tsp =  newHarmonyCity[index];
    TSPCITY tspTemp={0,0,0,0};
    double randRange = tspWorst.Fitness*bw;  //define range of distance to adjust city,get ajusted city nearby;
    
    for (i=0;i<tspdata->n;i++)
    {
        distance=sqrt((tsp.X-tspdata->x[i])*(tsp.X-tspdata->x[i]) + (tsp.Y-tspdata->y[i])*(tsp.Y-tspdata->y[i])) +0.5;//get current distance from tsp to tspdata[i]
        if (distance<randRange&&distance!=0) //it is Ok to do ajustment
        {
            tspTemp.CityIndex = i;
            tspTemp.X         = tspdata->x[i];
            tspTemp.Y         = tspdata->y[i];
            adjustHarmonycity[NearCitys] =  tspTemp;  //build adjustment array Harmony Citys;
            NearCitys++;
        }
    }
    
    if (NearCitys==0)
    {
        return tsp;  //no just for current ,because of no near city contains;
    }
    
    RandomNear= RandomCityOfIndex(NearCitys);
    free(adjustHarmonycity);
    return adjustHarmonycity[RandomNear];
}

/**@Function:PrintHarmoryMemory
 * @comments:prints Current solution of Harmony Memory
 * @Author:  Bill Gu
 */
void PrintHarmoryMemory()
{
    int i,j;
    for (i=0;i<HMS;i++)
    {
        printf("HarmoyMemory[%d]\n",i+1);
        for (j=0;j<CITYS;j++)
        {
            //printf("CityIndex[%d],LocationX[%f],LocationY[%f]==>",gHarmonyMemory[i][j].CityIndex+1,gHarmonyMemory[i][j].X,gHarmonyMemory[i][j].Y);
            printf("%d ",gHarmonyMemory[i][j].CityIndex+1);
        }
        //printf("CityIndex[%d], LocationX[%f], LocationY[%f]",gHarmonyMemory[i][0].CityIndex+1,gHarmonyMemory[i][0].X,gHarmonyMemory[i][0].Y);
        printf("%d \n",gHarmonyMemory[i][0].CityIndex+1);
        printf("total Distance = %d\n", gHarmonyMemory[i][CITYS].Fitness);
    }
    printf("the best  solution is: %d\n",tspBest.CityIndex+1);
    printf("the worst solution is: %d\n",tspWorst.CityIndex+1);
}

/**@Function: PrintOnlyBestAndWorst
 * @comments: prints the best and worst solution in HM
 * @Author:  Bill Gu
 */
void PrintOnlyBestAndWorst()
{
    printf("The Best Harmony[%d]\n",tspBest.CityIndex+1);
    for (int j=0;j<CITYS;j++)
    {
        printf("%d ",gHarmonyMemory[tspBest.CityIndex][j].CityIndex+1);
    }
    printf("%d \n",gHarmonyMemory[tspBest.CityIndex][0].CityIndex+1);
    printf("total Distance = %d\n", gHarmonyMemory[tspBest.CityIndex][CITYS].Fitness);
    printf("The Worst Harmony[%d]\n",tspWorst.CityIndex+1);
    for (int j=0;j<CITYS;j++)
    {
        printf("%d ",gHarmonyMemory[tspWorst.CityIndex][j].CityIndex+1);
    }
    printf("%d \n",gHarmonyMemory[tspWorst.CityIndex][0].CityIndex+1);
    printf("total Distance = %d\n", gHarmonyMemory[tspWorst.CityIndex][CITYS].Fitness);
}


/**@Function: ReturnBestHSSolution
 * @Params: Vdata *vdata
 * @Params: TSPdata *tspdata
 * @comments: return the best solution after calculated by Harmony Search Algorithm and pass this to bestsol vector of vdata
 * @author:  Bill Gu
 */
void ReturnBestHSSolution(Vdata *vdata,TSPdata *tspdata)
{
    for(int i=0;i<CITYS;i++)
    {
        vdata->bestsol[i] = gHarmonyMemory[tspBest.CityIndex][i].CityIndex;
    }
    for(int tmp=CITYS;tmp<tspdata->n;tmp++)
    {
        vdata->bestsol[tmp] = -1;
    }
    
}

/**@Function: updateGeneralFitness
 * @Params:   TSPCITY fitnessCity, Current fitnessCity
 * @comments: stores the Current best or worst City Distance
 * @author:  Bill Gu
 */
void updateGeneralFitness(TSPCITY fitnessCity)
{
    if (tspBest.Fitness == 0)
    {
        tspBest = fitnessCity;
    }
    
    if (tspWorst.Fitness == 0)
    {
        tspWorst = fitnessCity;
    }
    
    //find fitnessCity shorter than tspBest then implement
    if (fitnessCity.Fitness<tspBest.Fitness)
    {
        tspBest = fitnessCity;
    }
    
    //find fitnessCity longer than tspWorst then implement
    if (fitnessCity.Fitness>tspWorst.Fitness)
    {
        tspWorst = fitnessCity;
    }
}


/**@Function:updateFitnessHistory
 * @comments:updates the Worst Tsp Form the latest Harmony Memory, update the fitness information
 * @author:  Bill Gu
 */
void updateFitnessHistory()
{
    int i;
    TSPCITY tmpworst = {0,0,0,0};
    for (i=0;i<HMS;i++)
    {
        //if current fitness is longer than worst than do this
        if (gHarmonyMemory[i][CITYS].Fitness>tspWorst.Fitness)
        {
            tmpworst.CityIndex = i;
            tmpworst.Fitness = gHarmonyMemory[i][CITYS].Fitness;
            tspWorst = tmpworst;
        }
    }
}

/**@Function:calculateFitnessForRoute
 * @comments:updates Global Harmony Memory,store the Best result,this will be calculated through each iteration when a new path has been found
 * @author:  Bill Gu
 */
void calculateFitnessForRoute()
{
    int          i =  0;
    int      dBest =  0;  //store the distance of newHarmony;
    
    //calculate NewHarmony Distance
    TSPCITY  tspCityResult = {0,0,0,-1};
    for (i=1; i<CITYS; i++)
    {
        dBest += distNew(i);
    }
    dBest += ((int)(sqrt((newHarmonyCity[CITYS-1].X-newHarmonyCity[0].X)*(newHarmonyCity[CITYS-1].X-newHarmonyCity[0].X) + (newHarmonyCity[CITYS-1].Y-newHarmonyCity[0].Y)* (newHarmonyCity[CITYS-1].Y-newHarmonyCity[0].Y) ) + 0.5));
    tspCityResult.Fitness = dBest;
    UpdateHMMemory(-1,CITYS,tspCityResult);
    //printf("the worst solution[%d]",tspWorst.CityIndex);
    //printf("the best  solution[%d]",tspBest.CityIndex);
    //check Current Distance is shorter than the worst ;
    if (dBest < tspWorst.Fitness)
    {
        //replace the worst solution
        for (i=0;i<CITYS;i++)
        {
            gHarmonyMemory[tspWorst.CityIndex][i] = newHarmonyCity[i];
        }
        tspWorst.Fitness = dBest;
        UpdateHMMemory(tspWorst.CityIndex,CITYS,newHarmonyCity[CITYS]);
        updateFitnessHistory();
    }
}

//Initialization of HMMemory

/**@Function:HMParameterSetting
 * @Params: TSPdata *tspdata
 * @comments:initializes the algorithm parameters (HMCR,PAR,BW)
 * @author:  Bill Gu
 */
void HMParameterSetting(double ConsideringRate, double AdjustingRate, double bandWidth)
{
    superHMCR = ConsideringRate;
    superPAR  = AdjustingRate;
    superBW   = bandWidth;
}

/**@Function:InitHMMemory
 * @Params: TSPdata *tspdata
 * @comments:Initialize the harmony memory HM with random values drawn from tspdata
 * @author:  Bill Gu
 */
void InitHMMemory(TSPdata *tspdata)
{
    //initialize the Global Harmony memory from TSPdata
    int j;
    int i;
//    int startpoint = RandomCityOfIndex(tspdata->n);
//    newHarmonyCity[0].X         = tspdata->x[startpoint];
//    newHarmonyCity[0].Y         = tspdata->y[startpoint];
//    newHarmonyCity[0].CityIndex = startpoint;
    
    for(i=0;i<HMS;i++)
    {
        int      dBest =0;
        TSPCITY  tspCityResult = {0,0,0,-1};
        //Initialze all Cities;
//        gHarmonyMemory[i][0].X         = tspdata->x[startpoint];
//        gHarmonyMemory[i][0].Y		   = tspdata->y[startpoint];
//        gHarmonyMemory[i][0].CityIndex = startpoint;
        for(j=0;j<CITYS;j++)
        {
            int iChooseCity = RandomCityOfIndex(tspdata->n); 	//get randon City's Index,Choose One City;
            while(CheckCityIndexExist(i,iChooseCity,j))           //Check this City does not exist;
            {
                iChooseCity = RandomCityOfIndex(tspdata->n);
            }
            gHarmonyMemory[i][j].X         = tspdata->x[iChooseCity];
            gHarmonyMemory[i][j].Y		   = tspdata->y[iChooseCity];
            gHarmonyMemory[i][j].CityIndex = iChooseCity;
            //this is not First Solution Route , Calculating the distance of the near city of first solution
            if(j!=0){
                dBest += dist(i,j); //consider Current is Best Solution
            }
        }
        dBest += ((int)(sqrt((gHarmonyMemory[i][CITYS-1].X-gHarmonyMemory[i][0].X)*(gHarmonyMemory[i][CITYS-1].X-gHarmonyMemory[i][0].X) + (gHarmonyMemory[i][CITYS-1].Y-gHarmonyMemory[i][0].Y)* (gHarmonyMemory[i][CITYS-1].Y-gHarmonyMemory[i][0].Y) ) + 0.5));
        // add the distance to startingpoint
        tspCityResult.CityIndex  =  i;
        tspCityResult.Fitness    =  dBest;
        UpdateHMMemory(i, CITYS, tspCityResult); //save the best Solution;
        //updateGeneralFitness(tspCityResult);     //update best and worst Distance
    }
    tspWorstSetting();
    tspBestSetting();
    //printf("the worst solution[%d]",tspWorst.CityIndex);
    //printf("the best  solution[%d]",tspBest.CityIndex);
    //PrintHarmoryMemory();
    //Harmony Memory has been Created;
}
/**@Function: InitialFirstLine
 * @Params: TSPdata *tspdata, get data from tspdata
 * @comment: generate a initial line of Harmony Memory and stored them in last line (HMS+1)
 * @author: Bill Gu
*/
void InitialFirstLine(TSPdata *tspdata)
{
    int j;
    for(j=0; j<tspdata->n; j++)
    {
        int iChooseCity = RandomCityOfIndex(tspdata->n);
        while(CheckCityIndexExist(HMS, iChooseCity, j))
        {
            iChooseCity = RandomCityOfIndex(tspdata->n);
        }
        gHarmonyMemory[HMS][j].X         = tspdata->x[iChooseCity];
        gHarmonyMemory[HMS][j].Y         = tspdata->y[iChooseCity];
        gHarmonyMemory[HMS][j].CityIndex = iChooseCity;
    }
    
}

/**@Function:GenerateHMByHMCRAndPAR
 * @Params:  int iter, the termination point indicator to define how many times will be done in HS algorithm, however, there is also a cpu_time checker violate which will also break out the loop
 * @Params:  TSPdata *tspdata
 * @Params:  int timelimit, limit the time consumed in HS, another termination point indicator
 * @comments: using HMCR And PAR to Generator a new Cityinfo Harmony
 * @author:  Bill Gu
 */
void GenerateHMByHMCRAndPAR(int iter, int timelimit, TSPdata *tspdata, Vdata *vdata)
{
    int iCurrentIter = 0;
    int i=0;
    double pPAR;
    TSPCITY tmpcity;
    
    //InitHMMemory(tspdata);
    while(iCurrentIter < iter)
    {
        //printf("iter[%d]begin\n",iCurrentIter);
//        dynPAR  = PAR  + (double)iCurrentIter*(double)(0.98/iter);
//        dynHMCR = HMCR + (double)iCurrentIter*(double)(0.98/iter);
            for(i=0;i<CITYS;i++)
            {
                int counter = 0;
                double pHMCR = RandomDoubleNext();
//                if(cpu_time() - vdata->starttime>0.2*timelimit)
//                {
//                    superHMCR = 0.5;
//                }
//                if(cpu_time() - vdata->starttime>0.5*timelimit)
//                {
//                    superHMCR = 0.99;
//                }
//                if(cpu_time() - vdata->starttime>0.75*timelimit)
//                {
//                    superHMCR = 0.999;
//                }
                if (pHMCR<superHMCR){
                    //If rand less than HMCR, then initializes the random memories of harmony above a set of harmonies
                    // When HMCR's Value is bigger ,it is  good for  local contract, otherwise it is good for groups
                    TSPINDEX tspIndexG = memoryConsideration(i);
                    for(counter=0;counter<10;counter++){
                        if(CheckCityIndexExist(-1,gHarmonyMemory[tspIndexG.row][tspIndexG.col].CityIndex, i)==1){
                            tspIndexG = memoryConsideration(i);
                        }
                        if(CheckCityIndexExist(-1,gHarmonyMemory[tspIndexG.row][tspIndexG.col].CityIndex, i)==0){
                            newHarmonyCity[i] = gHarmonyMemory[tspIndexG.row][tspIndexG.col];
                            break;
                        }
                    }
                    if(counter==10){
                        newHarmonyCity[i] = RandomCityGeneration(tspdata, i);
                    }
                    pPAR = RandomDoubleNext();
                    //If rand2 less than above-initialization PAR, just need to initialize trim bandwidth BW to adjust the harmonies, get a set of new harmony
//                    if(cpu_time() - vdata->starttime>0.2*timelimit)
//                    {
//                        superPAR = 0.2;
//                    }
//                    if(cpu_time() - vdata->starttime>0.5*timelimit)
//                    {
//                        superPAR = 0.5;
//                    }
//                    if(cpu_time() - vdata->starttime>0.75*timelimit)
//                    {
//                        superPAR = 0.9;
//                    }
                    if (pPAR<superPAR)
                    {
                        //printf("pPAR<PAR\n");
                        //do some adjustment here,
                        //x1(new,j)=x(new,j)*BW  adjustment with BW
                        tmpcity = newHarmonyCity[i];
                        newHarmonyCity[i] = pitchAdjustment(superBW,i,tspdata);
//                        if(CheckCityIndexExist(-1, newHarmonyCity[i].CityIndex, i)==1)
//                        {
//                            newHarmonyCity[i] = tmpcity;
//                        }
                        /*newly wrote*/
                        for(counter=0;counter<5;counter++){
                            if(CheckCityIndexExist(-1,newHarmonyCity[i].CityIndex, i)==1){
                                newHarmonyCity[i] = pitchAdjustment(superBW,i,tspdata);
                            }
                            if(CheckCityIndexExist(-1,newHarmonyCity[i].CityIndex, i)==0){
                                break;
                            }
                        }
                        if(counter==5){
                            newHarmonyCity[i] = tmpcity;
                        }
                        /*end here*/
                        //pitch tuning bandwidth BW: it says that be removed from memory tuning a set of harmonies with a certain probability, here is this adjustment.
                    }
                }else
                    {
                        newHarmonyCity[i] = RandomCityGeneration(tspdata, i);
                    }
            }
        calculateFitnessForRoute();
        iCurrentIter++;
        if(cpu_time() - vdata->starttime > timelimit)
        {
            tspBestSetting();
            break;
        }
    }
    tspBestSetting();
}


/**@Function:GeneratorNewHarmonyMemory
 * @Params: TSPdata *tspdata
 * @comments: using HMCR And PAR to Generator a new Cityinfo Harmony
 * @author:  Bill Gu
 */
void GeneratorNewHarmonyMemory(TSPdata *tspdata)
{
    //printf("GeneratorNewHarmonyMemory\n");
    
    int dBest =0;
    int j;
    TSPCITY  tspCityResult = {0,0,0,-1};
    for (j=0;j<CITYS;j++)
    {
        int iChooseCity = RandomCityOfIndex(tspdata->n); 	//get randon City's Index,Choose One City;
        while(CheckCityIndexExist(-1,iChooseCity,j))           //Check this City is not exist;
        {
            iChooseCity = RandomCityOfIndex(tspdata->n);
        }
        newHarmonyCity[j].X         = tspdata->x[iChooseCity];
        newHarmonyCity[j].Y		    = tspdata->y[iChooseCity];
        newHarmonyCity[j].CityIndex = iChooseCity;
        //this is not First Solution Route , Calculating the distance of the near city of first solution
        if(j!=0){
            dBest += distNew(j); //consider Current is Best Solution
        }
    }
    dBest += ((int)(sqrt((newHarmonyCity[CITYS-1].X-newHarmonyCity[0].X)*(newHarmonyCity[CITYS-1].X-newHarmonyCity[0].X) + (newHarmonyCity[CITYS-1].Y-newHarmonyCity[0].Y)* (newHarmonyCity[CITYS-1].Y-newHarmonyCity[0].Y) ) + 0.5));
    tspCityResult.Fitness =  dBest; //save the best Solution;
    UpdateHMMemory(-1, CITYS, tspCityResult);
}

/**@Function:UpdateHMMemory
 * @Params:  int i, Golbal Harmony Solution Index
 * @Params:  int j, Golbal Harmony City Index
 * @Params:  TSPCITY cityInfo, Golbal Harmony CityINfo
 * @comments: if the newly improvised harmony is better than the worst harmony in HM then replace it.
 * @author:  Bill Gu
 */
void UpdateHMMemory(int i,int j,TSPCITY cityInfo)
{
    if (i == -1)
    {
        newHarmonyCity[j]    = cityInfo;
    }
    else
    {
        gHarmonyMemory[i][j] = cityInfo;
    }
}

/**@Function:RandomCityGeneration
 * @Params:  TSPdata *tspdata
 * @Params:  int index, the column number of HM array
 * @comments:generate a city randomly
 * @author:  Bill Gu
 */

TSPCITY RandomCityGeneration(TSPdata *tspdata, int index)
{
    TSPCITY newCity;
    int iChooseCity = RandomCityOfIndex(tspdata->n); 	//get random City's Index,Choose One City;
    while(CheckCityIndexExist(-1,iChooseCity,index))           //Check this City is not exist;
    {
        iChooseCity = RandomCityOfIndex(tspdata->n);
    }
    newCity.X         = tspdata->x[iChooseCity];
    newCity.Y		  = tspdata->y[iChooseCity];
    newCity.CityIndex = iChooseCity;
    
    return newCity;
}

//After HS Convergence Local Search Adjustment

/**@Function: executionHybrid
 * @Params:  int timelimit, the total timelimit for entire algorithm
 * @Params:  Vdata *vdata, pass starttime to hybrid
 * @comments: to execute the hybrid algothrim after the area being reduced by HM
 * @author:  Bill Gu
 */
void executionHybrid(int timelimit, Vdata *vdata)
{
    int *path;
    int max_id;
    int counter;
    
    //Install the SIGINT/SIGTERM signal handlers:
    //install_sig_handlers();
    
    //Get command line options:
    //get_options(argc, argv);
    
    //Read input, get list of cities:
    if(verbose)
        printf("Reading input...\n");
    
    for(counter=0; counter<CITYS; counter++)
    {
        cities[counter]     = malloc(sizeof(struct city));
        cities[counter]->id = counter;
        cities[counter]->x  = (int) gHarmonyMemory[tspBest.CityIndex][counter].X;
        cities[counter]->y  = (int) gHarmonyMemory[tspBest.CityIndex][counter].Y;
    }
    num_cities = CITYS;
    //Initialize variables to hold paths:
    best_path = malloc((num_cities) * sizeof(int));
    path = malloc((num_cities) * sizeof(int));
    
    //Get simple list of city ids into the working path:
    max_id = get_list_of_cities(path);
    
    //Get matrix of distances between cities:
//    if(verbose)
//        printf("Calculating distances...\n");
    calc_distances(max_id);
    nearest_neighbor(path, num_cities);
    hybrid(timelimit, path, num_cities, vdata);
    //two_opt(path, num_cities);
    //anneal(path, num_cities);
    //Call nearest_neighbor algorithm to get a good first approximation:
//    if(verbose)
//        printf("Calling nearest neighbor algorithm...\n");
//    nearest_neighbor(path, num_cities);
//    
//    //Unless nearest neighbor is being used alone,
//    //call another algorithm to improve the answer:
//    if(!use_nearest_neighbor) {
//        
//        //Simulated Anneal:
//        if(use_anneal) {
//            if(verbose)
//                printf("Calling anneal...\n");
//            anneal(path, num_cities);
//        }
//        
//        //Two-opt:
//        else if(use_two_opt) {
//            if(verbose)
//                printf("Calling two-opt...\n");
//            two_opt(path, num_cities);
//        }
//        
//        //Default: Hybrid algorithm
//        else {
//            if(verbose)
//                printf("Calling hybrid algorithm...\n");
//            hybrid(path, num_cities);
//        }
//    }
    
    //Print solution:
    //print_solution();
    
    //Free memory allocated for distances matrix:
    free_distances();
    //    cities[7]->id = 33;
    //    printf("%d\n", num_cities);
    //    printf("%d\n", cities[7]->id);
    
    //return EXIT_SUCCESS;
}

/*Test Function*/
/*Ignore Please*/
void GreedyGeneration(TSPdata *tspdata, Vdata *vdata, int timelimit)
{
    int *path;
    int max_id;
    int counter;
    num_cities = tspdata->n;
    for(int i=0; i<99999999; i++){
        GeneratorNewHarmonyMemory(tspdata);
        for(counter=0; counter<tspdata->n; counter++)
        {
            cities[counter]     = malloc(sizeof(struct city));
            cities[counter]->id = counter;
            cities[counter]->x  = (int) newHarmonyCity[counter].X;
            cities[counter]->y  = (int) newHarmonyCity[counter].Y;
        }
    
        best_path = malloc((num_cities+1) * sizeof(int));
        path = malloc((num_cities+1) * sizeof(int));
    
        max_id = get_list_of_cities(path);
    
        calc_distances(max_id);
        Greedy(path, num_cities, tspdata);
        DataPassNew(vdata, tspdata);
        calculateFitnessForRoute();
        if(cpu_time() - vdata->starttime > timelimit)
        {
            tspBestSetting();
            break;
        }
    }
    free_distances();
}

/**@Function: executionGreedy
 * @Params: TSPdata *tspdata, get data from tspdata
 * @Params: Vdata *vdata, use to pass data
 * @comments: execute Greedy Algorithm to generate the first several lines in HS, others will be generated randomly
 * @author:  Bill Gu
 */
void executionGreedy(TSPdata *tspdata, Vdata *vdata, int fuko)
{
    int *path;
    int max_id;
    int counter;
    int i,j;
    
    num_cities = tspdata->n;
    
    InitialFirstLine(tspdata);
    for(counter=0; counter<tspdata->n; counter++)
    {
        cities[counter]     = malloc(sizeof(struct city));
        cities[counter]->id = counter;
        cities[counter]->x  = (int) gHarmonyMemory[HMS][counter].X;
        cities[counter]->y  = (int) gHarmonyMemory[HMS][counter].Y;
    }
    //Initialize variables to hold paths:
    best_path = malloc((num_cities+1) * sizeof(int));
    path = malloc((num_cities+1) * sizeof(int));
    
    //Get simple list of city ids into the working path:
    max_id = get_list_of_cities(path);
    
    calc_distances(max_id);
    Greedy(path, num_cities, tspdata);
        
    
    TSPCITY  tspCityResult = {0,0,0,-1};
    DataPass(0,vdata,tspdata);
    tspCityResult.CityIndex = 0;
    tspCityResult.Fitness   = DataPass(0, vdata, tspdata);
    UpdateHMMemory(0, CITYS, tspCityResult);
    //printf("%d\n", gHarmonyMemory[0][CITYS].Fitness);
    
    for(i=1;i<HMS;i++)
    {
        int      dBest =0;
        TSPCITY  tspCityResult = {0,0,0,-1};

        for(j=0;j<CITYS;j++)
        {
            int iChooseCity = RandomCityOfIndex(tspdata->n); 	//get randon City's Index,Choose One City;
            while(CheckCityIndexExist(i,iChooseCity,j))           //Check this City does not exist;
            {
                iChooseCity = RandomCityOfIndex(tspdata->n);
            }
            gHarmonyMemory[i][j].X         = tspdata->x[iChooseCity];
            gHarmonyMemory[i][j].Y		   = tspdata->y[iChooseCity];
            gHarmonyMemory[i][j].CityIndex = iChooseCity;
            //this is not First Solution Route , Calculating the distance of the near city of first solution
            if(j!=0){
                dBest += dist(i,j); //consider Current is Best Solution
            }
        }
        dBest += ((int)(sqrt((gHarmonyMemory[i][CITYS-1].X-gHarmonyMemory[i][0].X)*(gHarmonyMemory[i][CITYS-1].X-gHarmonyMemory[i][0].X) + (gHarmonyMemory[i][CITYS-1].Y-gHarmonyMemory[i][0].Y)* (gHarmonyMemory[i][CITYS-1].Y-gHarmonyMemory[i][0].Y) ) + 0.5));
        // add the distance to startingpoint
        tspCityResult.CityIndex  =  i;
        tspCityResult.Fitness    =  dBest;
        UpdateHMMemory(i, CITYS, tspCityResult); //save the best Solution;
    }
    tspWorstSetting();
    tspBestSetting();
    free_distances();
    
    for(;;)
    {
        InitialFirstLine(tspdata);
        for(counter=0; counter<tspdata->n; counter++)
        {
            cities[counter]     = malloc(sizeof(struct city));
            cities[counter]->id = counter;
            cities[counter]->x  = (int) gHarmonyMemory[HMS][counter].X;
            cities[counter]->y  = (int) gHarmonyMemory[HMS][counter].Y;
        }
        //Initialize variables to hold paths:
        best_path = malloc((num_cities+1) * sizeof(int));
        path = malloc((num_cities+1) * sizeof(int));
        
        //Get simple list of city ids into the working path:
        max_id = get_list_of_cities(path);
        
        calc_distances(max_id);
        Greedy(path, num_cities, tspdata);
        
        TSPCITY  tspCityResult = {0,0,0,-1};
        DataPass(2,vdata,tspdata);
        tspCityResult.CityIndex = 2;
        tspCityResult.Fitness   = DataPass(2, vdata, tspdata);
        UpdateHMMemory(2, CITYS, tspCityResult);
        
        if(gHarmonyMemory[2][CITYS].Fitness<gHarmonyMemory[0][CITYS].Fitness)
        {
            TSPCITY  tspCityResult = {0,0,0,-1};
            DataPass(0,vdata,tspdata);
            tspCityResult.CityIndex = 0;
            tspCityResult.Fitness   = DataPass(0, vdata, tspdata);
            UpdateHMMemory(0, CITYS, tspCityResult);
            //printf("%d\n", gHarmonyMemory[0][CITYS].Fitness);
        }
        
        free_distances();
        
        if(cpu_time() - vdata->starttime > fuko)
        {
            tspBest.CityIndex = 0;
            tspBest.Fitness   = gHarmonyMemory[0][CITYS].Fitness;
            break;
        }
    }
    
}

/**@Function: InitReversion
 * @Params: TSPdata *tspdata, get data from tspdata
 * @Params: Vdata *vdata, use to pass data
 * @comments: reverse the intialization
 * @author:  Bill Gu
 */
void InitReversion(TSPdata *tspdata, Vdata *vdata)
{
    int *path;
    int max_id;
    int counter;
    int i,j;
    
    num_cities = tspdata->n;
    
    for(i=0;i<9;i++)
    {
        InitialFirstLine(tspdata);
        for(counter=0; counter<tspdata->n; counter++)
        {
            cities[counter]     = malloc(sizeof(struct city));
            cities[counter]->id = counter;
            cities[counter]->x  = (int) gHarmonyMemory[HMS][counter].X;
            cities[counter]->y  = (int) gHarmonyMemory[HMS][counter].Y;
        }
        //Initialize variables to hold paths:
        best_path = malloc((num_cities+1) * sizeof(int));
        path = malloc((num_cities+1) * sizeof(int));
        
        //Get simple list of city ids into the working path:
        max_id = get_list_of_cities(path);
        
        calc_distances(max_id);
        Greedy(path, num_cities, tspdata);
        
        
        TSPCITY  tspCityResult = {0,0,0,-1};
        DataPass(i,vdata,tspdata);
        tspCityResult.CityIndex = i;
        tspCityResult.Fitness   = DataPass(i, vdata, tspdata);
        UpdateHMMemory(i, CITYS, tspCityResult);
    }
    
    for(i=9;i<HMS;i++)
    {
        int      dBest =0;
        TSPCITY  tspCityResult = {0,0,0,-1};
        
        for(j=0;j<CITYS;j++)
        {
            int iChooseCity = RandomCityOfIndex(tspdata->n); 	//get randon City's Index,Choose One City;
            while(CheckCityIndexExist(i,iChooseCity,j))           //Check this City does not exist;
            {
                iChooseCity = RandomCityOfIndex(tspdata->n);
            }
            gHarmonyMemory[i][j].X         = tspdata->x[iChooseCity];
            gHarmonyMemory[i][j].Y		   = tspdata->y[iChooseCity];
            gHarmonyMemory[i][j].CityIndex = iChooseCity;
            //this is not First Solution Route , Calculating the distance of the near city of first solution
            if(j!=0){
                dBest += dist(i,j); //consider Current is Best Solution
            }
        }
        dBest += ((int)(sqrt((gHarmonyMemory[i][CITYS-1].X-gHarmonyMemory[i][0].X)*(gHarmonyMemory[i][CITYS-1].X-gHarmonyMemory[i][0].X) + (gHarmonyMemory[i][CITYS-1].Y-gHarmonyMemory[i][0].Y)* (gHarmonyMemory[i][CITYS-1].Y-gHarmonyMemory[i][0].Y) ) + 0.5));
        // add the distance to startingpoint
        tspCityResult.CityIndex  =  i;
        tspCityResult.Fitness    =  dBest;
        UpdateHMMemory(i, CITYS, tspCityResult); //save the best Solution;
    }
    tspWorstSetting();
    tspBestSetting();
    free_distances();
}



//Nearest Neighbor Algorithm (Greedy Algorithm)


/**@Function: nearest_neighbor
 * @Params: int *path, contains the list of cities to build the path from.  At completion, contains the newly created path
 * @Params: int len, the length of the path or the number of cities
 * @comments: creates a path by progressively selecting the nearest neighbor to the last city added to the path (essentially Greedy Algorithm)
 * @author:  Bill Gu
 */
void nearest_neighbor(int *path, int len)
{
    int i, dst;
    
    dst = 0;
    for(i=0; i<len-1; i++)
    {
        dst += swap_closest(&path[i], (len-i));
    }
    
    dst += get_distance(path[len-1], path[0]);
    set_best(dst, path);
}

/**@Function: Greedy
 * @Params: int *path, contains the list of cities to build the path from. At completion, contains the newly created path
 * @Params: int len, the length of the path; here is the number of all cities instead of the number of cities needed to be passed
 * @Params: TSPdata *tspdata, get data from tspdata
 * @comments: the greedy search algorithm, always find next city nearest
 * @author: Bill Gu
*/
void Greedy(int *path, int len, TSPdata *tspdata)
{
    int i, dst;
    
    dst = 0;
    for(i=0; i<tspdata->min_node_num; i++)
    {
        dst += swap_closest(&path[i], (len-i));
    }
    
    dst += get_distance(path[tspdata->min_node_num], path[0]);
    best_distance = dst;
    copy_array(best_path, path, tspdata->min_node_num);
}


/**@Function: swap_closest
 * @Params: int *remaining, pointer to a (section of) a list of cities
 * @Params: int num_remaining, the number of cities remaining before the end of the list
 * @return: int, the distance from the city in the first position in the list to the city in the second position in the list
 * @comments: swaps into the second position in the list the city that is closest to the city in the first position in the list
 * @author:  Bill Gu
 */
int swap_closest(int *remaining, int num_remaining)
{
    int i, cur, best;
    
    cur  = remaining[0];
    best = 1;
    
    for(i=2; i<num_remaining; i++)
    {
        if(get_distance(cur, remaining[i]) < get_distance(cur, remaining[best]))
        {
            best = i;
        }
    }
    
    swap(1, best, remaining);
    return get_distance(cur, remaining[1]);
}


//2 OPT Algorithm:


/**@Function: two_opt
 * @Params: int *path, the path to improve upon
 * @Params: int len, the length of the path or the number of cities
 * @Params: int *path, the path to perfrom the swap on
 * @comments: iteratively improves on a given path by swapping 2 edges at a time in order to create a shorter path, the improved path is left at the location specified by the path pointer
 * @author:  Bill Gu
 */
void two_opt(int *path, int len)
{
    int i, j, dist;
    
    for(i=1; i<len; i++)
    {
        for(j=i; j<len; j++)
        {
            dist = two_opt_dist(best_distance, i, j, path, len);
            if(dist < best_distance)
            {
                if(debug)
                {
                    printf("Two-opt found new path with distance: %d\n", dist);
                }
                two_opt_swap(i, j, path);
                set_best(dist, path);
                i=1;
                break;
            }
        }
    }
}



/**@Function: two_opt_swap
 * @Params: int i, the first index of the section to reverse
 * @Params: int j, the last index of the section to reverse
 * @Params: int *path, the path to perfrom the swap on
 * @comments: performs a 2 opt swap by reversing the section of the path found between indicies i and j
 * @author:  Bill Gu
 */
void two_opt_swap(int i, int j, int *path)
{
    while (i < j)
    {
        swap(i++, j--, path);
    }
}



/**@Function: two_opt_dist
 * @Params: int old_dist, the distance of the path being changed
 * @Params: int i, the first index of the section to reverse
 * @Params: int j, the last index of the section to reverse
 * @Params: int *path, the path that the swap would be performed on
 * @Params: int len, the length of the path or the number of cities
 * @return: int, the length of the path that would result from the associated 2 opt swap
 * @comments: returns the distance of the path that would result from the associated 2-opt swap without actually performing the swap
 * @author:  Bill Gu
 */
int two_opt_dist(int old_dist, int i, int j, int *path, int len)
{
    int new_dist;
    
    if(j == len-1)
    {
        new_dist = old_dist - (get_distance(path[i-1], path[i]) + get_distance(path[j], path[0]));
        new_dist += get_distance(path[i-1], path[j]) + get_distance(path[i], path[0]);
    }
    else
    {
        new_dist = old_dist - (get_distance(path[i-1], path[i]) + get_distance(path[j], path[j+1]));
        new_dist += get_distance(path[i-1], path[j]) + get_distance(path[i], path[j+1]);
    }
    
    return new_dist;
}


//SA(Simulated Annealing):


/**@Function: anneal
 * @Params: int *path, the path to perform simulated annealing on
 * @Params: int len, the length of the path or the number of cities
 * @comments: SA simulated annealing
 * @author:  Bill Gu
 */
void anneal(int *path, int len)
{
    int i, j, dst, swp_dst, attempt;
    double temp;
    
    //Seed random number generator:
    //srand((unsigned)time(NULL));
    
    //Get the current path's distance:
    dst = calc_path_dist(path, len);
    
    temp = START_TEMP;
    attempt = 0;
    while(attempt < SATISFIED)
    {
        
        //Try a random 2-opt swap:
        i = (rand() % (len-1)) + 1;
        j = i + (rand() % (len-i));
        swp_dst = two_opt_dist(dst, i, j, path, len);
        
        //If the result is acceptable:
        if(anneal_accept(swp_dst, dst, temp))
        {
            if(debug)
            {
                printf("Anneal: temp: %f, old path: %d, new path : %d", temp, dst, swp_dst);
                if(swp_dst > dst)
                    printf("\t < escape local optimum");
                printf("\n");
            }
            
            //Update the path and it's distance:
            two_opt_swap(i, j, path);
            dst = swp_dst;
            
            //Update the global best dst/path, if necessary:
            if(dst < best_distance)
                set_best(dst, path);
            
            //Reset attempt counter:
            attempt = 0;
        }
        
        //If the result is unacceptable:
        else
            attempt++;
        if(debug)
            printf("Decline #%d\n", attempt);
        
        //Decrease the temperature:
        temp = change_temp(temp);
    }
}



/**@Function: anneal_accept
 * @Params: int new_dst, the distance of the new path
 * @Params: int old_dst, the distance of the old path
 * @return: int, 1 if the move is accepted, 0 if not
 * @comments: accepts or rejects a proposed change from a path with old_dst to a path with new_dst, given temperature
 * @author:  Bill Gu
 */
int anneal_accept(int new_dst, int old_dst, double temp)
{
    double prob, q;
    
    if(new_dst == old_dst)
        return 0;
    
    prob = exp((old_dst - new_dst)/temp);
    q = rand() / (double) RAND_MAX;
    
    if(q < prob)
        return 1;
    else
        return 0;
}



/**@Function: change_temp
 * @Params: double old_temp, the temperature to decrease
 * @return: double, the new, decreased remperature
 * @comments: takes a temperature and returns a decreased temperature
 * @author:  Bill Gu
 */
double change_temp(double old_temp)
{
    if(old_temp > 0.01)
    {
        return old_temp*DELTA_TEMP;
    }
    return old_temp;
}


//Super Hybrid Algorithm:


/**@Function: hybrid
 * @Params: int timelimit, termination criterion indictor, violate then break out
 * @Params: int *path, the path to improve
 * @Params: int len, the length of the path or the number of cities
 * @Params: Vdata *vdata
 * @comments: a combination of 2 opt and simulated annealing algorithms
 * @author:  Bill Gu
 */
void hybrid(int timelimit, int *path, int len, Vdata *vdata)
{
    int i, j, dst, swp_dst, change, best_change, term_cnt;
    double temp;
    
    //Get the current path's distance:
    dst = calc_path_dist(path, len);
    
    term_cnt = 0;
    temp = 0.01;
    
    do {
        change = 0;
        best_change = 0;
        
        //Iterate through endpoints to be swapped:
        for(i=1; i<len; i++)
        {
            for(j=i; j<len; j++)
            {
                
                //Get the path distance of the tentative two-opt swap
                //(without actually performing the swap):
                swp_dst = two_opt_dist(dst, i, j, path, len);
                
                //Check to see whether the new distance is acceptable:
                if(anneal_accept(swp_dst, dst, temp))
                {
                    
                    //Print debug information:
                    if(debug)
                    {
                        printf("Hybrid: temp: %f, old path: %d, new path : %d", temp, dst, swp_dst);
                        if(swp_dst > dst)
                            printf("\t < up");
                        printf("\n");
                    }
                    
                    //Make the swap:
                    two_opt_swap(i, j, path);
                    dst = swp_dst;
                    change = 1;
                    
                    //If necessary, update the running best path/distance:
                    if(dst < best_distance)
                    {
                        set_best(dst, path);
                        best_change = 1;
                        term_cnt=0;
                    }
                }
            }
        }
        
        //Increment the termination counter if there's been no change
        //in the best path throughout the entire last loop:
        if(!best_change)
        {
            term_cnt++;
        }
        
        //Boost the temperature if there was no change
        //at all throughout the entire last loop
        //(i.e. it's stuck in a local optimum):
        if(!change)
        {
            
            //The closer the algorithm is to termination,
            //the larger the temperature boost (getting desperate):
            
            temp = get_max(START_TEMP,  END_TEMP * ((double)term_cnt/SATISFIED));
        }
        
        //If there has been some change, put it on ice
        //(this will get down it back down to a local optimum):
        else {
            temp = 0.01;
        }
        if(cpu_time() - vdata->starttime > timelimit)
        {
            break;
        }
        
    } while(term_cnt<SATISFIED);

}




//Utilities


/**@Function: set_best
 * @Params: int distance, the distance of the new path
 * @Params: int *path, the new path
 * @comments: sets the static variables which hold the best path found thus far
 * @author:  Bill Gu
 */
void set_best(int distance, int *path)
{
    sigset_t set;
    
    sigemptyset(&set);
    sigaddset(&set, SIGTERM);
    sigaddset(&set, SIGINT);
    sigprocmask(SIG_BLOCK, &set, NULL);
    
    if(verbose)
        printf("New best path found: %d\n", distance);
    best_distance = distance;
    copy_array(best_path, path, num_cities);
    
    sigprocmask(SIG_UNBLOCK, &set, NULL);
}



/**@Function: swap
 * @Params: int i, the First index
 * @Params: int j, the Second index
 * @comments: swaps elements at indicies i and j in array
 * @author:  Bill Gu
 */
void swap(int i, int j, int *array)
{
    int temp;
    
    temp = array[i];
    array[i] = array[j];
    array[j] = temp;
}



/**@Function: copy_array
 * @Params: int *to, the array to copy elements to
 * @Params: int *from, the array to copy elements from
 * @comments: copies the elements of one array into another array
 * @author:  Bill Gu
 */
void copy_array(int *to, int *from, int len)
{
    int i=0;
    
    for(i=0; i<len; i++) {
        to[i] = from[i];
    }
}



/**@Function: get_max
 * @Params: double a, the First element
 * @Params: double b, the Second element
 * @return: double, the max of the two elements
 * @comments: returns the max of the two values passed in
 * @author:  Bill Gu
 */
double get_max(double a, double b)
{
    if(a>b)
        return a;
    else
        return b;
}


//INPUT/OUTPUT of the super hybrid algorithm, for debugging use:


/**@Function: get_options
 * @Params: int argc, the length of the argument vector
 * @Params: char **argv, the argument vector
 * @comments: I build this function to get the command line options and sets the appropriate static variables representing them. This function is used in debugging because I want to test all the algorithms including the 2 opt and adjusted SA in order to make best choice and build the hybrid algorithm properly.
 * @author:  Bill Gu
 */
void get_options(int argc, char **argv)
{
    char opt;
    
    while((opt = getopt(argc, argv, "adf:hntv")) != -1)
    {
        switch(opt)
        {
            case 'a':
                use_anneal = 1;
                break;
            case 'n':
                use_nearest_neighbor = 1;
                break;
            case 'f':
                in_file = 1;
                out_file = 1;
                strncpy(in_filename, optarg, WORD_MAX);
                snprintf(out_filename, WORD_MAX, "%s%s", in_filename, ".tour");
                break;
            case 't':
                use_two_opt = 1;
                break;
            case 'v':
                verbose = 1;
                break;
            case 'd':
                verbose = 1;
                debug = 1;
                break;
            case 'h':
            default:
                printf("Usage: %s -[adntv] -[f filename]\n", argv[0]);
                printf("Algorithms:\n");
                printf("\t-Default: super hybrid algorithm\n");
                printf("\t-n: Nearest Neighbor (only)\n");
                printf("\t-t: Two-opt\n");
                printf("\t-a: Simulated Anneal\n");
                printf("Display modes:\n");
                printf("\t-v: Verbose (minor progress messages)\n");
                printf("\t-d: Debug (lots of detailed messages)\n");
                printf("Input/Output:\n");
                printf("\t-f: Specify file to use as input/source file\n");
                printf("\t    Note: this will result in a output file named [input file].tour\n");
                
                exit(EXIT_SUCCESS);
        }
    }
}



/**@Function:get_list_of_cities
 * @Params:  int *list, city list
 * @return:  max_id, max city id
 * @comments:gets a list of all cities into the array specified by list it is useful for setting up the initial path
 * @author:  Bill Gu
 */
int get_list_of_cities(int *list)
{
    int i, max_id;
    max_id = cities[0]->id;
    
    for(i=0; i<num_cities; i++)
    {
        list[i] = cities[i]->id;
        if(list[i] > max_id)
        {
            max_id = list[i];
        }
    }
    return max_id;
}


/**@Function: print_distances
 * @comments: prints all the distances between all the cities
 * @author:  Bill Gu
 */
void print_distances()
{
    int i, j;
    for(i=0; i<num_cities; i++)
    {
        for(j=i+1; j<num_cities; j++)
        {
            print_distance(cities[i]->id, cities[j]->id);
        }
    }
}



/**@Function: print_distance
 * @Params: int i, the id of the first city
 * @Paramas: int j, the id of the second city
 * @comments: prints the distance between the two cities which are represented by i and j (IDs)
 * @author:  Bill Gu
 */
void print_distance(int i, int j)
{
    printf("Distance between %d and %d: %d\n", i, j, get_distance(i, j));
}



/**@Function: print_cities
 * @comments: prints a list of all cities and their cordinates
 * @author:  Bill Gu
 */
void print_cities()
{
    int i;
    for(i=0; i<num_cities; i++)
    {
        print_city(cities[i]);
    }
}


/**@Function: print_city
 * @Params: city *c, the city to be printed
 * @comments: prints a single city's id and coordinates
 * @author:  Bill Gu
 */
void print_city(city * c)
{
    printf("City: %d, X: %d, Y: %d\n", c->id, c->x, c->y);
}



/**@Function: print_solution
 * @comments: prints the solution to stdout or a file for intermediate checking of hybrid algorithm, it is no need showed in the entire program
 * @author:  Bill Gu
 */
void print_solution()
{
    int i;
    //    FILE * f;
    //
    //    if(out_file)
    //        f = fopen(out_filename, "w");
    //    else
    //        f = stdout;
    
    printf("%d\n", best_distance);
    for(i=0; i<num_cities; i++)
    {
        printf("%d\n", gHarmonyMemory[tspBest.CityIndex][best_path[i]].CityIndex);
    }
    
    //    if(out_file)
    //        fclose(f);
}


/**@Function: calc_distances
 * @Params: int max_id, the max id of all cities
 * @comments: calculates the distances between all cities, storing them in the static distances matrix
 * @author:  Bill Gu
 */
void calc_distances(int max_id)
{
    int i, j;
    unsigned long sum;
    
    sum=0;
    distances = malloc((max_id+1) * sizeof(int *));
    
    for(i=0; i<num_cities; i++)
    {
        distances[cities[i]->id] = malloc((max_id+1) * sizeof(int));
        
        for(j=i+1; j<num_cities; j++)
        {
            sum += distances[cities[i]->id][cities[j]->id] = calc_distance(cities[i], cities[j]);
        }
    }
    avg_distance = (sum/pow(max_id, 2));
}


/**@Function: calc_distance
 * @Params: city *a, The first city
 * @Params: city *b, The second city
 * @return: the distance rounded to the nearest integer, between the two cities
 * @comments: calculates the distances between all cities, storing them in the static distances matrix
 * @author:  Bill Gu
 */
int calc_distance(city *a, city *b)
{
    int x_dif, y_dif;
    double distance;
    
    x_dif = a->x - b->x;
    y_dif = a->y - b->y;
    
    distance = sqrt(pow((double) x_dif, 2.0) + pow((double) y_dif, 2.0));
    
    return (int) round(distance);
}


/**@Function: get_distance
 * @comments: returns the distance between cities with ids i and j
 * @Params: int i, The first city's id
 * @Params: int j, The second city's id
 * @return: the distance between the two cities
 * @author:  Bill Gu
 */
int get_distance(int i, int j)
{
    if(i<j)
        return distances[i][j];
    else
        return distances[j][i];
}

/**@Function: calc_path_dist
 * @Params: int *path, the path needed to be calculated
 * @Params: int len, the length of the path or the number of cities
 * @comments: distance calculation
 * @author:  Bill Gu
 */
int calc_path_dist(int * path, int len)
{
    int i, dist;
    
    dist = 0;
    for(i=0; i<len-1; i++)
    {
        dist += get_distance(path[i], path[i+1]);
    }
    dist += get_distance(path[i], path[0]);
    return dist;
}

/**@Function: free_distance
 * @comments: to free the matrix of distances allocated by calc_distance
 * @author:  Bill Gu
 */
void free_distances()
{
    int i;
    
    for(i=0; i<num_cities; i++)
    {
        free(distances[(cities[i]->id)]);
    }
    free(distances);
}


//SIGNAL HANDLERS:
/**@Function: sig_handler
 * @Params:  int sig, the signal received
 * @comments: Signal handler for SIGINT or SIGTERM signals.
 * @comments: Prints solution before exiting
 * @author:  Bill Gu
 */
void sig_handler(int sig)
{
    if(verbose)
        printf("Received signal %d: exiting...\n", sig);
    print_solution();
    exit(EXIT_SUCCESS);
}

/**@Function: install_sig_handlers
 * @comments: installs the signal handler for SIGINT and SIGTERM
 * @author:  Bill Gu
 */
void install_sig_handlers()
{
    struct sigaction siga;
    
    siga.sa_handler = sig_handler;
    sigemptyset(&siga.sa_mask);
    siga.sa_flags = 0;
    
    sigaction(SIGTERM, &siga, NULL);
    sigaction(SIGINT, &siga, NULL);
}

/**@Function:ReturnBestSolution
 * @Params:  Vdata *vdata, send data to vdata->bestsol in the main
 * @Params:  TSPdata *tspdata
 * @comments: return best solution find so far pass it to bestsol
 * @author:  Bill Gu
 */

void ReturnBestSolution(Vdata *vdata,TSPdata *tspdata)
{
        for(int i=0;i<CITYS;i++)
        {
            vdata->bestsol[i] = gHarmonyMemory[tspBest.CityIndex][best_path[i]].CityIndex;
        }
        for(int tmp=CITYS;tmp<tspdata->n;tmp++)
        {
            vdata->bestsol[tmp] = -1;
        }
}



/**@Function:DataPass
 * @Params: int linenumber
 * @Params: Vdata *vdata
 * @Params: TSPdata *tspdata
 * @return: dBest, the distance of each line after calculated
 * @comments: pass the data back to the Harmony Search Memory (after adjusted by hybrid algorithm) if necessary
 * @author:  Bill Gu
 */
int DataPass(int linenumber, Vdata *vdata, TSPdata *tspdata)
{
    int j;
    int dBest = 0;
    for(j=0;j<CITYS;j++)
    {
        gHarmonyMemory[linenumber][j].X = tspdata->x[gHarmonyMemory[HMS][best_path[j]].CityIndex];
        gHarmonyMemory[linenumber][j].Y = tspdata->y[gHarmonyMemory[HMS][best_path[j]].CityIndex];
        gHarmonyMemory[linenumber][j].CityIndex = gHarmonyMemory[HMS][best_path[j]].CityIndex;
        if(j!=0){
            dBest += dist(linenumber,j);
        }
    }
    dBest += ((int)(sqrt((gHarmonyMemory[linenumber][CITYS-1].X-gHarmonyMemory[linenumber][0].X)*(gHarmonyMemory[linenumber][CITYS-1].X-gHarmonyMemory[linenumber][0].X) + (gHarmonyMemory[linenumber][CITYS-1].Y-gHarmonyMemory[linenumber][0].Y)* (gHarmonyMemory[linenumber][CITYS-1].Y-gHarmonyMemory[linenumber][0].Y) ) + 0.5));
    gHarmonyMemory[linenumber][CITYS].Fitness = dBest;
    
    return dBest;
}

/*Test Function*/
/*it may be no use*/
void DataPassNew(Vdata *vdata, TSPdata *tspdata)
{
    for(int i=0;i<CITYS;i++)
    {
        newHarmonyCity[i].X = tspdata->x[newHarmonyCity[best_path[i]].CityIndex];
        newHarmonyCity[i].Y = tspdata->y[newHarmonyCity[best_path[i]].CityIndex];
        newHarmonyCity[i].CityIndex = newHarmonyCity[best_path[i]].CityIndex;
    }
}



