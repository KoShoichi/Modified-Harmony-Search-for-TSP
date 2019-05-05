//
//  hstsp.h
//  HStest
//
//  Created by 顾肖一 on 2017/05/17.
//  Copyright © 2017年 Bill Gu. All rights reserved.
//

#ifndef hstsp_h
#define hstsp_h

#define MAX_STR 1024
typedef struct{
    char     name[MAX_STR];         /* name of the instance */
    int      n;                     /* number of nodes */
    double   *x;                    /* x-coordinates of nodes */
    double   *y;                    /* y-coordinates of nodes */
    int      min_node_num;          /* minimum number of nodes the solution contains */
} TSPdata;      /* data of TSP instance */

typedef struct{
    double X;                          //City x coordinate
    double Y;						   //City y coordinate
    int    Fitness;                    //Fitness Solution
    int    CityIndex;                  //While this city has been Chosen;
} TSPCITY;

typedef struct {
    double        timebrid;       /* the time before reading the instance data */
    double        starttime;      /* the time the search started */
    double        endtime;        /* the time the search ended */
    int           *bestsol;       /* the best solution found so far */
    /* NEVER MODIFY THE ABOVE FOUR VARIABLES. */
    /* You can add more components below. */
    
} Vdata;                /* various data often necessary during the search */

typedef struct {
    int    timelim;              /* the time limit for the algorithm in secs. */
    int    givesol;              /* give a solution (1) or not (0) */
    int    outformat;            /* 1: output the computed tour in TSPLIB format;
                                  0: do not output it */
    char   tourfile[MAX_STR];    /* the output file of computed tour */
    /* NEVER MODIFY THE ABOVE VARIABLES.  */
    /* You can add more components below. */
    
} Param;                /* parameters */

//Information for a single city, as read in:
typedef struct city{
    int id;
    int x;
    int y;
} city;

void PrintHarmoryMemory();
void PrintOnlyBestAndWorst();
void InitHMMemory(TSPdata *tspdata);
void GenerateHMByHMCRAndPAR(int iter,int timelimit,TSPdata *tspdata, Vdata *vdata);
void ReturnBestHSSolution(Vdata *vdata,TSPdata *tspdata);
void HMParameterSetting(double ConsideringRate, double AdjustingRate, double bandWidth);
void getMinNode(TSPdata *tspdata);
TSPCITY RandomCityGeneration(TSPdata *tspdata, int i);
void executionHybrid(int timelimit, Vdata *vdata);
void ReturnBestSolution(Vdata *vdata,TSPdata *tspdata);
int DataPass(int linenumber,Vdata *vdata,TSPdata *tspdata);
void executionGreedy(TSPdata *tspdata, Vdata *vdata, int fuko);
void GreedyGeneration(TSPdata *tspdata, Vdata *vdata, int timelimit);
void InitReversion(TSPdata *tspdata, Vdata *vdata);

#endif /* hstsp_h */
