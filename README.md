//Copyright @ Xiaoyi (Bill) GU
//Master of Information Systems
//School of Engineering
//University of Melbourne, Level 5, 173, Wilson Ave, Parkville VIC 3052
//For any problems or academic use please contact                                  
“xiaoyig@student.unimelb.edu.au”
//Please do not publicize it without permission

*****************************************************************************
*                                                                           *
*                 Modified Harmony Search Algorithm for TSP                 *
*                                        *
*****************************************************************************

Harmony Search (HS) is a novel intelligent algorithm firstly designed by Geem, Kim and Loganathan (2001). This is a solver based on HS and it is modified by introducing greedy algorithm for early initialization, which can be proved to obtain a much faster convergence speed. In order to ensure the intensification and diversification of Harmony Memory, and to refrain from being tendentious towards similar solutions during early period, it makes the starting-point of greedy paths distinctly. Decent combination of local search has also be utilized to avoid the uncertainty caused by randomness as far as possible. 


To compile:

type “make” to compile

To execute:

using stdin to read file
type “cat filename | ./main” to read file
give parameter “timelim” to control execution time (default is 300s without typing in)

execute with user-defined time limitation
type “cat filename | ./main timelim 200” (redefine time limitation as 200s)


*****************************************************************************
Reference:
Geem, ZW, Kim, JH & Loganathan, GV 2001, ‘A new heuristic optimization algorithm: harmony search’, Simulation, vol. 76, no.2, pp. 60-68.


