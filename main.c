#include <time.h>
#include "ladd.c"

int main(argc, argv)

    int argc;
char *argv[];

{
    PreprocData *pd = walloc(1,sizeof(PreprocData));
    //no mode 2, so it is just randomly picking up children
    pd->P = 25; // P% percent of total passes in mode 1, the remaining in mode 2, 300 is mode 3
    pd->total = 10; //number of total passes
    pd->K = atoi(argv[1]);
    pd->meUpbound = 1; // one side with <= this number is a biased cut, the node in this side will be surely used to connect
    int runType = 0; //0 out algorithm , 1/2/3/4 are heurisitcs in sadd-master-multi, validating in treecut-master-checkpair



    loadGraphData(pd); //load graph data from standard input

	struct timespec time_start={0,0},time_end={0,0};

	clock_gettime(CLOCK_REALTIME,&time_start);
    initPreprocData(pd); //init data structure
	clock_gettime(CLOCK_REALTIME,&time_end);
    double tm1 = (10e9*time_end.tv_sec +time_end.tv_nsec - 10e9*time_start.tv_sec - time_start.tv_nsec)/10e9;
	
	clock_gettime(CLOCK_REALTIME,&time_start);
    printf("c K=%ld, total pass=%d, percent=%d%%\n", pd->K,pd->total,pd->P);

    if(runType == 0){
        
        cType cntSmallDeg = heuristicSmallwDegConn(pd);
        
        preProc(pd); // preproc by traversing the graph for pd->total times
		clock_gettime(CLOCK_REALTIME,&time_end);
        double tm2 = (10e9*time_end.tv_sec +time_end.tv_nsec - 10e9*time_start.tv_sec - time_start.tv_nsec)/10e9;
		
		clock_gettime(CLOCK_REALTIME,&time_start);

        buildMCT(pd);

		clock_gettime(CLOCK_REALTIME,&time_end);

        double tm3 = (10e9*time_end.tv_sec +time_end.tv_nsec - 10e9*time_start.tv_sec - time_start.tv_nsec)/10e9;

        //calcuRandomPairs(1,pd); // randomly choose 100 node pairs and calcu their min-cut and output
    
        printf("the final total time: %lf (%lf %lf %lf)\n",(tm1+tm2+tm3),tm1,tm2,tm3);

        //conn deg 1 nodes, a quick step that effectively for K=2 from 1


        //then conn all candi edge head by tail
        cType cnt2 = 0;
        cnt2 = cycleConnShieldEdges(pd);

        printf("c edge adding: %ld conn deg<K pairs, %ld more cyclely added edges\n",cntSmallDeg,cnt2);
        outputGraph(pd);

    }

  

    exit(0);

    
}
