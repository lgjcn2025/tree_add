#include <time.h>
#include "ladd.c"

int main(argc, argv)

    int argc;
char *argv[];

{
    PreprocData *pd = walloc(1,sizeof(PreprocData));
    pd->P = 25;
    pd->total = 10;
    pd->K = atoi(argv[1]);
    pd->meUpbound = 1;
    int runType = 0;



    loadGraphData(pd);

	struct timespec time_start={0,0},time_end={0,0};

	clock_gettime(CLOCK_REALTIME,&time_start);
    initPreprocData(pd);
	clock_gettime(CLOCK_REALTIME,&time_end);
    double tm1 = (10e9*time_end.tv_sec +time_end.tv_nsec - 10e9*time_start.tv_sec - time_start.tv_nsec)/10e9;
	
	clock_gettime(CLOCK_REALTIME,&time_start);
    printf("c K=%ld, total pass=%d, percent=%d%%\n", pd->K,pd->total,pd->P);

    if(runType == 0){
        
        cType cntSmallDeg = heuristicSmallwDegConn(pd);
        
        preProc(pd);
		clock_gettime(CLOCK_REALTIME,&time_end);
        double tm2 = (10e9*time_end.tv_sec +time_end.tv_nsec - 10e9*time_start.tv_sec - time_start.tv_nsec)/10e9;
		
		clock_gettime(CLOCK_REALTIME,&time_start);

        buildMCT(pd);

		clock_gettime(CLOCK_REALTIME,&time_end);

        double tm3 = (10e9*time_end.tv_sec +time_end.tv_nsec - 10e9*time_start.tv_sec - time_start.tv_nsec)/10e9;

    
        printf("the final total time: %lf (%lf %lf %lf)\n",(tm1+tm2+tm3),tm1,tm2,tm3);

        
        cType cnt2 = 0;
        cnt2 = cycleConnShieldEdges(pd);

        printf("c edge adding: %ld conn deg<K pairs, %ld more cyclely added edges\n",cntSmallDeg,cnt2);
        outputGraph(pd);

    }

  

    exit(0);

    
}
