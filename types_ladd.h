




typedef unsigned long long int excessType; 


typedef unsigned long cType;
typedef unsigned int sType;
typedef unsigned char aType;

typedef 
   struct edgeProp
{

   cType endNode;
   cType cap;
   cType w; 
   cType avgCV; 
   long tmp;
   long cost;

   struct edgeProp* rev; 

}edgeP;


typedef  
   struct nodeProp
{
   edgeP* edges;
   cType maxEdges;
   cType nIdx;
   cType totalCap;
   
   sType* orderedEdges;

} nodeP;

typedef
struct NodePropExtra_
{
   cType fa;
   cType dep;
   long cv;
   short s;


} NodePropExtra;



typedef
struct NodePropArr_{

   cType * pfa;
   cType * pdep;
   long * pcv;
   long * poh; 
               
   long * pcof; 
   short * ps;
   cType * pacc_upid; 
   long* pacc_upmincv; 
   cType * pacc_pos_upmincv; 
   long ** pacc_cut_adjust; 
   aType ** pacc_cut_adjust_sign; 
   cType * pacc_jointid; 
   
   
} NodePropArr;

typedef
struct GraphData_{

  long N, M;
  nodeP *nodes;

} GraphData;








typedef
struct RandomData_{
  cType *randNums;
  cType randNumIdx;
  cType maxLen;
} RandomData;



typedef
struct PreprocData_{

  
  GraphData *gd;

  
  NodePropArr* allResults;

  
  RandomData* rd;

  
  cType *gpfa;
  cType *gpdep;
  long *gpcv;
  long *gpoh;
  long *gpnc;
  short *gps;

  long *gpaccmcv;
  cType *gprn;
  cType *gprc;
  long **gcutadj;
  aType **gcutadjsign;

  

  cType *roots; 

  cType *nodeTreeId; 

  int mode; 

  int P; 
  int total; 
  int SPAN_LEN; 
  int LEVEL_NUM;

  cType K;

  cType meUpbound; 

} PreprocData;

