
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <values.h>
#include <math.h>

#include "types_ladd.h"   
#include "parser_ladd.c"  
#include "timer.c"         

#define MAX_LONG LONG_MAX

#define min(s, t) ((s) < (t) ? (s) : (t))
#define max(s, t) ((s) > (t) ? (s) : (t))









void *walloc(unsigned int num, unsigned int size)
{
  void *ptr = calloc(num, size);
  assert(ptr != NULL);
  return ptr;
}



StochasData* initrand2(StochasData *rd )
{

  rd->randNumIdx = 0;

  for (int i = 0; i < rd->maxLen; i++)
  {
    rd->randNums[i] = (cType)rand();
  }
  
  return rd;

}

StochasData* initrand(cType len)
{
  StochasData *rd = walloc(1, sizeof(StochasData));
  rd->maxLen = len;
  rd->randNums = (cType *)walloc(rd->maxLen, sizeof(cType));
  
  initrand2(rd);

  return rd;
}


cType mrand(StochasData* rd)
{
  if(rd->randNumIdx >= rd->maxLen){
		initrand2(rd);
  }	  
  return rd->randNums[rd->randNumIdx++];
}


void HeapAdjustDown(sType *idx, edgeP * edges ,int start,int end)  
{  
    sType tempIdx = idx[start];  

    int i = 2*start+1;      
    
    
    
    while(i<=end)  
    {  
        if(i+1<=end && edges[idx[i+1]].tmp > edges[idx[i]].tmp )    
            i++;  

        if(edges[idx[i]].tmp <= edges[tempIdx].tmp )   
            break;  

        idx[start] = idx[i];

        start = i;  
        i = 2*start+1;  
    }  

    idx[start] = tempIdx;  
}  
  
void HeapSort(sType *idx, edgeP * edges, int len)  
{  

    int i;  
    for(i=(len-1)/2;i>=0;i--){  
        HeapAdjustDown(idx,edges,i,len-1);  
    }

    for(i=len-1;i>0;i--)
    {  
        
        sType temp = idx[i];  
        idx[i] = idx[0];  
        idx[0] = temp;  

        HeapAdjustDown(idx,edges,0,i-1);  
    }  

}  

#define MARK_CUT 9
#define MARK_CUT_SHIELD 31
long total_cut_edge = 0;
long total_cut_edge_lessK = 0;
  

void deOrderLinkByStochas(nodeP *np,PreprocData *pd)
{

  cType cnt = np->nIdx;

  sType *idxs = np->orderedLinks;
  edgeP *pedges = np->edges;

  for (int i = 0; i < cnt; i++)
  {
    
    pedges[i].tmp = rand()%1000;
    if(pedges[i].w != MARK_CUT){
      pedges[i].tmp += 10000; 
    } 
  }

    assert(cnt<4 || idxs[2]!=idxs[3]);
    HeapSort(idxs,pedges,cnt);
    assert(cnt<4 || idxs[2]!=idxs[3]);  
}


void deOrderLinkByStochasCap(nodeP *np,PreprocData *pd)
{

  cType cnt = np->nIdx;

  sType *idxs = np->orderedLinks;
  edgeP *pedges = np->edges;

  for (int i = 0; i < cnt; i++)
  {
    
    pedges[i].tmp = 1000-pedges[i].cap+1; 
  }

    assert(cnt<4 || idxs[2]!=idxs[3]);
    HeapSort(idxs,pedges,cnt);
    assert(cnt<4 || idxs[2]!=idxs[3]);  
}


void aOrderLinkByAvgCV(nodeP *np,PreprocData *pd)
{
  if (np->nIdx == 0)
  {
    return;
  }
  int cnt = np->nIdx;

  sType *idxs = np->orderedLinks;
  edgeP *pedges = np->edges;

  for (int i = 0; i < cnt; i++)
  {
    edgeP *pp = pedges +i;
    long acv = pp->avgCV;
    if(acv == 0 ){
      pp->tmp = MAX_LONG;
    }
    else{
      pedges[i].tmp = mrand(pd->rd) % acv; 
    }
  }

    HeapSort(idxs,pedges,cnt);
}


void aOrderLinkByDegree(nodeP *np,PreprocData *pd)
{
  if (np->nIdx == 0)
  {
    return;
  }
  int cnt = np->nIdx;
  nodeP* nodes = pd->gd->nodes;
  sType *idxs = np->orderedLinks;
  edgeP *pedges = np->edges;

  for (int i = 0; i < cnt; i++)
  {
    edgeP *pp = pedges +i;
    cType zn = pp->endNode;
    pedges[i].tmp =  -1 * (nodes+zn)->nIdx;
  }

    HeapSort(idxs,pedges,cnt);
}


























long gRoot = 0;





















  


































































































void calcuTotalCap(PreprocData *pd)
{
  nodeP *nodes = pd->gd->nodes;

  for (cType curN = 1; curN <= pd->gd->N; curN++)
  {
    nodeP *np = nodes + curN;
    edgeP *pedges = np->edges;
    int cnt = np->nIdx;
    np->totalCap = 0;
    for (int ni = 0; ni < cnt; ni++)
    {
      edgeP *eh = pedges + ni;
      np->totalCap += eh->cap;
    }
  }
}

void preMark(PreprocData *pd)
{
  
  nodeP* nodes = pd->gd->nodes;

  for(cType curN =1; curN<=pd->gd->N; curN++){

      nodeP *np = nodes + curN;
      edgeP *pedges = np->edges;
      int cnt = np->nIdx;
      if (cnt < pd->K)
      {
        
        for (int ni = 0; ni < cnt; ni++)
        {
          

          edgeP *eh = pedges + ni;
          if (eh->w != MARK_CUT)
          {
            eh->rev->w = eh->w = MARK_CUT; 
            total_cut_edge++;
            total_cut_edge_lessK ++;
          }
        }
      }
  }
 

}


void collectCutLinks_step1(cType curN, PreprocData *pd)
{
  
  nodeP* nodes = pd->gd->nodes;
  assert(curN <= pd->gd->N && curN >= 1);

  
  assert((nodes + curN)->nIdx > 0);

  pd->gpnc[curN] = 1;

  short *curS = pd->gps + curN;

  assert(*curS == 0);


  long *curCV = pd->gpcv + curN;
  cType *curDep = pd->gpdep + curN;

  *curS = 1;
  *curCV = 0;
  nodeP *np = nodes + curN;
  edgeP *pedges = np->edges;
  int cnt = np->nIdx;

  if (pedges == NULL)
  {
    *curS = 2;
    return;
  }

  if (np->orderedLinks == NULL)
  {
    np->orderedLinks = (sType *)walloc(cnt + 1, sizeof(sType));
    for (int i = 0; i < cnt; i++)
    {
      np->orderedLinks[i] = i;
    }
  }
  sType *idxs = np->orderedLinks;

  

  deOrderLinkByStochas(np,pd);

  
  
    
  
  
  
  
  
  
  
  
  

  

  
  for (int ni = 0; ni < cnt; ni++)
  {
    

    edgeP *eh = pedges + idxs[ni];
    cType zn = eh->endNode;

    assert(zn != 0);
    assert(zn != curN);

    short zs = pd->gps[zn];

    assert(!(pd->gpfa[zn] == curN && zs == 2));

    
    if (zs == 1) 
    {
      
      *curCV += 1;
      pd->gpcv[zn] -= 1;

    }
    else if (zs == 0) 
    {
      
      

      pd->gpfa[zn] = curN;
      pd->gpdep[zn] = *curDep + 1;

      
      collectCutLinks_step1(zn,pd);

      
      
      assert(pd->gpdep[zn] == pd->gpdep[curN] + 1);
      
      pd->gpnc[curN] += pd->gpnc[zn];

      *curCV += pd->gpcv[zn];




		 
         
      
      
      
      
      
      
      

    }
    else
    {
      
      
      assert(pd->gpdep[curN] < pd->gpdep[zn]);
      
    }

  }

  
  
  
  
  
  
  
  

  
  
  
  
  if(*curCV == 0 && curN != gRoot){
    printf("cv=0 error: (cnt is %d) curN is %ld (gRoot is %ld), cv is %ld, dep is %ld \n",cnt,curN, gRoot,*curCV,pd->gpdep[curN]);
    assert(1==2);
  }
  
  
 

  
  
  

  
    
  assert(pd->gpdep[curN] == 0 || *curCV >0);

  
  *curS = 2;

}


cType gMe[1100000];
void collectCutLinks_step2(cType curN, PreprocData *pd)
{


  nodeP* nodes = pd->gd->nodes;

  
  cType *curDep = pd->gpdep + curN;

  cType y = pd->gpfa[curN];
  cType v = pd->gprn[curN] = pd->gprn[y];
  cType u;
  if(pd->gpcv[y] < pd->K){
    v = pd->gprn[curN] = y;
  }

  cType preRc = 0; 
  if(pd->gpcv[curN] < pd->K){
    u = curN,  preRc = pd->gprc[v], pd->gprc[v] = curN;
  }
  else{
    u = v;
  }


  nodeP *np = nodes + curN;
  edgeP *pedges = np->edges;
  int cnt = np->nIdx;
  sType *idxs = np->orderedLinks;

  for (int ni = 0; ni < cnt; ni++)
  {
    

    edgeP *eh = pedges + idxs[ni];
    cType zn = eh->endNode;
    if(cnt < pd->K){
      eh->rev->w = eh->w = MARK_CUT;
      gMe[curN] = 1;
      gMe[zn] = min(gMe[zn], pd->gd->N-1);
    }

    
    if(pd->gpdep[zn] < pd->gpdep[u] ){
      
      if(eh->w != MARK_CUT){
        eh->rev->w = eh->w = MARK_CUT; 
        total_cut_edge ++;
      }

      cType s = pd->gprn[zn];
      s = pd->gprc[s];

      if(s == zn){s = pd->gprc[s];}

      if(u > 0){ gMe[curN] = min(gMe[curN], pd->gpnc[u]);}

      gMe[zn] = min(gMe[zn], pd->gd->N - pd->gpnc[s]); 
    }
    else if(pd->gpdep[zn] == *curDep + 1){
      
      collectCutLinks_step2(zn,pd);
    }
    

  }

  if(preRc != 0){
   pd->gprc[v] = preRc;
  }
}


/*

算法整体步骤：相对之前treem算法的改进
(1)预处理
	
	sum_f = 0
	对f的每个为访问的child n:
		在遍历n前，在f上设置oh_f并置0
		在遍历计算中，每次访问祖先，除了计算mc，还更新oh_f (加上就可以)
		n返回后
			此时知道mc_n, oh_f(oh_f就是n这一支连到f的边的值，要包含f_n父子的边)
			这时计算n这一支去掉对f的mc的影响cof_n = oh_f - (mc_n-oh_f) = 2* oh_f - mc_n
				
			如果cof_n是负值 
				就加到sum_f上,即sum_f 加上 负值，sum_f指的是f的所有减少割值的子n去掉，总共减少的割值
				如果不减少，这个n就不去掉
	
	得到mc2_f = mc_f + sum_f 
				
(2)计算时: solve 和 build的时候
	
	每次向上回溯，n回溯到f时，如果cof_n是负值，说明如果要保留n这一支，目前f最优割就包含n，即此时f用于计算的mc应该取(mc2_f + -cof_n)，即加上n这一支减掉的值
  现在问题来了：
    f的mc2和访问哪一支有关系，buildAcc预处理咋做？
    这样就意味着，不同的底层上来，每个节点的mc还不一样，导致 节点段 中最小值还不一样
    buildAcc记录的是向上的，所以可以记录的呀

*/


void markCut(cType curN, PreprocData *pd)
{
  
  nodeP* nodes = pd->gd->nodes;
  assert(curN <= pd->gd->N && curN >= 1);

  
  assert((nodes + curN)->nIdx > 0);

  short *curS = pd->gps + curN;

  assert(*curS == 0);

  long *curCV = pd->gpcv + curN;
  cType *curDep = pd->gpdep + curN;

  *curS = 1;
  *curCV = 0;
  nodeP *np = nodes + curN;
  edgeP *pedges = np->edges;
  int cnt = np->nIdx;

  if (pedges == NULL)
  {
    *curS = 2;
    return;
  }

  if (np->orderedLinks == NULL)
  {
    np->orderedLinks = (sType *)walloc(cnt + 1, sizeof(sType));
    for (int i = 0; i < cnt; i++)
    {
      np->orderedLinks[i] = i;
    }
  }

  long cap;
  sType *idxs = np->orderedLinks;

  if (pd->mode == 1)
  {
    deOrderLinkByStochasCap(np,pd); 
  }
  else if (pd->mode == 2)
  {
    aOrderLinkByAvgCV(np,pd);
  }
  else if (pd->mode == 3){
    
    aOrderLinkByDegree(np,pd);
  }

  

  
  
  cType fa = curN;
  for(int i=1; i<=pd->LEVEL_NUM; i++ ){
    if(pd->gpdep[fa] == 0){
      break;
    }
    
    fa = pd->gpfa[fa];
    
    
    pd->gcutadj[curN][i] = pd->gpoh[fa]; 
  }


  for (int ni = 0; ni < cnt; ni++)
  {
    

    edgeP *eh = pedges + idxs[ni];
    cType zn = eh->endNode;

    assert(zn != 0);
    assert(zn != curN);

    short zs = pd->gps[zn];

    assert(!(pd->gpfa[zn] == curN && zs == 2));

    
    if (zs == 1) 
    {
      
      cap = eh->cap;
      *curCV += cap;
      pd->gpcv[zn] -= cap;
      pd->gpoh[zn] += cap;
    }
    else if (zs == 0) 
    {
      
      

      pd->gpfa[zn] = curN;
      pd->gpdep[zn] = *curDep + 1;

      
      markCut(zn,pd);

      
      
      assert(pd->gpdep[zn] == pd->gpdep[curN] + 1);
      
      *curCV += pd->gpcv[zn];




		 
         
      
      
      
      
      
      
      

    }
    else
    {
      
      
      assert(pd->gpdep[curN] < pd->gpdep[zn]);
      
    }

  }

  
  
  
  
  
  
  
  fa = curN;
  cType faBelow = 0; 
  int actual_ln = -1;
  for(int i=1; i<=pd->LEVEL_NUM; i++ ){
    
  
    if(pd->gpdep[fa] == 0){ 
	
      actual_ln = i-1;
      break;
    }

    fa = pd->gpfa[fa];
    
    pd->gcutadj[curN][i] = pd->gpoh[fa] - pd->gcutadj[curN][i];
    
    faBelow += pd->gcutadj[curN][i];
    
    
    pd->gcutadj[curN][i] = faBelow - (*curCV - faBelow);    
  
  }

  
  if(actual_ln < 0){
    actual_ln = pd->LEVEL_NUM;
  }

  
  

  
  
  
  

  

  /*
	下面记录curN对不同层的祖先，是否移除以及移除对祖先cv的改变
	如果actual_ln=1，说明curN的父亲就是root
  */
  
  
  
  long tempCumSum = 0;
  pd->gcutadj[curN][0] = 0; 

  
  if(actual_ln == pd->LEVEL_NUM){
	  if(pd->gcutadj[curN][actual_ln] < 0){
		  pd->gcutadjsign[curN][actual_ln] = 1; 
	  }
  }
	
  
  for (int i = (actual_ln < pd->LEVEL_NUM ? actual_ln : actual_ln - 1); i >= 0; i--) 
  {
	
	
	tempCumSum = 0;

    
    for (int ni = 0; ni < cnt; ni++)
    {
      edgeP *eh = pedges + idxs[ni];
      cType zn = eh->endNode;
      assert(zn != 0);
      assert(zn != curN);

      if (pd->gpfa[zn] != curN)
      {
        continue;
      }

      if (pd->gcutadj[zn][i + 1] < 0)
      {
		  if(curN == 820812){ 
			printf("curN %ld child zn %ld  val %ld sign %d \n",curN,zn, pd->gcutadj[zn][i+1],pd->gcutadjsign[zn][i+1]);
		  }
        tempCumSum += pd->gcutadj[zn][i + 1];
		
      }
    }

	if(curN == 820812){ 
		printf("i %d, level_cumsum[i + 1] %ld\n",i,tempCumSum);
	}
	
    if (pd->gcutadj[curN][i] < tempCumSum) 
    {
      
      
      
      pd->gcutadjsign[curN][i] = 1; 
      
    }
    else
    {


      
      
      
      pd->gcutadj[curN][i] = tempCumSum;
      pd->gcutadjsign[curN][i] = 0; 
      
    }

  }


  
  
  
  
  if(*curCV == 0 && curN != gRoot){
    printf("cv=0 error: (cnt is %d) curN is %ld (gRoot is %ld), cv is %ld, dep is %ld \n",cnt,curN, gRoot,*curCV,pd->gpdep[curN]);
    assert(1==2);
  }

  

  if(pd->mode == 1){

    for (int ni = 0; ni < cnt; ni++)
    {
      

      edgeP *eh = pedges + idxs[ni];
      cType zn = eh->endNode;
      

      assert(zn != 0);
      assert(zn != curN);
      short zs = pd->gps[zn];
    

      
      if (zs == 1 && pd->gpdep[zn] != *curDep - 1)
      {
          cType weight = eh-> w;
          if(eh->avgCV == 0){
            eh->avgCV = MAX_LONG;
          }
          eh->avgCV = min(eh->avgCV, *curCV);
          eh->w = weight+1;
          
          edgeP *reh = eh->rev;

          if(reh->avgCV == 0){
            reh->avgCV = MAX_LONG;
          }

          weight = reh-> w;
          reh->avgCV = min(reh->avgCV, *curCV);
          reh->w = weight+1;

      }
      
    }
  }

    
  assert(pd->gpdep[curN] == 0 || *curCV >0);

  
  *curS = 2;

}


/*
 (2)求解时：
      深度大的出发节点直接cv2，
      另外一个不是这个的祖先时，另外一个也可以用cv2
      [这个先不优化有点复杂]如果t是s祖先，
        如果cv2正好割掉下面，直接可以用cv2
        如果不是，也可以计算割掉对应树后这个t的割

      PS: t是s祖先，可以用cv2,t不是s的祖先也可以cv2，两个cv2都可以用
      问题是，如果t是s祖先，t的cv2对应的割并不能切断s这条支线，这样就不对了
        如果t是s祖先，其实我们算的是去掉这个支后的t的最小割
      那就简化：
        只要t不是s祖先，就可以用t的cv2

      后面的逐个逼近，需要用cv'
      
*/

void traceUp2LN(NodePropArr* pnp, cType startDep, cType curN, cType curDep, long *adj, cType LN)
{

  for (int i = 1; i <= LN; i++)
  {
    if ( curDep < (unsigned int)i || (curDep+LN) <= (startDep+i) )
    {
      break; 
    }

    if (pnp->pacc_cut_adjust_sign[curN][i] == 1)
    {
      adj[startDep+i-(curDep)] = pnp->pacc_cut_adjust[curN][i]; 
    }
  }
}




long isOptimalExclude(NodePropArr* pnp, cType v, cType v_top, cType LN_v, cType *pDep, cType *pFa, cType LN){
  
  cType minDep = 0;
  
  while(v != LN_v){
    
    
    if(pDep[v] - pDep[LN_v] <= LN && pnp->pacc_cut_adjust_sign[v][pDep[v] - pDep[LN_v]]==1){
      
      minDep = pDep[v];    
    }
    v = pFa[v];
  }

  

  
  
  
  
  
  
  /*(2) LN_s(含)及之上的祖先，则需要看不含一方的最高点，如果不含的最高点在LN_s(不含)以下，可以认为不含，因为去掉t不影响s。
        如果在之上就麻烦了，去掉t也去掉s了，这个就不好办了

  */
  
  
  if(minDep > 0){
    if(minDep <= pDep[v_top]){
      return 2; 
    }
    else{
      return 1;
    }
  }

  return 0;
}



long solveMaxFlowAccVER4(PreprocData *pd, cType root, long minCandi, NodePropArr* pnp, cType s, cType t, int SPAN_LEN, aType LN)
{
  cType *pDep = pnp->pdep;
  long *pCV = pnp->pcv;
  cType *pFa = pnp->pfa;
  
  
  
  

  
  assert(s != t);
  if (pDep[s] < pDep[t])
  {
    cType tmp = s;
    s = t;
    t = tmp;
  }

  printf("\nstart: s %ld dep %ld, t %ld dep %ld\n",s,pDep[s],t,pDep[t]);

  assert(pDep[s] >= pDep[t]);

  
  long *adjs = (long *)walloc(pd->gd->N+1, sizeof(long));
  long *adjt = (long *)walloc(pd->gd->N+1, sizeof(long));
  
  memset(adjs, 0, (pd->gd->N+1) *  sizeof(long));
  memset(adjt, 0, (pd->gd->N+1) *  sizeof(long));

  long mcv = MAX_LONG;

  long startDeps = pDep[s];
  long startDept = pDep[t];
	
	
  while(pDep[s] > pDep[t]){

    
    for (int i = 1; i <= LN; i++)
    {
      if (pnp->pacc_cut_adjust_sign[s][i] == 1)
      {
        adjs[ startDeps - (pDep[s]-i) ] = pnp->pacc_cut_adjust[s][i]; 
      }
    }

    

	
	
    mcv = min(mcv, pCV[s] + pnp->pacc_cut_adjust[s][0] - adjs[startDeps - pDep[s]]); 
	
/*	
	if(mcv < temp){
		printf("update using node s %ld , before mcv %ld\n", s,temp); 	
		printf("after mcv %ld\n", mcv); 
	}
*/
    s = pFa[s];

  }

  assert(pDep[s] == pDep[t]);

  if(s == t){
    goto end;
  }

  
  while(s != t){
    for (int i = 1; i <= LN; i++)
    {
      if (pnp->pacc_cut_adjust_sign[s][i] == 1)
      {
        adjs[ startDeps - (pDep[s]-i) ] = pnp->pacc_cut_adjust[s][i]; 
      }

      if (pnp->pacc_cut_adjust_sign[t][i] == 1)
      {
        adjs[ startDept - (pDep[t]-i) ] = pnp->pacc_cut_adjust[t][i]; 
      }


    }

    
    mcv = min(mcv, pCV[s] + pnp->pacc_cut_adjust[s][0] - adjs[startDeps - pDep[s]]); 
    mcv = min(mcv, pCV[t] + pnp->pacc_cut_adjust[t][0] - adjs[startDept - pDep[t]]); 

    s = pFa[s];
    t = pFa[t];
  }

end:
  free(adjs);
  free(adjt);
  return mcv;


}














































































































































































long parse_all_cost = 0;
void loadGraphData(PreprocData *pd){
  ;
  pd->gd = walloc(1,sizeof(GraphData));  
  printf("c\nc hi_treem version 0.9\n");
  printf("c Copyright C by nsyncw, nsyncw@gmail.com\nc\n");

  parse(&(pd->gd->N), &(pd->gd->M), &(pd->gd->nodes), &parse_all_cost);

  printf("c nodes:       %10ld\nc arcs:        %10ld\nc parse_all_cost %ld\n", pd->gd->N, pd->gd->M, parse_all_cost);
}


void initPreprocData(PreprocData *pd){
  printf("c initPreProcData\n");
  pd->rd = NULL;
  pd->SPAN_LEN = (int)(sqrt(pd->gd->N));

  pd->roots = (cType *)walloc(pd->total + 2, sizeof(cType));
  pd->allResults = walloc(pd->total+2, sizeof(NodePropArr));

  NodePropArr * allResults = pd->allResults;
  cType len = pd->gd->N + 2;

  calcuTotalCap(pd);

  int LN = pd->LEVEL_NUM+1;
  for (int i = 0; i < pd->total; i++)
  {
    
    
    
    

    allResults[i].pfa = (cType *)walloc(len, sizeof(cType));
    allResults[i].pdep = (cType *)walloc(len, sizeof(cType));
    allResults[i].pcv = (long *)walloc(len, sizeof(long));
    allResults[i].poh = (long *)walloc(len, sizeof(long));
    allResults[i].pcof = (long *)walloc(len, sizeof(long));
    allResults[i].ps = (short *)walloc(len, sizeof(short));
    allResults[i].pacc_upid = (cType *)walloc(len, sizeof(cType));
    allResults[i].pacc_upmincv = (long *)walloc(len, sizeof(long));
    allResults[i].pacc_pos_upmincv = (cType *)walloc(len, sizeof(cType));
    

    long* ptr = (long *)walloc(len*LN, sizeof(long));
    memset(ptr, 0, len * LN * sizeof(long));
    allResults[i].pacc_cut_adjust = (long **)walloc(len, sizeof(long*));
    for(int j = 0; j<len; j++){
      allResults[i].pacc_cut_adjust[j] = ptr+j*LN;
    }


    aType* ptr2 = (aType *)walloc(len*LN, sizeof(aType));
    memset(ptr2, 0, len * LN * sizeof(aType));
    allResults[i].pacc_cut_adjust_sign = (aType **)walloc(len, sizeof(aType*));
    for(int j = 0; j<len; j++){
      allResults[i].pacc_cut_adjust_sign[j] = ptr2+j*LN;
    }    

    memset(allResults[i].pfa, 0, len * sizeof(cType));
    memset(allResults[i].pdep, 0, len * sizeof(cType));
    memset(allResults[i].pcv, 0, len * sizeof(long));
    memset(allResults[i].poh, 0, len * sizeof(long));
    memset(allResults[i].pcof, 0, len * sizeof(long));
    memset(allResults[i].ps, 0, len * sizeof(short));
    memset(allResults[i].pacc_upid, 0, len * sizeof(cType));
    memset(allResults[i].pacc_upmincv, 0, len * sizeof(long));
    memset(allResults[i].pacc_pos_upmincv, 0, len * sizeof(cType));
    
    
    

  }  
}



void preProc(PreprocData *pd){
  printf("c preProc \n");
  double tm;
  double totalProcTime = 0;
  NodePropArr *allResults = pd->allResults;
  
  cType root;

  cType len = pd->gd->N + 2;
  long *apply_adj = (long *)walloc(len, sizeof(long));
  cType *depth_map = (cType *)walloc(len, sizeof(cType));

  preMark(pd);
  printf("c total_cut_edge %ld after premark \n",total_cut_edge);

  memset(gMe, 100000, len*sizeof(cType));



  for (int ipass = 0; ipass < pd->total; ipass++)
  {

    if(pd->rd != NULL){
      free(pd->rd->randNums);
      free(pd->rd);
      pd->rd = NULL;
    }
    pd->rd = initrand(pd->gd->M*2);    
    
    pd->gpfa = allResults[ipass].pfa;
    pd->gpdep = allResults[ipass].pdep;
    pd->gpcv = allResults[ipass].pcv;
    pd->gpoh = allResults[ipass].poh;
    pd->gpnc = allResults[ipass].pcof;
    pd->gps = allResults[ipass].ps;
    pd->gprc = allResults[ipass].pacc_upid; 
    pd->gpaccmcv = allResults[ipass].pacc_upmincv; 
    pd->gprn = allResults[ipass].pacc_pos_upmincv; 
    pd->gcutadj = allResults[ipass].pacc_cut_adjust;
    pd->gcutadjsign = allResults[ipass].pacc_cut_adjust_sign;

    if (pd->P == 300)
    {
      pd->mode = 3;
    }
    else
    {
      pd->mode = ipass < pd->P * pd->total / 100 ? 1 : 2;
    }
	

    
    
    
    
    
    
    root = (rand() % pd->gd->N)+1;
    
    pd->roots[ipass] = root;
    pd->gpdep[root] = 0;
    printf("pass %d before markCut: root is %ld\n",ipass,root);
    fflush(stdout);
    
    
    tm = timer();
    gRoot = root;
    
    collectCutLinks_step1(root,pd);

    pd->gprn[root] = 0;

    collectCutLinks_step2(root,pd);

    pd->gpcv[root] = MAX_LONG;
    pd->gpoh[root] = MAX_LONG;


    
    /*
    printf("c before buildAcc\n");
	  fflush(stdout);

    memset(apply_adj, 0, len *  sizeof(long));
    memset(depth_map, 0, len *  sizeof(cType));
    buildAcc(pd, root, root, 0,MAX_LONG, apply_adj,depth_map);
    printf("c after buildAcc\n");
    fflush(stdout);
    */
    

    totalProcTime += timer() - tm;
    printf("c proctime for onepass: %10.06f\n", timer() - tm);
    if (ipass % 10 == 0)
    {
      printf("c the %d passes\n", ipass);
    }

    
    
    
    
    
    
    
    
    

  }

  free(apply_adj);
  free(depth_map);

  printf("c preprocess times %10.6f\n", totalProcTime);
  printf("total_cut_edge %ld \n", total_cut_edge);

}


void HeapAdjustDown2(sType *idx, edgeP ** edges ,int start,int end)  
{  
    sType tempIdx = idx[start];  

    int i = 2*start+1;      
    
    
    
    while(i<=end)  
    {  
        if(i+1<=end && edges[idx[i+1]]->cost > edges[idx[i]]->cost )    
            i++;  

        if(edges[idx[i]]->cost <= edges[tempIdx]->cost )   
            break;  

        idx[start] = idx[i];

        start = i;  
        i = 2*start+1;  
    }  

    idx[start] = tempIdx;  
}  
  
void HeapSort2(sType *idx, edgeP ** edges, int len)  
{  

    int i;  
    for(i=(len-1)/2;i>=0;i--){  
        HeapAdjustDown2(idx,edges,i,len-1);  
    }

    for(i=len-1;i>0;i--)
    {  
        
        sType temp = idx[i];  
        idx[i] = idx[0];  
        idx[0] = temp;  

        HeapAdjustDown2(idx,edges,0,i-1);  
    }  

} 



edgeP** estack = NULL;
int stackTop = 0;
void buildMCT_cost0Tree(cType curN, PreprocData *pd, cType treeId)
{


  nodeP* nodes = pd->gd->nodes;
  nodeP *np = nodes + curN;
  edgeP *pedges = np->edges;
  int cnt = np->nIdx;
  pd->gpaccmcv[curN] = treeId;

  for (int ni = 0; ni < cnt; ni++)
  {
    

    edgeP *eh = pedges + ni;
    cType zn = eh->endNode;

    if(eh->w == MARK_CUT){ 
      
      if(pd->gpaccmcv[zn] > 0 && pd->gpaccmcv[zn] != treeId){
        
        estack[stackTop++] = eh;
        
      }
    }
    else{
      if(pd->gpaccmcv[zn] == 0){
        
        buildMCT_cost0Tree(zn,pd,treeId);
      }
      else{
        
      }
    }

 
    

  }

}

sType* idxs = NULL;

   long usedLinkCount = 0;
   long usedLinkCountDegree = 0; 
   long totalLinkCountDegree = 0; 
   long totalCost = 0;


void buildMCT(PreprocData *pd){
    
    
    long M = pd->gd->M;
    long N = pd->gd->N;
    long len = M + 2;
  

    estack = (edgeP**)walloc(len, sizeof(edgeP*)); 
    memset(estack, 0, len * sizeof(edgeP*));
    stackTop = 0;

    idxs = (sType*)walloc(len,sizeof(sType));
    memset(idxs,0,len*sizeof(sType));  
    
    
    
    
        
    
    for(int i=1; i<=N; i++){
      if(pd->gpaccmcv[i] > 0){
        
      }
      else{
        buildMCT_cost0Tree(i,pd, i); 
      }
    }
    
    printf("stackTop is %d \n",stackTop);
    
    for(int i=0; i<stackTop; i++){
      idxs[i] = i;
    }

    
    
    HeapSort2(idxs, estack, stackTop);  
    
    /*
      先根据映射树，确认两边节点所属于的最终的ID
        此时同时更新为id直接映射过去
      如果属于同个ID,则continue下一个
        否则可以合并树，得到仍然是树

        

    */    

    for(int k=0; k<stackTop; k++){
      int i = idxs[k]; 
      
      

      cType n1 = estack[i]->endNode;
      cType n2 = estack[i]->rev->endNode;
      cType t1 = n1;
      

      while(t1 != pd->gpaccmcv[t1]){
        assert(t1 > 0);
        t1 = pd->gpaccmcv[t1];
      }

      cType t2 = n2;
      while(t2 != pd->gpaccmcv[t2]){
        assert(t2 > 0);
        t2 = pd->gpaccmcv[t2];
      }

      if(pd->gpaccmcv[n1] > 0) { pd->gpaccmcv[pd->gpaccmcv[n1]] = t1;}
      if(pd->gpaccmcv[n2] > 0) { pd->gpaccmcv[pd->gpaccmcv[n2]] = t1;}

      pd->gpaccmcv[n1] = t1;
      pd->gpaccmcv[n2] = t2;

      if((pd->gd->nodes + n1)->nIdx < pd->K || (pd->gd->nodes + n2)->nIdx < pd->K){
        totalLinkCountDegree++;
      }

      if(t1 != t2){
        pd->gpaccmcv[t2] = t1;    
        usedLinkCount ++;
        totalCost += estack[i] ->cost;   

        if(estack[i]->w != (MARK_CUT_SHIELD + 1)){ 
          estack[i]->w = estack[i]->rev->w = MARK_CUT_SHIELD;   
        }

        
        if((pd->gd->nodes + n1)->nIdx < pd->K || (pd->gd->nodes + n2)->nIdx < pd->K){
          usedLinkCountDegree++;
        }
      }
      else{
        
      }

    }

    printf("c K %ld, total, chosen, cost is %ld(%ld) %ld(%ld) %ld;  all lessK neighbor edge %ld\n",pd->K,total_cut_edge,totalLinkCountDegree,usedLinkCount,usedLinkCountDegree,totalCost,total_cut_edge_lessK);

}



void printGraphWithShieldedLinks(PreprocData *pd){
  printf("p max %ld %ld\n",pd->gd->N, pd->gd->M);
  int cap = 1;
  for(cType curN = 1; curN<=pd->gd->N; curN++){
    
    for(int ni=0; ni<(pd->gd->nodes + curN)->nIdx;ni++){
      cType n2 = (pd->gd->nodes + curN)->edges[ni].endNode;
      if(n2 > curN){
        cap = 1;
        if((pd->gd->nodes + curN)->edges[ni].w == MARK_CUT_SHIELD){
          
          cap = 100;
        }
        printf("a %ld %ld %d\n",curN,n2,cap);
      }      
    }
  }
}

/*

===========================random  n1000000dx.inp
=====d=2
total_cut_edge 1471331
stackTop is 1471331
c total, chosen, cost is 1471331(1441581) 675677(665856) 2460821;  all lessK neighbor edge 1441581
the final total time: 239.535999 (0.128000 233.223993 6.184006)


=====d=3
total_cut_edge 886294
stackTop is 886294
c total, chosen, cost is 886294(882600) 290249(289258) 906531;  all lessK neighbor edge 882600
the final total time: 325.743995 (0.180000 322.895996 2.667999)



=====d=4
total_cut_edge 321770
stackTop is 321770
c total, chosen, cost is 321770(321542) 93404(93344) 274425;  all lessK neighbor edge 321542
the final total time: 403.859996 (0.220000 402.388012 1.251984)



========================scale free  n1000000sfxc1.inp
===d=2

total_cut_edge 1710193
stackTop is 1710193
c total, chosen, cost is 1710193(1710193) 848310(848310) 3160140;  all lessK neighbor edge 1710193
the final total time: 207.780006 (0.132000 202.936000 4.712006)




===d=3
total_cut_edge 1736095
stackTop is 1736095
c total, chosen, cost is 1736095(1736095) 569082(569082) 1661318;  all lessK neighbor edge 1736095
the final total time: 305.743993 (0.148000 301.103988 4.492004)



===d=4
total_cut_edge 1240152
stackTop is 1240152
c total, chosen, cost is 1240152(1240152) 310038(310038) 785834;  all lessK neighbor edge 1240152
the final total time: 379.723990 (0.084000 377.363989 2.276001)





*/


long harrIdx[5][5][10];
edgeP* harrLink[5][5][10][500000];
void heuristicPreprocess(PreprocData *pd, int gType, int gDegree){
  long toCheckTotalCost = 0;
  long toCheckTotalLink = 0;
  long totalTmp1 = 0;
  long totalCostArr[2][3];
  totalCostArr[0][0] = 2460821;
  totalCostArr[0][1] = 906531;
  totalCostArr[0][2] = 274425;
  totalCostArr[1][0] = 3160140;
  totalCostArr[1][1] = 1661318;
  totalCostArr[1][2] = 785834;
  totalCost = totalCostArr[gType-1][gDegree-2];
  memset(harrIdx, 0, 250 *  sizeof(long));
  memset(harrLink, 0, 250 * 100000*  sizeof(edgeP*));
  for(cType curN = 1; curN<=pd->gd->N; curN++){
    
    int curCnt = (pd->gd->nodes + curN)->nIdx;
    totalTmp1 += curCnt;
    
    for(int ni=0; ni<curCnt; ni++){
      toCheckTotalLink ++;      
  
      cType n2 = (pd->gd->nodes + curN)->edges[ni].endNode;
      
      if(n2 > pd->gd->N){
        printf("error n2 is %ld\n",n2);
        exit(0);
      }

      if(n2 > curN){
        int b = (pd->gd->nodes + n2)->nIdx;
        edgeP* ep = (pd->gd->nodes + curN)->edges+ni;
        ep->w = ep->rev->w = 0;

        int a = curCnt;
        if(a > b){
          int tmp = a;
          a = b;
          b = tmp;
          ep = ep->rev;
        }

        a = min(a, pd->K);
        b = min(b, pd->K);
        long cost = ep->cost;

        toCheckTotalCost += cost;
        

        harrLink[a-1][b-1][cost-1][harrIdx[a-1][b-1][cost-1]] = ep;
        harrIdx[a-1][b-1][cost-1] += 1;

      }      
    }

  }

  printf("c toCheckTotalCost %ld, toCheckTotalLink %ld, totalTmp1 %ld\n",toCheckTotalCost,toCheckTotalLink,totalTmp1);

  long test_all_cost = 0;
  for (int i = 1; i <= pd->K; i++)
  {
    for (int j = i; j <= pd->K; j++)
    {
      for (int k = 0; k < 10; k++)
      {
        test_all_cost += harrIdx[i-1][j-1][k]*(k+1);
      }
    }
  }

  printf("c test_all_cost %ld, parse_all_cost %ld\n", test_all_cost, parse_all_cost);
  if(test_all_cost != parse_all_cost || toCheckTotalLink/2 != pd->gd->M){
    printf("test_all_cost != parse_all_cost or toCheckTotalLink != M\n");
    fflush(stdout);
    exit(0);
  }

}

    
    /*
    (1,2)(2,2)(2,3)(3,3)
    */






































/*
  无非就是可能选大一点的
  随机咋实现？排序不能用，太耗时
  (1)范围就是<5的所有，先看代价是否合适？
  (2)怎么挑选？随机去除？
  随机+交换方法，一个索引记录队尾，随机结果如果不是队尾和队尾交换

*/

/*

*/






































































































































































































































































void calcuStochasPairs(int numOfPairs, PreprocData *pd){
  double totalTime = 0;
  long mv = MAX_LONG;

  double curTime = 0;
  cType ns, nt;

  if(pd->rd != NULL){
    free(pd->rd);
    pd->rd = NULL;
  }
  pd->rd = initrand(pd->gd->M*2);  
  ns = nt  = 1;  

  for (int ipair = 0; ipair < numOfPairs;)
  {


    if (ns != nt)
    {
      
      mv = min((pd->gd->nodes+ns)->totalCap, (pd->gd->nodes+nt)->totalCap);
      
      curTime = timer();
      for (int j = 0; j < pd->total; j++)
      {
        cType root = pd->roots[j];
        
        long tmp = solveMaxFlowAccVER4(pd, root,mv, pd->allResults+j, ns, nt,pd->SPAN_LEN,pd->LEVEL_NUM);
        
        if (mv > tmp)
        {
          mv = tmp;
        }
        
      }
      curTime = timer() - curTime;
      totalTime += curTime;
      ipair++;
      printf("c hi_treem_res(n,s,mflow,tm) %lu %lu %12.01f %12.06f1\n", ns, nt, 1.0 * mv, curTime);
    }
	
	    
    ns = 1 + ((mrand(pd->rd) * mrand(pd->rd)) % (pd->gd->N));
    nt = 1 + ((mrand(pd->rd) * mrand(pd->rd)) % (pd->gd->N));
  }

  printf("c run ok! average time %10.6f\n", totalTime / numOfPairs);


}










void connTwoNodes(PreprocData *pd, cType head, cType tail){

    if(head==999999 || tail == 999999) {
      printf("conn %ld to %ld\n",head,tail);
    }
    nodeP* np = pd->gd->nodes +head;
    adjustLinksIfNecessary(np);

    edgeP* newedge = np->edges+np->nIdx;
    newedge->endNode = tail;
    newedge->w = 1;
    newedge->cap = 1;
    np->nIdx ++;

    if(np->orderedLinks != NULL) {
      free(np->orderedLinks);
      np->orderedLinks = (sType *)walloc(np->nIdx + 1, sizeof(sType));
    }
    edgeP* newedge_prev = newedge;
    
    np = pd->gd->nodes + tail;
    adjustLinksIfNecessary(np);

    newedge = np->edges + np->nIdx;
    newedge->endNode = head;
    newedge->w = 1;
    newedge->cap = 1;
    np->nIdx ++;

    if(np->orderedLinks != NULL) {
      free(np->orderedLinks);
      np->orderedLinks = (sType *)walloc(np->nIdx + 1, sizeof(sType));
    }

    newedge_prev->rev = newedge;
    newedge->rev = newedge_prev;  

    pd->gd->M++;

}



cType _justConCandiPair(PreprocData *pd, cType n1, cType n2){

  if(n1==999999 && n2==1000000){ printf(" here try conn 999999 with 1000000\n");}
  if(n2==999999 && n1==1000000){ printf(" here try conn 999999 with 1000000\n");}

  if(n1 > 0 && n2 > 0 && n1 != n2 && n1 <= pd->gd->N && n2 <=pd->gd->N)
  {  
      
      
      
      
      for (int ni = 0; ni < (pd->gd->nodes+n2)->nIdx; ni++)
      {
        

        cType zn = ((pd->gd->nodes+n2)->edges + ni)->endNode;
        if (n2 > 999997) printf("check n2 %ld n1 %ld zn %ld\n",n2,n1,zn);
        if(zn == n1){
          return 0;
        }

      }


      connTwoNodes(pd, n1, n2);
      return 1;
  }
  else{
      return 0;
  }
  
}

cType heuristicSmallwDegConn(PreprocData *pd){
	
	cType cnt = 0;
  cType nodeStack[1100000];
  cType nodeStackIdx = 0;




  cType totalSmallDegCount = 0;
  cType totalLessUpBoundCount = 0;

	for(cType curN=1; curN <= pd->gd->N; curN++){
    if(pd->gd->nodes[curN].nIdx < pd->K){
      gMe[curN] = 1;
      totalSmallDegCount++;
    }
    else{
      gMe[curN] = 100;
    }

    

		if(gMe[curN] <= pd->meUpbound ){
    
      
      totalLessUpBoundCount++;
      for (int ni = 0; ni < (pd->gd->nodes+curN)->nIdx; ni++)
      {
        

        edgeP* eh = ((pd->gd->nodes+curN)->edges + ni);
        
        eh->w = eh->rev->w = MARK_CUT_SHIELD+1;

          

        
        

          
        
        
        
        

      }

			if(nodeStackIdx == 0){ nodeStack[nodeStackIdx++] = curN;}
      else{
        
				cType ret = _justConCandiPair(pd,nodeStack[nodeStackIdx-1],curN);
        if(ret > 0){
          nodeStackIdx--;
          cnt ++;
        }
				else{
          
          nodeStack[nodeStackIdx++] = curN;
        }
			}
      assert(nodeStackIdx < 1100000);
		}
	}

    cType rnodeStack[1100000];
    cType rnodeStackIdx = 0;

	if(nodeStackIdx > 0){



    for(cType i = nodeStackIdx -1; i>0; i--){
      if(nodeStack[i] == 0) { continue; }
      
      cType j=0;
      for(; j<i; j++){
        if(nodeStack[j] == 0) { continue; }
    		cType ret = _justConCandiPair(pd, nodeStack[i], nodeStack[j]);	
        if(ret > 0){
          nodeStack[j] = 0;
          cnt ++;
          nodeStack[i] = 0;
          break;
        }
      }

      if(j<i){
        
      }
      else{
        
        rnodeStack[rnodeStackIdx++] = i;
      }
    }

    if(nodeStack[0] > 0){
      rnodeStack[rnodeStackIdx++] = nodeStack[0];
    }

    assert(rnodeStackIdx < 1100000);
	}

  printf("c %ld paired added edges using given small deg nodes\n", cnt);

  if(rnodeStackIdx > 0){
    while(rnodeStackIdx > 0){
      cType n = rnodeStack[--rnodeStackIdx];

 		  while(_justConCandiPair(pd, n ,1+rand()%pd->gd->N) <=0){
        
      }	
      cnt++;

    }
   
  }




  printf("c small deg count: %ld, paired count less than upbound count %ld \n",totalSmallDegCount,totalLessUpBoundCount);




	
	return cnt;
}

int isNeighbor(PreprocData *pd, cType n1, cType n2){
      for (int ni = 0; ni < (pd->gd->nodes+n2)->nIdx; ni++)
      {
        
        cType zn = ((pd->gd->nodes+n2)->edges + ni)->endNode;
        if(zn == n1){ return 1;}  

      }  

      return 0;
}


cType cycleConnShieldLinks(PreprocData *pd){
  cType segArr[1000000][2]; 
  int segArrIdx = 0;
  sType nodeSign[1100000];
  memset(nodeSign, 0,1100000*sizeof(sType));

  sType tcnt = 0;
  
  while(stackTop >0){
    edgeP* eh = estack[--stackTop];
    if(eh->w != MARK_CUT_SHIELD){ continue; }

    tcnt ++;

    eh->w = MARK_CUT_SHIELD+1;
    eh->rev->w = MARK_CUT_SHIELD+1;

    cType curN = eh->endNode;
    
    while(1==1){
      int i;
      for(i=0; i<pd->gd->nodes[curN].nIdx; i++){
        if(pd->gd->nodes[curN].edges[i].w == MARK_CUT_SHIELD){
          
          tcnt++;
          pd->gd->nodes[curN].edges[i].w = MARK_CUT_SHIELD +1;
          pd->gd->nodes[curN].edges[i].rev->w = MARK_CUT_SHIELD +1;
          
          curN = pd->gd->nodes[curN].edges[i].endNode;
          break;
        }
      }   

      if(i >= pd->gd->nodes[curN].nIdx ){
        
        break;
      }
      else{
        
      }

    }

    cType end1 = curN;
    
    curN = eh->rev->endNode;
    while(1==1){
      int i;
      for(i=0; i<pd->gd->nodes[curN].nIdx; i++){
        if(pd->gd->nodes[curN].edges[i].w == MARK_CUT_SHIELD){
          
          tcnt++;
          pd->gd->nodes[curN].edges[i].w = MARK_CUT_SHIELD +1;
          pd->gd->nodes[curN].edges[i].rev->w = MARK_CUT_SHIELD +1;
          curN = pd->gd->nodes[curN].edges[i].endNode;
          break;
        }
      }   

      if(i >= pd->gd->nodes[curN].nIdx ){
        
        break;
      }
      else{
        
      }

    }

    cType end2 = curN;

    
    
    

    

    segArr[segArrIdx][0] = end1;
    segArr[segArrIdx][1] = end2;
    segArrIdx++;

    
    


  }

  assert(tcnt == usedLinkCount);
  if(segArrIdx < 2){
    printf("c warning: segArrIdx %d <2 \n", segArrIdx);
  }


  

  
  cType fCnt = 0;
  

  for(int i = 0 ; i < segArrIdx; i++){
    cType n1 = segArr[i][0];
    cType n2 = segArr[i][1];
    cType outNode;

    while(1==1){
      outNode = 1 + rand()%pd->gd->N;
      if(outNode == n1 || outNode == n2 || isNeighbor(pd, n1, outNode) == 1 || isNeighbor(pd,n2,outNode) == 1){
        continue;
      }

      break;
    }

    connTwoNodes(pd,n1,outNode );
    connTwoNodes(pd,n2,outNode );

    fCnt +=2;
  }

  
  
  
  
  

  assert(fCnt == segArrIdx*2);

  return fCnt;

}



void outputGraph(PreprocData *pd){
  char fname[1024];
  sprintf(fname,"n%ldm%ldc1.inp",pd->gd->N, pd->gd->M);
  FILE *file = fopen(fname,"w");
  fprintf(file,"p max %ld %ld\n",pd->gd->N, pd->gd->M); 
  for(cType curN=1; curN < pd->gd->N; curN++){
    for(int i=0; i<pd->gd->nodes[curN].nIdx; i++){
      if(pd->gd->nodes[curN].edges[i].endNode > curN){
        
        fprintf(file,"a %ld %ld 1\n",curN,pd->gd->nodes[curN].edges[i].endNode);
      }
    }
  }

  fclose(file);
	printf("save to file %s\n",fname);
}
