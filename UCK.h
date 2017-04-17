#ifndef _MCMF_H_
#define _MCMF_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <utility>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>
#include <set>
#include <vector>
#include <stack>
#include <map>
#include <queue>
#include <deque>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>
#include <sys/types.h>
#include <unistd.h>
#include <functional>
#include <fstream>
#include <cassert>

using namespace std;

typedef long long int excessType,priceType;

class MCMF_YCW
{
        public:
                class ARC;
                //节点
                class NODE {
                        public:
                                excessType excessNodeVar;
                                priceType priceNodeVar;
                                ARC *firstNodeVar,*currentNodeVar,*suspendedNodeVar;
                                NODE *qNextNodeVar,*bNextNodeVar,*bPrevNodeVar;
                                long rankNodeVar, inpNodeVar;

                        public:
                                void set_excess( excessType excess) { excessNodeVar = excess; }
                                void dec_excess( long delta) { excessNodeVar -= delta; }
                                void inc_excess( long delta) { excessNodeVar += delta; }
                                void set_price( priceType price) { priceNodeVar = price; }
                                void dec_price( long delta) { priceNodeVar -= delta; }
                                void set_first( ARC *first) { firstNodeVar = first; }
                                void set_current( ARC *current) { currentNodeVar = current; }
                                void inc_current() { currentNodeVar ++; }
                                void set_suspended( ARC *suspended) { suspendedNodeVar = suspended; }
                                void set_q_next( NODE *q_next) { qNextNodeVar = q_next; }
                                void set_b_next( NODE *b_next) { bNextNodeVar = b_next; }
                                void set_b_prev( NODE *b_prev) { bPrevNodeVar = b_prev; }
                                void set_rank( long rank) { rankNodeVar = rank; }
                                void set_inp( long inp) { inpNodeVar = inp; }
                                excessType excess() { return excessNodeVar; }
                                priceType price() { return priceNodeVar; }
                                ARC *first() { return firstNodeVar; }
                                void dec_first() { firstNodeVar --; }
                                void inc_first() { firstNodeVar ++; }
                                ARC *current() { return currentNodeVar; }
                                ARC *suspended() { return suspendedNodeVar; }
                                NODE *q_next() { return qNextNodeVar; }
                                NODE *b_next() { return bNextNodeVar; }
                                NODE *b_prev() { return bPrevNodeVar; }
                                long rank() { return rankNodeVar; }
                                long inp() { return inpNodeVar; }
                };

                //弧定义
                class ARC {
                        public:
                                long rezCapacityARCVar; // 剩余容量
                                priceType costARCVar;       // 费用
                                NODE *headARCVar;
                                ARC *sisterARCVar;      // 相对弧
                        public:
                                void set_rez_capacity( long rez_capacity) { rezCapacityARCVar = rez_capacity; }
                                void dec_rez_capacity( long delta) { rezCapacityARCVar -= delta; }
                                void inc_rez_capacity( long delta) { rezCapacityARCVar += delta; }
                                void set_cost( priceType cost) { costARCVar = cost; }
                                void multiply_cost( priceType mult) { costARCVar *= mult; }
                                void set_head( NODE *head) { headARCVar = head; }
                                void set_sister( ARC *sister) { sisterARCVar = sister; }
                                long rez_capacity() { return rezCapacityARCVar; }
                                priceType cost() { return costARCVar; }
                                NODE *head() { return headARCVar; }
                                ARC *sister() { return sisterARCVar; }
                };



                //
                class BUCKET {
                        private:
                                NODE *pFirstBucketVar;
                        public:
                                BUCKET( NODE *p_first) : pFirstBucketVar(p_first) {}
                                BUCKET() {}
                                void set_p_first( NODE *p_first) { pFirstBucketVar = p_first; }
                                NODE *p_first() { return pFirstBucketVar; }
                };

                //计算总费用需要的属性
                int netStates,consumeStates,totalNeed;
                map<int,pair<int,int> > serverLevel;//保存服务器档次信息
                map<int, int> serverAllPrice;
                map<int, int> getServerAllPrice();
                vector<int> netStateDeployPrice;//网络节点部署成本
                bool isOutputResult = false;

                long nodeNumYCWVar,arcNumYCWVar,*capArrayYCWVar;   // 容量数组
                NODE *nodesArrayYCWVar,*sentinelNode,*excqFirst,*excqLast;
                ARC *arcsArrayYCWVar,*sentinelArc;

                BUCKET *bucketsArrayYCWVar,*lBucket;
                long _linf;
                int timeForPriceIn;

                priceType epsilonYCWVar,_dn,priceMinYCWVar,_mmc;
                double fScaleYCWVar,cutOffFactor,cutOn,cutOff;
                excessType totalExcess;

                int flagPrice,flagUpdt,sncMaxYCWVar;

                ARC dArcYCWVar;
                NODE dOneNode,*_dummy_node,*dTwoNode;

                long nRelYCWVar,nRefYCWVar,nSrcYCWVar,nPushYCWVar,nRelabelYCWVar,nDischargeYCWVar,nRefineYCWVar,nUpdateYCWVar,nScanYCWVar,nPrscanYCWVar,nPrscanOneYCWVar;
                long nPrscanTwoYCWVar,nBadPricein,nBadRelabelYCWVar,nPrefineYCWVar;

                bool noZeroCycles,checkSolution,compDualsYCWVar,costRestartYCWVar,printAnsYCWVar;
                long long int *nodeBalanceYCWVar;

                long nodeMinYCWVar,nodeMaxYCWVar,*arcFirstYCWVar,*arcTailYCWVar,posCurrent;
                ARC *arcCurrent,*arcNew,*arcTmp;
                priceType maxCost;
                excessType totalPYCWVar,totalNYCWVar;

                NODE *iNode,*jNode;

                //需要计算总费用的构造函数
                MCMF_YCW( long num_nodes, long num_arcs,int netStates,int consumeStates,int totalNeed,map<int,pair<int,int> > serverLevel,vector<int> netStateDeployPrice)
                    :flagPrice(0),flagUpdt(0),nPushYCWVar(0),nRelabelYCWVar(0),nDischargeYCWVar(0),nRefineYCWVar(0),nUpdateYCWVar(0),nScanYCWVar(0),nPrscanYCWVar(0),
                nPrscanOneYCWVar(0),nPrscanTwoYCWVar(0),nBadPricein(0),nBadRelabelYCWVar(0),nPrefineYCWVar(0),noZeroCycles(false),checkSolution(false),compDualsYCWVar(false),
                costRestartYCWVar(false),printAnsYCWVar(true){

                        this->nodeNumYCWVar = num_nodes;
                        this->arcNumYCWVar = num_arcs;

                        this->netStates = netStates;
                        this->consumeStates = consumeStates;
                        this->totalNeed = totalNeed;

                        //服务器档次和网络节点部署费用
                        this->serverLevel = serverLevel;
                        this->netStateDeployPrice = netStateDeployPrice;


                        allocate_arrays();
                }
                ~MCMF_YCW() {}

                void allocate_arrays();
                void deallocate_arrays();
                void set_arc( long tail_node_id, long head_node_id,
                                          long low_bound, long up_bound, priceType cost);
                void set_supply_demand_of_node( long id, long excess);
                void pre_processing();
                void cs2_initialize();
                void up_node_scan( NODE *i);
                void price_update();
                int relabel( NODE *i);
                void discharge( NODE *i);
                int price_in();
                int refine();
                int price_refine();
                void compute_prices();
                void price_out();
                int update_epsilon();
                int check_feas();
                int check_cs();
                int check_eps_opt();
                void init_solution();
                void cs_cost_reinit();
                void cs2_cost_restart( double *objective_cost);
                void print_solution(string& result);
                void print_graph();
                void finishup( double *objective_cost);
                void cs2( double *objective_cost);
                int greenTea();
                int addServerAndDeployPrice(double *cost);
                pair<int,int> determineDC(int serverOutput);
                //确定存储每个服务器的档次
                void storeServerGrade();

                void increase_flow( NODE *i, NODE *j, ARC *a, long df) {
                        i->dec_excess( df);
                        j->inc_excess( df);
                        a->dec_rez_capacity( df);
                        a->sister()->inc_rez_capacity( df);
                }
                bool time_for_update() {
                        return ( nRelYCWVar > nodeNumYCWVar * 0.4 + nSrcYCWVar * 30);
                }
                void reset_excess_q() {
                        for ( ; excqFirst != NULL; excqFirst = excqLast ) {
                                excqLast = excqFirst->q_next();
                                excqFirst->set_q_next( sentinelNode);
                        }
                }
                bool out_of_excess_q( NODE *i) { return ( i->q_next() == sentinelNode); }
                bool empty_excess_q() { return ( excqFirst == NULL); }
                bool nonempty_excess_q() { return ( excqFirst != NULL); }
                void insert_to_excess_q( NODE *i) {
                        if ( nonempty_excess_q() ) {
                                excqLast->set_q_next( i);
                        } else {
                                excqFirst = i;
                        }
                        i->set_q_next( NULL);
                        excqLast = i;
                }
                void insert_to_front_excess_q( NODE *i) {
                        if ( empty_excess_q() ) {
                                excqLast = i;
                        }
                        i->set_q_next( excqFirst);
                        excqFirst = i;
                }
                void remove_from_excess_q( NODE *i) {
                        i = excqFirst;
                        excqFirst = i->q_next();
                        i->set_q_next( sentinelNode);
                }
                bool empty_stackq() { return empty_excess_q(); }
                bool nonempty_stackq() { return nonempty_excess_q(); }
                void reset_stackq() { reset_excess_q(); }
                void stackq_push( NODE *i) {
                        i->set_q_next( excqFirst);
                        excqFirst = i;
                }
                void stackq_pop( NODE *i) {
                        remove_from_excess_q( i);
                }
                void reset_bucket( BUCKET *b) { b->set_p_first( dTwoNode); }
                bool nonempty_bucket( BUCKET *b) { return ( (b->p_first()) != dTwoNode); }
                void insert_to_bucket( NODE *i, BUCKET *b) {
                        i->set_b_next( b->p_first() );
                        b->p_first()->set_b_prev( i);
                        b->set_p_first( i);
                }
                void get_from_bucket( NODE *i, BUCKET *b) {
                        i = b->p_first();
                        b->set_p_first( i->b_next());
                }
                void remove_from_bucket( NODE *i, BUCKET *b) {
                        if ( i == b->p_first() ) {
                                b->set_p_first( i->b_next());
                        } else {
                                i->b_prev()->set_b_next( i->b_next());
                                i->b_next()->set_b_prev( i->b_prev());
                        }
                }
                void update_cut_off() {
                        if ( nBadPricein + nBadRelabelYCWVar == 0) {
                                cutOffFactor = pow( (double)nodeNumYCWVar, 0.75 );
                                cutOffFactor = cutOffFactor > 12 ? cutOffFactor:12;
                                cutOff = cutOffFactor * epsilonYCWVar;
                                cutOn = cutOff * 0.8;
                        } else {
                                cutOffFactor *= 4;
                                cutOff = cutOffFactor * epsilonYCWVar;
                                cutOn = cutOff * 0.8;
                        }
                }
                void exchange( ARC *a, ARC *b) {
                        if ( a != b) {
                                ARC *sa = a->sister();
                                ARC *sb = b->sister();
                                long d_cap;

                                dArcYCWVar.set_rez_capacity( a->rez_capacity());
                                dArcYCWVar.set_cost( a->cost());
                                dArcYCWVar.set_head( a->head());

                                a->set_rez_capacity( b->rez_capacity());
                                a->set_cost( b->cost());
                                a->set_head( b->head());

                                b->set_rez_capacity( dArcYCWVar.rez_capacity());
                                b->set_cost( dArcYCWVar.cost());
                                b->set_head( dArcYCWVar.head());

                                if ( a != sb) {
                                        b->set_sister( sa);
                                        a->set_sister( sb);
                                        sa->set_sister( b);
                                        sb->set_sister( a);
                                }

                                d_cap = capArrayYCWVar[ a - arcsArrayYCWVar];
                                capArrayYCWVar[ a - arcsArrayYCWVar] = capArrayYCWVar[ b - arcsArrayYCWVar];
                                capArrayYCWVar[ b - arcsArrayYCWVar] = d_cap;
                        }
                }
};
#endif
