#include "mcmf.h"

void MCMF_YCW::allocate_arrays() {

        nodesArrayYCWVar = (NODE*) calloc ( nodeNumYCWVar+2,   sizeof(NODE) );
        arcsArrayYCWVar = (ARC*)  calloc ( 2*arcNumYCWVar+1, sizeof(ARC) );
        capArrayYCWVar = (long*) calloc ( 2*arcNumYCWVar,   sizeof(long) );

        arcTailYCWVar = (long*) calloc ( 2*arcNumYCWVar,   sizeof(long) );
        arcFirstYCWVar = (long*) calloc ( nodeNumYCWVar+2,   sizeof(long) );

        for ( NODE *in = nodesArrayYCWVar; in <= nodesArrayYCWVar + nodeNumYCWVar; in ++ ) {
                in->set_excess( 0);
        }

        // (2) resets;
        posCurrent = 0;
        arcCurrent = arcsArrayYCWVar; 
        nodeMaxYCWVar = 0;
        nodeMinYCWVar = nodeNumYCWVar;
        maxCost = 0;
        totalPYCWVar = totalNYCWVar = 0;
}

void MCMF_YCW::deallocate_arrays() {
        if ( arcsArrayYCWVar) free ( arcsArrayYCWVar );
        if ( dTwoNode) delete dTwoNode;
        if ( capArrayYCWVar) free ( capArrayYCWVar );
        if ( bucketsArrayYCWVar) free ( bucketsArrayYCWVar );
        if ( checkSolution == true) free ( nodeBalanceYCWVar );
        if ( nodesArrayYCWVar) {
                nodesArrayYCWVar = nodesArrayYCWVar - nodeMinYCWVar;
                free ( nodesArrayYCWVar );
        }
}

void MCMF_YCW::set_arc( long tail_node_id, long head_node_id,
                        long low_bound, long up_bound, // up_bound is basically capacity;
                        priceType cost) {
        if ( up_bound < 0 ) {
                up_bound = 0x7fffffff;
        }

        // no of arcs incident to node i is placed in _arc_first[i+1]
        arcFirstYCWVar[tail_node_id + 1] ++;
        arcFirstYCWVar[head_node_id + 1] ++;
        iNode = nodesArrayYCWVar + tail_node_id;
        jNode = nodesArrayYCWVar + head_node_id;

        // store information about the arc
        arcTailYCWVar[posCurrent]   = tail_node_id;
        arcTailYCWVar[posCurrent+1] = head_node_id;
        arcCurrent->set_head( jNode );
        arcCurrent->set_rez_capacity( up_bound - low_bound );
        capArrayYCWVar[posCurrent] = up_bound;
        arcCurrent->set_cost( cost );

        arcCurrent->set_sister( arcCurrent + 1 );
        ( arcCurrent + 1 )->set_head( nodesArrayYCWVar + tail_node_id );
        ( arcCurrent + 1 )->set_rez_capacity( 0 );
        capArrayYCWVar[posCurrent+1] = 0;
        ( arcCurrent + 1 )->set_cost( -cost );
        ( arcCurrent + 1 )->set_sister( arcCurrent );

        iNode->dec_excess( low_bound );
        jNode->inc_excess( low_bound );

        // searching for minimum and maximum node
        if ( head_node_id < nodeMinYCWVar ) nodeMinYCWVar = head_node_id;
        if ( tail_node_id < nodeMinYCWVar ) nodeMinYCWVar = tail_node_id;
        if ( head_node_id > nodeMaxYCWVar ) nodeMaxYCWVar = head_node_id;
        if ( tail_node_id > nodeMaxYCWVar ) nodeMaxYCWVar = tail_node_id;

        if ( cost < 0 ) cost = -cost;
        if ( cost > maxCost && up_bound > 0 ) maxCost = cost;

        // prepare for next arc to be added;
        arcCurrent += 2;
        posCurrent += 2;
}

void MCMF_YCW::set_supply_demand_of_node( long id, long excess) {
        // set supply and demand of nodes; not used for transhipment nodes;
        (nodesArrayYCWVar + id)->set_excess( excess);
        if ( excess > 0) totalPYCWVar += excess;
        if ( excess < 0) totalNYCWVar -= excess;
}

void MCMF_YCW::pre_processing() {
        // called after the arcs were just added and before run_cs2();
        // ordering arcs - linear time algorithm;
        long i;
        long last, arc_num, arc_new_num;;
        long tail_node_id;
        NODE *head_p;
        ARC *arc_new, *arc_tmp;
        long up_bound;
        priceType cost; // arc cost;
        excessType cap_out; // sum of outgoing capacities
        excessType cap_in; // sum of incoming capacities

        // first arc from the first node
        ( nodesArrayYCWVar + nodeMinYCWVar )->set_first( arcsArrayYCWVar );

        for ( i = nodeMinYCWVar + 1; i <= nodeMaxYCWVar + 1; i ++ ) {
                arcFirstYCWVar[i] += arcFirstYCWVar[i-1];
                ( nodesArrayYCWVar + i )->set_first( arcsArrayYCWVar + arcFirstYCWVar[i] );
        }

        // scanning all the nodes except the last
        for ( i = nodeMinYCWVar; i < nodeMaxYCWVar; i ++ ) {

                last = ( ( nodesArrayYCWVar + i + 1 )->first() ) - arcsArrayYCWVar;

                for ( arc_num = arcFirstYCWVar[i]; arc_num < last; arc_num ++ ) {
                        tail_node_id = arcTailYCWVar[arc_num];

                        while ( tail_node_id != i ) {

                                arc_new_num = arcFirstYCWVar[tail_node_id];
                                arcCurrent = arcsArrayYCWVar + arc_num;
                                arc_new = arcsArrayYCWVar + arc_new_num;

                                head_p = arc_new->head();
                                arc_new->set_head( arcCurrent->head() );
                                arcCurrent->set_head( head_p );

                                up_bound          = capArrayYCWVar[arc_new_num];
                                capArrayYCWVar[arc_new_num] = capArrayYCWVar[arc_num];
                                capArrayYCWVar[arc_num]     = up_bound;

                                up_bound = arc_new->rez_capacity();
                                arc_new->set_rez_capacity( arcCurrent->rez_capacity() );
                                arcCurrent->set_rez_capacity( up_bound) ;

                                cost = arc_new->cost();
                                arc_new->set_cost( arcCurrent->cost() );
                                arcCurrent->set_cost( cost );

                                if ( arc_new != arcCurrent->sister() ) {
                                        arc_tmp = arc_new->sister();
                                        arc_new->set_sister( arcCurrent->sister() );
                                        arcCurrent->set_sister( arc_tmp );

                                        arcCurrent->sister()->set_sister( arcCurrent );
                                        arc_new->sister()->set_sister( arc_new );
                                }

                                arcTailYCWVar[arc_num] = arcTailYCWVar[arc_new_num];
                                arcTailYCWVar[arc_new_num] = tail_node_id;

                                // we increase arc_first[tail_node_id]
                                arcFirstYCWVar[tail_node_id] ++ ;

                                tail_node_id = arcTailYCWVar[arc_num];
                        }
                }
                // all arcs outgoing from  i  are in place
        }

        for ( NODE *ndp = nodesArrayYCWVar + nodeMinYCWVar; ndp <= nodesArrayYCWVar + nodeMaxYCWVar; ndp ++ ) {
                cap_in  =   ( ndp->excess() );
                cap_out = - ( ndp->excess() );
                for ( arcCurrent = ndp->first(); arcCurrent != (ndp+1)->first();
                        arcCurrent ++ ) {
                        arc_num = arcCurrent - arcsArrayYCWVar;
                        if ( capArrayYCWVar[arc_num] > 0 ) cap_out += capArrayYCWVar[arc_num];
                        if ( capArrayYCWVar[arc_num] == 0 )
                                cap_in += capArrayYCWVar[ arcCurrent->sister() - arcsArrayYCWVar ];
                }
        }

        // adjustments due to nodes' ids being between _node_min - _node_max;
        nodeNumYCWVar = nodeMaxYCWVar - nodeMinYCWVar + 1;
        nodesArrayYCWVar = nodesArrayYCWVar + nodeMinYCWVar;

        // () free internal memory, not needed anymore inside CS2;
        free ( arcFirstYCWVar );
        free ( arcTailYCWVar );
}

void MCMF_YCW::cs2_initialize() {
        NODE *i; // current node
        ARC *a; // current arc
        ARC *a_stop;
        BUCKET *b; // current bucket
        long df;

        fScaleYCWVar = (long) 12.0;
        sentinelNode = nodesArrayYCWVar + nodeNumYCWVar;
        sentinelArc  = arcsArrayYCWVar + arcNumYCWVar;

        for ( i = nodesArrayYCWVar; i != sentinelNode; i ++ ) {
                i->set_price( 0);
                i->set_suspended( i->first());

                i->set_q_next( sentinelNode);
        }

        sentinelNode->set_first( sentinelArc);
        sentinelNode->set_suspended( sentinelArc);


        for ( i = nodesArrayYCWVar; i != sentinelNode; i ++ ) {
                for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++ ) {
                        if ( a->cost() < 0) {
                                if ( ( df = a->rez_capacity()) > 0) {
                                        increase_flow( i, a->head(), a, df);
                                }
                        }
                }
        }

        _dn = nodeNumYCWVar + 1;
        if ( noZeroCycles == true) { // NO_ZERO_CYCLES
                _dn = 2 * _dn;
        }

        for ( a = arcsArrayYCWVar; a != sentinelArc; a ++ ) {
                a->multiply_cost( _dn);
        }

        if ( noZeroCycles == true) { // NO_ZERO_CYCLES
                for ( a = arcsArrayYCWVar; a != sentinelArc; a ++ ) {
                        if ((a->cost() == 0) && (a->sister()->cost() == 0)) {
                                a->set_cost( 1);
                                a->sister()->set_cost( -1);
                        }
                }
        }
        _mmc = maxCost * _dn;

        _linf = (long) (_dn * ceil(fScaleYCWVar) + 2);

        bucketsArrayYCWVar = (BUCKET*) calloc ( _linf, sizeof(BUCKET));


        lBucket = bucketsArrayYCWVar + _linf;

        dTwoNode = new NODE; // used as reference;

        for ( b = bucketsArrayYCWVar; b != lBucket; b ++ ) {
                reset_bucket( b);
        }

        epsilonYCWVar = _mmc;
        if ( epsilonYCWVar < 1) {
                epsilonYCWVar = 1;
        }

        priceMinYCWVar = -(0x7fffffffffffffffLL);

        cutOffFactor = 1.5 * pow( (double)nodeNumYCWVar, 0.44);

        cutOffFactor = cutOffFactor > 12  ?  cutOffFactor : 12;

        nRefYCWVar = 0;

        flagPrice = 0;

        _dummy_node = &dOneNode;

        excqFirst = NULL;

}

void MCMF_YCW::up_node_scan( NODE *i) {
        NODE *j; // opposite node
        ARC *a; // (i, j)
        ARC *a_stop; // first arc from the next node
        ARC *ra; // (j, i)
        BUCKET *b_old; // old bucket contained j
        BUCKET *b_new; // new bucket for j
        long i_rank;
        long j_rank; // ranks of nodes
        long j_new_rank;
        priceType rc; // reduced cost of (j, i)
        priceType dr; // rank difference

        nScanYCWVar ++;

        i_rank = i->rank();

        // scanning arcs;
        for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++ ) {

                ra = a->sister();

                if ( ra->rez_capacity() > 0 ) {
                        j = a->head();
                        j_rank = j->rank();

                        if ( j_rank > i_rank ) {
                                if ( ( rc = j->price() + ra->cost() - i->price() ) < 0 ) {
                                        j_new_rank = i_rank;
                                } else {
                                        dr = rc / epsilonYCWVar;
                                        j_new_rank = ( dr < _linf ) ? i_rank + (long)dr + 1 : _linf;
                                }

                                if ( j_rank > j_new_rank ) {
                                        j->set_rank( j_new_rank);
                                        j->set_current( ra);

                                        if ( j_rank < _linf ) {
                                                b_old = bucketsArrayYCWVar + j_rank;
                                                // REMOVE_FROM_BUCKET( j, b_old );
                                                if ( j == ( b_old -> p_first() ) )
                                                        b_old ->set_p_first( j -> b_next() );
                                                else
                                                {
                                                        ( j -> b_prev() )->set_b_next( j -> b_next() );
                                                        ( j -> b_next() )->set_b_prev( j -> b_prev() );
                                                }
                                        }

                                        b_new = bucketsArrayYCWVar + j_new_rank;
                                        insert_to_bucket( j, b_new );
                                }
                        }
                }
        }

        i->dec_price( i_rank * epsilonYCWVar);
        i->set_rank( -1);
}

void MCMF_YCW::price_update() {
        register NODE *i;
        excessType remain;
        BUCKET *b; 
        priceType dp; 

        nUpdateYCWVar ++;

        for ( i = nodesArrayYCWVar; i != sentinelNode; i ++ ) {
                if ( i->excess() < 0 ) {
                        insert_to_bucket( i, bucketsArrayYCWVar );
                        i->set_rank( 0);
                } else {
                        i->set_rank( _linf);
                }
        }

        remain = totalExcess;
        if ( remain < 0.5 ) return;

        for ( b = bucketsArrayYCWVar; b != lBucket; b ++ ) {

                while ( nonempty_bucket( b) ) {
                        i=(b -> p_first() );
                        b ->set_p_first( i -> b_next() );
                        up_node_scan( i );

                        if ( i ->excess() > 0 ) {
                                remain -= ( i->excess());
                                if ( remain <= 0 ) break;
                        }
                }
                if ( remain <= 0 ) break;
        }

        if ( remain > 0.5 ) flagUpdt = 1;

        dp = ( b - bucketsArrayYCWVar ) * epsilonYCWVar;

        for ( i = nodesArrayYCWVar; i != sentinelNode; i ++ ) {

                if ( i->rank() >= 0 ) {
                        if ( i->rank() < _linf ) {
                                if ( i == ( ( bucketsArrayYCWVar + i->rank()) -> p_first() ) )
                                        ( bucketsArrayYCWVar + i->rank()) ->set_p_first( i -> b_next() );
                                else
                                {
                                        ( i -> b_prev() )->set_b_next( i -> b_next() );
                                        ( i -> b_next() )->set_b_prev( i -> b_prev() );
                                }
                        }
                        if ( i->price() > priceMinYCWVar ) {
                                i->dec_price( dp);
                        }
                }
        }
}

int MCMF_YCW::relabel( NODE *i) {
        register ARC *a; 
        register ARC *a_stop; 
        register ARC *a_max; 
        register priceType p_max; 
        register priceType i_price; 
        register priceType dp; 

        p_max = priceMinYCWVar;
        i_price = i->price();

        a_max = NULL;

        for ( a = i->current() + 1, a_stop = (i + 1)->suspended(); a != a_stop; a ++ ) {

                if ( (a->rez_capacity() > 0) && ( (dp = (a->head()->price() - a->cost())) > p_max ) ) {
                        if ( i_price < dp ) {
                                i->set_current( a);
                                return ( 1);
                        }
                        p_max = dp;
                        a_max = a;
                }
        }

        for ( a = i->first(), a_stop = i->current() + 1; a != a_stop; a ++ ) {
                if ( (a->rez_capacity() > 0) && ( (dp = (a->head()->price() - a->cost())) > p_max ) ) {
                        if ( i_price < dp ) {
                                i->set_current( a);
                                return ( 1);
                        }
                        p_max = dp;
                        a_max = a;
                }
        }

        if ( p_max != priceMinYCWVar ) {
                i->set_price( p_max - epsilonYCWVar);
                i->set_current( a_max);
        } else { 
                if ( i->suspended() == i->first() ) {
                        if ( i->excess() == 0 ) {
                                i->set_price( priceMinYCWVar);
                        } else {
                                if ( nRefYCWVar == 1 ) {
                                        return -1;
                                } else {
                                        return -1;
                                }
                        }
                } else { 
                        flagPrice = 1;
                }
        }

        nRelabelYCWVar ++;
        nRelYCWVar ++;
        return ( 0);
}

void MCMF_YCW::discharge( NODE *i) {
        register ARC *a;// an arc from i
        register NODE *j; // head of a
        register long df; // amoumt of flow to be pushed through a
        excessType j_exc; // former excess of j

        nDischargeYCWVar ++;

        a = i->current();
        j = a->head();

        if ( !(a->rez_capacity() > 0 && (i->price() + a->cost() < j->price())) ) {
                relabel( i );
                a = i->current();
                j = a->head();
        }

        while ( 1 ) {

                j_exc = j->excess();
                if ( j_exc >= 0 ) {

                        df =  i->excess() < a->rez_capacity()?i->excess():a->rez_capacity();
                        if ( j_exc == 0) nSrcYCWVar++;
                        increase_flow( i, j, a, df ); // INCREASE_FLOW
                        nPushYCWVar ++;

                        if ( out_of_excess_q( j ) ) {
                                insert_to_excess_q( j );
                        }
                } else { // j_exc < 0;

                        df = i->excess() < a->rez_capacity()?i->excess():a->rez_capacity();
                        increase_flow( i, j, a, df ); // INCREASE_FLOW
                        nPushYCWVar ++;

                        if ( j->excess() >= 0 ) {
                                if ( j->excess() > 0 ) {
                                        nSrcYCWVar ++;
                                        relabel( j );
                                        insert_to_excess_q( j );
                                }
                                totalExcess += j_exc;
                        } else {
                                totalExcess -= df;
                        }
                }

                if ( i->excess() <= 0) nSrcYCWVar --;
                if ( i->excess() <= 0 || flagPrice ) break;

                relabel( i );

                a = i->current();
                j = a->head();
        }

        i->set_current( a);
}

int MCMF_YCW::price_in() {
        NODE *i; // current node
        NODE *j;
        ARC *a; // current arc from i
        ARC *a_stop; // first arc from the next node
        ARC *b; // arc to be exchanged with suspended
        ARC *ra; // opposite to a
        ARC *rb; // opposite to b
        priceType rc; // reduced cost
        int n_in_bad; // number of priced_in arcs with negative reduced cost
        int bad_found; // if 1 we are at the second scan if 0 we are at the first scan
        excessType i_exc; // excess of i
        excessType df; // an amount to increase flow


        bad_found = 0;
        n_in_bad = 0;

restart:

        for ( i = nodesArrayYCWVar; i != sentinelNode; i ++ ) {

                for ( a = i->first() - 1, a_stop = i->suspended() - 1; a != a_stop; a -- ) {

                        rc = i->price() + a->cost() - a->head()->price();
                        if ( ( rc < 0) && ( a->rez_capacity() > 0) ) { // bad case;
                                if ( bad_found == 0 ) {
                                        bad_found = 1;
                                        update_cut_off();
                                        goto restart;
                                }
                                df = a->rez_capacity();
                                increase_flow( i, a->head(), a, df );

                                ra = a->sister();
                                j  = a->head();

                                i->dec_first();
                                b = i->first();
                                exchange( a, b );

                                if ( a < j->first() ) {
                                        j->dec_first();
                                        rb = j->first();
                                        exchange( ra, rb );
                                }

                                n_in_bad ++;
                        } else {
                                if ( ( rc < cutOn ) && ( rc > -cutOn ) ) {
                                        i->dec_first();
                                        b = i->first();
                                        exchange( a, b );
                                }
                        }
                }
        }


        if ( n_in_bad != 0 ) {

                nBadPricein ++;

                // recalculating excess queue;
                totalExcess = 0;
                nSrcYCWVar = 0;
                reset_excess_q();

                for ( i = nodesArrayYCWVar; i != sentinelNode; i ++ ) {
                        i->set_current( i->first());
                        i_exc = i->excess();
                        if ( i_exc > 0 ) { // i is a source;
                                totalExcess += i_exc;
                                nSrcYCWVar ++;
                                insert_to_excess_q( i );
                        }
                }

                insert_to_excess_q( _dummy_node );
        }

        if ( timeForPriceIn == 4)
                timeForPriceIn = 6;
        if ( timeForPriceIn == 2)
                timeForPriceIn = 4;

        return ( n_in_bad);
}

int MCMF_YCW::refine() {
        NODE *i; // current node
        excessType i_exc; // excess of i
        long np, nr, ns; // variables for additional print
        int pr_in_int; // current number of updates between price_in

        np = nPushYCWVar;
        nr = nRelabelYCWVar;
        ns = nScanYCWVar;

        nRefineYCWVar ++;
        nRefYCWVar ++;
        nRelYCWVar = 0;
        pr_in_int = 0;

        // initialize;
        totalExcess = 0;
        nSrcYCWVar = 0;
        reset_excess_q();

        timeForPriceIn = 2;

        for ( i = nodesArrayYCWVar; i != sentinelNode; i ++ ) {
                i->set_current( i->first());
                i_exc = i->excess();
                if ( i_exc > 0 ) { // i  is a source
                        totalExcess += i_exc;
                        nSrcYCWVar++;
                        insert_to_excess_q( i );
                }
        }

        if ( totalExcess <= 0 ) return -2;;

        // (2) main loop

        while ( 1 ) {

                if ( empty_excess_q() ) {
                        if ( nRefYCWVar > 1 ) {
                                pr_in_int = 0;
                                price_in();
                        }

                        if ( empty_excess_q() ) break;
                }

                i = excqFirst;
                excqFirst = i -> q_next();
                i ->set_q_next( sentinelNode );

                if ( i->excess() > 0 ) {
                        discharge( i );

                        if ( time_for_update() || flagPrice ) {
                                if ( i->excess() > 0 ) {
                                        insert_to_excess_q( i );
                                }

                                if ( flagPrice && ( nRefYCWVar > 1 ) ) {
                                        pr_in_int = 0;
                                        price_in();
                                        flagPrice = 0;
                                }

                                price_update();

                                while ( flagUpdt ) {
                                        if ( nRefYCWVar == 1 ) {
                                                return -1;
                                        } else {
                                                flagUpdt = 0;
                                                update_cut_off();
                                                nBadRelabelYCWVar ++;
                                                pr_in_int = 0;
                                                price_in();
                                                price_update();
                                        }
                                }
                                nRelYCWVar = 0;

                                if ( nRefYCWVar > 1 && (pr_in_int ++ > timeForPriceIn) ) {
                                        pr_in_int = 0;
                                        price_in();
                                }
                        }
                }
        }

        return -2;
}

int MCMF_YCW::price_refine() {
        NODE *i; // current node
        NODE *j; // opposite node
        NODE *ir; // nodes for passing over the negative cycle
        NODE *is;
        ARC *a; // arc (i,j)
        ARC *a_stop; // first arc from the next node
        ARC *ar;
        long bmax;            // number of farest nonempty bucket
        long i_rank;          // rank of node i
        long j_rank;         // rank of node j
        long j_new_rank;      // new rank of node j
        BUCKET *b;              // current bucket
        BUCKET *b_old;          // old and new buckets of current node
        BUCKET *b_new;
        priceType rc = 0; // reduced cost of a
        priceType dr; // ranks difference
        priceType dp;
        int cc;
        // return code: 1 - flow is epsilon optimal
        // 0 - refine is needed
        long df; // cycle capacity
        int nnc; // number of negative cycles cancelled during one iteration
        int snc; // total number of negative cycle cancelled

        nPrefineYCWVar ++;

        cc = 1;
        snc = 0;

        sncMaxYCWVar = 0;


        // (1) main loop
        // while negative cycle is found or eps-optimal solution is constructed
        while ( 1 ) {

                nnc = 0;
                for ( i = nodesArrayYCWVar; i != sentinelNode; i ++ ) {
                        i->set_rank( 0);
                        i->set_inp( 0);
                        i->set_current( i->first());
                }
                reset_stackq();

                for ( i = nodesArrayYCWVar; i != sentinelNode; i ++ ) {
                        if ( i->inp() == 2 ) continue;

                        i->set_b_next( NULL);

                        // deapth first search
                        while ( 1 ) {
                                i->set_inp(1);

                                // scanning arcs from node i starting from current
                                for ( a = i->current(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                                        if ( (a->rez_capacity() > 0) ) {
                                                j = a->head();
                                                if ( i->price() + a->cost() - j->price() < 0 ) {
                                                        if ( j->inp() == 0 ) { // fresh node  - step forward
                                                                i->set_current( a);
                                                                j->set_b_next( i);
                                                                i = j;
                                                                a = j->current();
                                                                a_stop = (j+1)->suspended();
                                                                break;
                                                        }

                                                        if ( j->inp() == 1 ) { // cycle detected
                                                                cc = 0;
                                                                nnc ++;
                                                                i->set_current( a);
                                                                is = ir = i;
                                                                df = 0x7fffffff;

                                                                while ( 1 ) {
                                                                        ar = ir->current();
                                                                        if ( ar->rez_capacity() <= df ) {
                                                                                df = ar->rez_capacity();
                                                                                is = ir;
                                                                        }
                                                                        if ( ir == j ) break;
                                                                        ir = ir->b_next();
                                                                }

                                                                ir = i;

                                                                while ( 1 ) {
                                                                        ar = ir->current();
                                                                        increase_flow( ir, ar->head(), ar, df);
                                                                        if ( ir == j ) break;
                                                                        ir = ir->b_next();
                                                                }

                                                                if ( is != i ) {
                                                                        for ( ir = i; ir != is; ir = ir->b_next() ) {
                                                                                ir->set_inp( 0);
                                                                        }
                                                                        i = is;
                                                                        a = is->current() + 1;
                                                                        a_stop = (is+1)->suspended();
                                                                        break;
                                                                }
                                                        }
                                                }
                                                // if j-color is BLACK - continue search from i
                                        }
                                } // all arcs from i are scanned

                                if ( a == a_stop ) {
                                        // step back
                                        i->set_inp( 2);
                                        nPrscanOneYCWVar ++;
                                        j = i->b_next();
                                        stackq_push( i );
                                        if ( j == NULL ) break;
                                        i = j;
                                        i->inc_current();
                                }

                        } // end of deapth first search
                } // all nodes are scanned


                // () no negative cycle
                // computing longest paths with eps-precision

                snc += nnc;
                if ( snc < sncMaxYCWVar ) cc = 1;
                if ( cc == 0 ) break;
                bmax = 0;

                while ( nonempty_stackq() ) {

                        nPrscanTwoYCWVar ++;
                        // STACKQ_POP( i );
                        i = excqFirst;
                        excqFirst = i -> q_next();
                        i ->set_q_next( sentinelNode );
                        i_rank = i->rank();
                        for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {

                                if (a->rez_capacity() > 0) {
                                        j  = a->head();
                                        rc = i->price() + a->cost() - j->price();

                                        if ( rc < 0 ) { // admissible arc;
                                                dr = (priceType) (( - rc - 0.5 ) / epsilonYCWVar);
                                                if (( j_rank = dr + i_rank ) < _linf ) {
                                                        if ( j_rank > j->rank() )
                                                                j->set_rank( j_rank);
                                                }
                                        }
                                }
                        } // all arcs from i are scanned

                        if ( i_rank > 0 ) {
                                if ( i_rank > bmax ) bmax = i_rank;
                                b = bucketsArrayYCWVar + i_rank;
                                insert_to_bucket( i, b );
                        }
                } // end of while-cycle: all nodes are scanned - longest distancess are computed;


                if ( bmax == 0 ) { // preflow is eps-optimal;
                        break;
                }


                for ( b = bucketsArrayYCWVar + bmax; b != bucketsArrayYCWVar; b -- ) {
                        i_rank = b - bucketsArrayYCWVar;
                        dp = i_rank * epsilonYCWVar;

                        while ( nonempty_bucket( b) ) {
                                // GET_FROM_BUCKET( i, b );
                                i=(b -> p_first() );
                                b ->set_p_first( i -> b_next() );
                                nPrscanYCWVar ++;

                                for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                                        if (a->rez_capacity() > 0) {
                                                j = a->head();
                                                j_rank = j->rank();
                                                if ( j_rank < i_rank ) {
                                                        rc = i->price() + a->cost() - j->price();
                                                        if ( rc < 0 ) {
                                                                j_new_rank = i_rank;
                                                        } else {
                                                                dr = rc / epsilonYCWVar;
                                                                j_new_rank = ( dr < _linf ) ? i_rank - ( (long)dr + 1 ) : 0;
                                                        }
                                                        if ( j_rank < j_new_rank ) {
                                                                if ( cc == 1 ) {
                                                                        j->set_rank( j_new_rank);
                                                                        if ( j_rank > 0 ) {
                                                                                b_old = bucketsArrayYCWVar + j_rank;
                                                                                // REMOVE_FROM_BUCKET( j, b_old );
                                                                                if ( j == ( b_old -> p_first() ) )
                                                                                        b_old ->set_p_first( j -> b_next() );
                                                                                else
                                                                                {
                                                                                        ( j -> b_prev() )->set_b_next( j -> b_next() );
                                                                                        ( j -> b_next() )->set_b_prev( j -> b_prev() );
                                                                                }
                                                                        }
                                                                        b_new = bucketsArrayYCWVar + j_new_rank;
                                                                        insert_to_bucket( j, b_new );
                                                                } else {
                                                                        df = a->rez_capacity();
                                                                        increase_flow( i, j, a, df );
                                                                }
                                                        }
                                                }
                                        } // end if opened arc
                                } // all arcs are scanned

                                i->dec_price( dp);

                        } // end of while-cycle: the bucket is scanned
                } // end of for-cycle: all buckets are scanned

                if ( cc == 0 ) break;

        } // end of main loop



        // (2) finish
        // if refine needed - saturate non-epsilon-optimal arcs;

        if ( cc == 0 ) {
                for ( i = nodesArrayYCWVar; i != sentinelNode; i ++) {
                        for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                                if ( i->price() + a->cost() - a->head()->price() < - epsilonYCWVar ) {
                                        if ( ( df = a->rez_capacity() ) > 0 ) {
                                                increase_flow( i, a->head(), a, df );
                                        }
                                }
                        }
                }
        }

        return ( cc );
}

void MCMF_YCW::compute_prices() {
        NODE *i; // current node
        NODE *j; // opposite node
        ARC *a; // arc (i,j)
        ARC *a_stop; // first arc from the next node
        long bmax; // number of farest nonempty bucket
        long i_rank; // rank of node i
        long j_rank; // rank of node j
        long j_new_rank; // new rank of node j
        BUCKET *b; // current bucket
        BUCKET *b_old; // old and new buckets of current node
        BUCKET *b_new;
        priceType rc; // reduced cost of a
        priceType dr; // ranks difference
        priceType dp;
        int cc; // return code: 1 - flow is epsilon optimal 0 - refine is needed

        nPrefineYCWVar ++;
        cc = 1;

        // (1) main loop
        // while negative cycle is found or eps-optimal solution is constructed
        while ( 1 ) {

                for ( i = nodesArrayYCWVar; i != sentinelNode; i ++) {
                        i->set_rank( 0);
                        i->set_inp( 0);
                        i->set_current( i->first());
                }
                reset_stackq();

                for ( i = nodesArrayYCWVar; i != sentinelNode; i ++ ) {
                        if ( i->inp() == 2 ) continue;

                        i->set_b_next( NULL);
                        // depth first search
                        while ( 1 ) {
                                i->set_inp( 1);

                                // scanning arcs from node i
                                for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                                        if (a->rez_capacity() > 0) {
                                                j = a->head();
                                                if ( i->price() + a->cost() - j->price() < 0 ) {
                                                        if ( j->inp() == 0 ) { // fresh node  - step forward
                                                                i->set_current( a);
                                                                j->set_b_next( i);
                                                                i = j;
                                                                a = j->current();
                                                                a_stop = (j+1)->suspended();
                                                                break;
                                                        }

                                                        if ( j->inp() == 1 ) { // cycle detected; should not happen
                                                                cc = 0;
                                                        }
                                                }
                                                // if j-color is BLACK - continue search from i
                                        }
                                } // all arcs from i are scanned

                                if ( a == a_stop ) {
                                        // step back
                                        i->set_inp( 2);
                                        nPrscanOneYCWVar ++;
                                        j = i->b_next();
                                        stackq_push( i );
                                        if ( j == NULL ) break;
                                        i = j;
                                        i->inc_current();
                                }

                        } // end of deapth first search
                } // all nodes are scanned


                // no negative cycle
                // computing longest paths

                if ( cc == 0 ) break;
                bmax = 0;

                while ( nonempty_stackq() ) {
                        nPrscanTwoYCWVar ++;
                        i = excqFirst;
                        excqFirst = i -> q_next();
                        i ->set_q_next( sentinelNode );
                        i_rank = i->rank();
                        for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                                if (a->rez_capacity() > 0) {
                                        j  = a->head();
                                        rc = i->price() + a->cost() - j->price();


                                        if ( rc < 0 ) {// admissible arc
                                                dr = - rc;
                                                if (( j_rank = dr + i_rank ) < _linf ) {
                                                        if ( j_rank > j->rank() )
                                                                j->set_rank( j_rank);
                                                }
                                        }
                                }
                        } // all arcs from i are scanned

                        if ( i_rank > 0 ) {
                                if ( i_rank > bmax ) bmax = i_rank;
                                b = bucketsArrayYCWVar + i_rank;
                                insert_to_bucket( i, b );
                        }
                } // end of while-cycle: all nodes are scanned - longest distancess are computed;

                if ( bmax == 0 ) {
                        break;
                }

                for ( b = bucketsArrayYCWVar + bmax; b != bucketsArrayYCWVar; b -- ) {
                        i_rank = b - bucketsArrayYCWVar;
                        dp = i_rank;

                        while ( nonempty_bucket( b) ) {
                                i=(b -> p_first() );
                                b ->set_p_first( i -> b_next() );
                                nPrscanYCWVar ++;

                                for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                                        if (a->rez_capacity() > 0) {
                                                j = a->head();
                                                j_rank = j->rank();
                                                if ( j_rank < i_rank ) {
                                                        rc = i->price() + a->cost() - j->price();

                                                        if ( rc < 0 ) {
                                                                j_new_rank = i_rank;
                                                        } else {
                                                                dr = rc;
                                                                j_new_rank = ( dr < _linf ) ? i_rank - ( (long)dr + 1 ) : 0;
                                                        }
                                                        if ( j_rank < j_new_rank ) {
                                                                if ( cc == 1 ) {
                                                                        j->set_rank( j_new_rank);
                                                                        if ( j_rank > 0 ) {
                                                                                b_old = bucketsArrayYCWVar + j_rank;
                                                                                if ( j == ( b_old -> p_first() ) )
                                                                                        b_old ->set_p_first( j -> b_next() );
                                                                                else
                                                                                {
                                                                                        ( j -> b_prev() )->set_b_next( j -> b_next() );
                                                                                        ( j -> b_next() )->set_b_prev( j -> b_prev() );
                                                                                }
                                                                        }
                                                                        b_new = bucketsArrayYCWVar + j_new_rank;
                                                                        insert_to_bucket( j, b_new );
                                                                }
                                                        }
                                                }
                                        } // end if opened arc
                                } // all arcs are scanned

                                i->dec_price( dp);

                        } // end of while-cycle: the bucket is scanned
                } // end of for-cycle: all buckets are scanned

                if ( cc == 0 ) break;

        } // end of main loop
}

void MCMF_YCW::price_out() {
        NODE *i; // current node
        ARC *a; // current arc from i
        ARC *a_stop; // first arc from the next node
        ARC *b; // arc to be exchanged with suspended
        double n_cut_off; // -cut_off
        double rc; // reduced cost

        n_cut_off = - cutOff;

        for ( i = nodesArrayYCWVar; i != sentinelNode; i ++) {
                for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {

                        rc = i->price() + a->cost() - a->head()->price();
                        if ( ( rc > cutOff && a->sister()->rez_capacity() <= 0 ) ||
                                ( rc < n_cut_off && a->rez_capacity() <= 0 ) ) { // suspend the arc

                                b = i->first();
                                i->inc_first();
                                exchange( a, b );
                        }
                }
        }
}

int MCMF_YCW::update_epsilon() {
        // decrease epsilon after epsilon-optimal flow is constructed;
        if ( epsilonYCWVar <= 1 ) return ( 1 );

        epsilonYCWVar = (priceType) (ceil ( (double) epsilonYCWVar / fScaleYCWVar ));
        cutOff = cutOffFactor * epsilonYCWVar;
        cutOn = cutOff * 0.8;

        return ( 0 );
}




int MCMF_YCW::check_eps_opt() {
        NODE *i;
        ARC *a, *a_stop;

        for ( i = nodesArrayYCWVar; i != sentinelNode; i ++) {
                for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {

                        if ( (a->rez_capacity() > 0) && ((i->price() + a->cost() - a->head()->price()) < - epsilonYCWVar) ) {
                                return ( 0);
                        }
                }
        }
        return(1);
}

void MCMF_YCW::init_solution() {
        ARC *a; // current arc (i,j)
        NODE *i; // tail of a
        NODE *j; // head of a
        long df; // residual capacity

        for ( a = arcsArrayYCWVar; a != sentinelArc; a ++ ) {
                if ( a->rez_capacity() > 0 && a->cost() < 0 ) {
                        df = a->rez_capacity();
                        i  = a->sister()->head();
                        j  = a->head();
                        increase_flow( i, j, a, df );
                }
        }
}

void MCMF_YCW::cs_cost_reinit() {
        if ( costRestartYCWVar == false)
                return;

        NODE *i; // current node
        ARC *a;          // current arc
        ARC *a_stop;
        BUCKET *b; // current bucket
        priceType rc, minc, sum;


        for ( b = bucketsArrayYCWVar; b != lBucket; b ++) {
                reset_bucket( b);
        }

        rc = 0;
        for ( i = nodesArrayYCWVar; i != sentinelNode; i ++) {
                rc = rc < i->price()?rc:i->price();
                i->set_first( i->suspended());
                i->set_current( i->first());
                i->set_q_next( sentinelNode);
        }

        // make prices nonnegative and multiply
        for ( i = nodesArrayYCWVar; i != sentinelNode; i ++) {
                i->set_price( (i->price() - rc) * _dn);
        }

        // multiply arc costs
        for (a = arcsArrayYCWVar; a != sentinelArc; a ++) {
                a->multiply_cost( _dn);
        }

        sum = 0;
        for ( i = nodesArrayYCWVar; i != sentinelNode; i ++) {
                minc = 0;
                for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                        if ( ((a->rez_capacity() > 0) && ((rc = i->price() + a->cost() - a->head()->price()) < 0)) )
                                minc =  epsilonYCWVar > -rc?epsilonYCWVar:-rc;
                }
                sum += minc;
        }

        epsilonYCWVar = ceil(sum / _dn);

        cutOffFactor = 1.5 * pow((double)nodeNumYCWVar, 0.44);

        cutOffFactor =  cutOffFactor > 12 ? cutOffFactor:12;

        nRefYCWVar = 0;

        nRefineYCWVar = nDischargeYCWVar = nPushYCWVar = nRelabelYCWVar = 0;
        nUpdateYCWVar = nScanYCWVar = nPrefineYCWVar = nPrscanYCWVar = nPrscanOneYCWVar =
                                               nBadPricein = nBadRelabelYCWVar = 0;

        flagPrice = 0;

        excqFirst = NULL;
}

void MCMF_YCW::cs2_cost_restart( double *objective_cost) {
        // restart after a cost update;
        if ( costRestartYCWVar == false)
                return;

        int cc; // for storing return code;

        cs_cost_reinit();

        cc = update_epsilon();

        if (cc == 0) {
                do { // scaling loop
                        while ( 1 ) {
                                if ( ! price_refine() )
                                        break;

                                if ( nRefYCWVar >= 1 ) {
                                        if ( price_in() )
                                                break;
                                }
                                if ((cc = update_epsilon ()))
                                        break;
                        }
                        if (cc) break;
                        refine();
                        if ( nRefYCWVar >= 1 ) {
                                price_out();
                        }
                        if ( update_epsilon() )
                                break;
                } while ( cc == 0 );
        }

        finishup( objective_cost );
}


void MCMF_YCW::finishup( double *objective_cost) {
        ARC *a; // current arc
        long na; // corresponding position in capacity array
        double obj_internal = 0; // objective
        priceType cs; // actual arc cost
        long flow; // flow through an arc
        NODE *i;

        // (1) NO_ZERO_CYCLES?
        if ( noZeroCycles == true) {
                for ( a = arcsArrayYCWVar; a != sentinelArc; a ++ ) {
                        if ( a->cost() == 1) {
                                assert( a->sister()->cost() == -1);
                                a->set_cost( 0);
                                a->sister()->set_cost( 0);
                        }
                }
        }

        // (2)
        for ( a = arcsArrayYCWVar, na = 0; a != sentinelArc ; a ++, na ++ ) {
                cs = a->cost() / _dn;
                if ( capArrayYCWVar[na]  > 0 && (flow = capArrayYCWVar[na] - a->rez_capacity()) != 0 )
                        obj_internal += (double) cs * (double) flow;
                a->set_cost( cs);
        }

        for ( i = nodesArrayYCWVar; i != sentinelNode; i ++) {
                i->set_price( (i->price() / _dn));
        }

        // (3) COMP_DUALS?
        if ( compDualsYCWVar == true) {
                compute_prices();
        }

        *objective_cost = obj_internal;
}

void MCMF_YCW::cs2( double *objective_cost) {
        // the main calling function;
        int cc = 0; // for storing return code;


        // (1) update epsilon first;
        update_epsilon();


        // (2) scaling loop;
        do {
//              refine();

                int ref = refine();
                if(ref==-1) {
                        *objective_cost = -1;
                        return ;
                }

                if ( nRefYCWVar >= 1 )
                        price_out();

                if ( update_epsilon() )
                        break;

                while (1) {
                        if ( ! price_refine() )
                                break;

                        if ( nRefYCWVar >= 1 ) {
                                if ( price_in() ) break;
                                if ( (cc = update_epsilon()) ) break;
                        }
                }
        } while ( cc == 0 );


        // (3) finishup;
        finishup( objective_cost );
}


int MCMF_YCW::greenTea() {


        double objective_cost;


        // (4) ordering, etc.;

        pre_processing();

        // () CHECK_SOLUTION?
        if ( checkSolution == true) {
                nodeBalanceYCWVar = (long long int *) calloc (nodeNumYCWVar+1, sizeof(long long int));
                for ( NODE *i = nodesArrayYCWVar; i < nodesArrayYCWVar + nodeNumYCWVar; i ++ ) {
                        nodeBalanceYCWVar[i - nodesArrayYCWVar] = i->excess();
                }
        }


        // (5) initializations;
        arcNumYCWVar = 2 * arcNumYCWVar;
        cs2_initialize(); // works already with 2*m;

        // (6) run CS2;
        cs2( &objective_cost );


        if(objective_cost<0){
                if(!isOutputResult){
                        deallocate_arrays();
                }


                return -1;
        }else{
                // cout<<"加总费用之前："<<objective_cost<<endl;
                addServerAndDeployPrice(&objective_cost);

                if(!isOutputResult){
                        // () cleanup;
                        deallocate_arrays();
                }
                return objective_cost;
        }


}

//通过服务器输出流量确定相应的档次价格
pair<int,int> MCMF_YCW::determineDC(int serverOutput){
    map<int, pair<int,int> >::iterator iter;


    for(iter = serverLevel.begin() ; iter != serverLevel.end(); iter++){
        // cout<<"iter->first:"<<iter->first<<endl;
        // cout<<"iter->second价格:"<<iter->second.second<<endl;
        if(serverOutput<=iter->first){
            return iter->second;
        }
    }
    // cout<<"不要走这里啊啊啊啊啊啊啊=================================="<<endl;
    return serverLevel.rbegin()->second;
}



//确定并加上不同档次服务器和网络节点部署的费用
int MCMF_YCW::addServerAndDeployPrice(double *cost){
    // cout<<"不加服务器和部署的费用："<<cost<<endl;

        int total = 0;

        NODE* srcNode = nodesArrayYCWVar + nodeNumYCWVar - 1;
                ARC* srcArc = srcNode->suspended();

                while(srcArc != (srcNode+1)->suspended()) {
                        if((capArrayYCWVar[ srcArc == NULL ? -1 : srcArc - arcsArrayYCWVar ] - srcArc->rez_capacity()) > 0) {
                                NODE* serverNode = srcArc->head();
                                ARC* a = serverNode->suspended();

                                int serverOutput = 0;

                                while(a != (serverNode+1)->suspended()) {

                                        if((capArrayYCWVar[ a == NULL ? -1 : a - arcsArrayYCWVar ] - a->rez_capacity()) > 0) {
                                                serverOutput += (capArrayYCWVar[ a == NULL ? -1 : a - arcsArrayYCWVar ] - a->rez_capacity());

                                        }

                                        ++a;
                                }
                                total += serverOutput;

                int serverFee = determineDC(serverOutput).second;
                // if(serverOutput>250)
                // cout<<"有服务器超过了250++++++++++++++++++++++++++++++++"<<endl;
                // cout<<"服务器："<<serverNode-_nodes<<"总输出："<<serverOutput<<"对应的服务器费用："<<serverFee<<endl;
                *cost += (serverFee + netStateDeployPrice[serverNode-nodesArrayYCWVar]);
//                cout << "^^^" <<serverNode - _nodes << endl;
                serverAllPrice[serverNode-nodesArrayYCWVar] = serverFee + netStateDeployPrice[serverNode-nodesArrayYCWVar];
                        }

                        srcArc++;
                }

}
map<int, int> MCMF_YCW::getServerAllPrice() {
    return this->serverAllPrice;
}


vector<MCMF_YCW::ARC*> edgePath;//存放一条路径上的边
vector<int> path[10000];
vector<int> pathFlow;
int pathNum = 0;
//服务器_档次map 用于生成路径
map<int,int> server_dc_map;
//确定存储每个服务器的档次
void MCMF_YCW::storeServerGrade(){

    // cout<<"storeServerGrade+++"<<endl;

    for(int i=0; i< pathNum ;i++){

        int serverIndex = path[i][1];

        map<int, int>::iterator iter = server_dc_map.find(serverIndex);
        if(iter!=server_dc_map.end()){
            // cout<<"iter->second加之前："<<iter->second;
            iter->second += pathFlow[i];
            // cout<<"iter->second加之后："<<iter->second<<"加了多少："<<outFlow[i]<<endl;
        }else{
            server_dc_map.insert(pair<int,int>(serverIndex,pathFlow[i]));
        }
    }
    map<int, int>::iterator iter;
    int totalFeed = 0;
    for(iter = server_dc_map.begin() ; iter != server_dc_map.end(); iter++){
        totalFeed += iter->second;
        // cout<<iter->first<<"攻击"<<iter->second<<endl;
        pair<int,int> p = determineDC(iter->second);
        iter->second = p.first;
    }
    // cout<<"总供给："<<totalFeed<<endl;
    // cout<<"总需求：k"<<endl;
}

void MCMF_YCW::print_solution(string& result) {
        bool  flag = false;
        for(int k=0;k<10000;k++){
                path[k].clear();
        }

        while(true) {
                NODE* u = nodesArrayYCWVar + nodeNumYCWVar - 1;
                long minFlowOfPath = 10000;
                edgePath.clear();
                while(!((u-nodesArrayYCWVar)>=netStates&&(u-nodesArrayYCWVar)<=(nodeNumYCWVar-2))) {
                        path[pathNum].push_back(u-nodesArrayYCWVar);
                        ARC* a = u->suspended();
                        while(a != (u+1)->suspended()) {
                                if((capArrayYCWVar[ a == NULL ? -1 : a - arcsArrayYCWVar ] - a->rez_capacity()) > 0) {
                                        minFlowOfPath = minFlowOfPath < (capArrayYCWVar[ a == NULL ? -1 : a - arcsArrayYCWVar ] - a->rez_capacity()) ? minFlowOfPath:(capArrayYCWVar[ a == NULL ? -1 : a - arcsArrayYCWVar ] - a->rez_capacity());

                                        edgePath.push_back(a) ;
                                        u = a->head();
                                        break;
                                }
                                ++a;
                        }
                        if(u == nodesArrayYCWVar + nodeNumYCWVar - 1) {
                                flag = true;
                                break;
                        }
                }
                if(flag) {
                        break;
                }
                path[pathNum].push_back(u-nodesArrayYCWVar);
                pathNum++;
                pathFlow.push_back(minFlowOfPath);

                //在路径上删除最小流量
                for(int k=0; k<edgePath.size(); k++) {
                        edgePath[k]->rezCapacityARCVar += minFlowOfPath;
                }
        }
        //确定当前路径的每个服务器的档次
        storeServerGrade();

        //打印路径

         char c[50];
         sprintf(c, "%d\n\n", pathNum);
         result += c;
         for(int i=0;i<pathNum;i++){
                for(int j=1;j<path[i].size();j++){
                        sprintf(c,"%d ",path[i][j]>=netStates?(path[i][j]-netStates):path[i][j]);
                        result += c;
                 }
                 sprintf(c,"%d %d\n",pathFlow[i],server_dc_map[path[i][1]]);
                 result += c;
         }
}
