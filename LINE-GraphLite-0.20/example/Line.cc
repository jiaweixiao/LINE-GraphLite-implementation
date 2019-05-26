/* 1, 2018E8013261005, Zhang Yunfei */

#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <map>

#include "GraphLite.h"

#define VERTEX_CLASS_NAME(name) LINE##name

#define EPS 1e-6
#define SAMPLE_NUM 16
#define NEG_NUM 10
#define VERTEX_DIM 200
#define ROOT 1

typedef float real;
typedef enum {START, END, NEG} vertex_type;//start,end,negative vertex
typedef double VertexWeit;
struct VertexMsg{
    int type; // message type: 1,2,3
    int64_t source_id; // start vertex of edge
    int64_t target_id;   // end vertex of edge

    /* Type 1 Message: all vertexes send vertex and edge data to root vertex  */
    /* ---------------------------------1------------------------------------ */
    VertexWeit edge_weight; // weight of edge
    /* ---------------------------------------------------------------------- */

    /* Type 2 Message: root vertex send sample edge and vertex data to related
     * vertexes */
    /* ---------------------------------2------------------------------------ */
    vertex_type v_type; //start,end,negative vertex
    int64_t neg_vid[NEG_NUM]; // negative sample vertexes array
    /* ---------------------------------------------------------------------- */

    /* Type 3 Message: vertex send vector of vertex to related vertexes */
    /* ---------------------------------3------------------------------------ */
    vertex_type source_vertex_type; // type of vertex which sends message
    vertex_type target_vertex_type; // type of vertex which receives message
    real vec_v[VERTEX_DIM];
    /* ---------------------------------------------------------------------- */

};
struct VertexVal{
    real *emb_vertex;
    real *emb_context;
};


class VERTEX_CLASS_NAME(InputFormatter): public InputFormatter {
public:
    int64_t getVertexNum() {
        unsigned long long n;
        sscanf(m_ptotal_vertex_line, "%lld", &n);
        m_total_vertex= n;
        return m_total_vertex;
    }
    int64_t getEdgeNum() {
        unsigned long long n;
        sscanf(m_ptotal_edge_line, "%lld", &n);
        m_total_edge= n;
        return m_total_edge;
    }
    int getVertexValueSize() {
        m_n_value_size = sizeof(VertexVal);
        return m_n_value_size;
    }
    int getEdgeValueSize() {
        m_e_value_size = sizeof(int);
        return m_e_value_size;
    }
    int getMessageValueSize() {
        m_m_value_size = sizeof(VertexMsg);
        return m_m_value_size;
    }
    void loadGraph() {
        unsigned long long last_vertex;
        unsigned long long from;
        unsigned long long to;
        VertexWeit weight = 0;
        
        struct VertexVal value = 0;
        // init embedding
        posix_memalign((void **)&(value.emb_vertex), 128, (long long)VERTEX_DIM * sizeof(real));
        posix_memalign((void **)&(value.emb_context), 128, (long long)VERTEX_DIM * sizeof(real));
        if (value.emb_vertex == NULL) { printf("Error: memory allocation failed\n"); exit(1); }
        if (value.emb_context == NULL) { printf("Error: memory allocation failed\n"); exit(1); }
        for (int b = 0; b < VERTEX_DIM; b++) {
            emb_vertex[b] = (rand() / (real)RAND_MAX - 0.5) / VERTEX_DIM;
            emb_context[b] = 0;
        }

        int outdegree = 0;
        
        const char *line= getEdgeLine();

        // Note: modify this if an edge weight is to be read
        //       modify the 'weight' variable

        sscanf(line, "%lld %lld %d", &from, &to, &weight);
        addEdge(from, to, &weight);

        last_vertex = from;
        ++outdegree;
        for (int64_t i = 1; i < m_total_edge; ++i) {
            line= getEdgeLine();

            // Note: modify this if an edge weight is to be read
            //       modify the 'weight' variable

            sscanf(line, "%lld %lld %d", &from, &to, &weight);
            if (last_vertex != from) {
                addVertex(last_vertex, &value, outdegree);
                last_vertex = from;
                outdegree = 1;
            } else {
                ++outdegree;
            }
            addEdge(from, to, &weight);
        }
        addVertex(last_vertex, &value, outdegree);
    }
};

class VERTEX_CLASS_NAME(OutputFormatter): public OutputFormatter {
public:
    void writeResult() {
        int64_t vid;
        int current_degree;
        char s[1024];

        for (ResultIterator r_iter; ! r_iter.done(); r_iter.next() ) {
            r_iter.getIdValue(vid, &current_degree);
            if (current_degree > 0){
                int n = sprintf(s, "%lld\n", (unsigned long long)vid);
                writeNextResLine(s, n);
            }
        }
    }
};

// An aggregator that records a int value which denote the number of vertex that send message
class VERTEX_CLASS_NAME(Aggregator): public Aggregator<int> {
private:
    int k_val; //input parameter K
public:
    void init() {
        m_global = k_val;
        m_local = 0;
    }
    void* getGlobal() {
        return &m_global;
    }
    void setGlobal(const void* p) {
        m_global = * (int *)p;
    }
    void* getLocal() {
        return &m_local;
    }
    int getKValue() {
        return k_val;
    }
    void setKValue(int k) {
        k_val = k;
    }
    void merge(const void* p) {
        m_global += * (int *)p;
    }
    void accumulate(const void* p) {
        m_local += * (int *)p;
    }
};

class VERTEX_CLASS_NAME(): public Vertex <VertexVal, VertexWeit, VertexMsg> {
private:
    int k; //input parameter K

    /* Parameters in ROOT */
    long long num_edges;

    // Parameters for edge sampling
    long long *alias;
    double *prob;

    std::vector<VertexWeit> *edge_weight;
    std::vector<int64_t> *edge_source_id, *edge_target_id;
    int64_t *vertex_hash_table, *neg_table;

public:
    void compute(MessageIterator* pmsgs) {

        int64_t vid = getVertexId();

        if (getSuperstep() == 0) {
            // all vertices send edge information to root vertex
            OutEdgeIterator outEdgeIterator = getOutEdgeIterator()
            VertexMsg msg;
            msg.type = 1;
            msg.source_vid = vid;
            for ( ; ! outEdgeIterator.done(); outEdgeIterator.next() ) {
                msg.target_vid = outEdgeIterator.target();
                msg.edge_weight = outEdgeIterator.getValue();      
                sendMessageTo(ROOT, msg);
            }
        } else if (getSuperstep() == 1){
            // root vertex initialize sample tables in superStep 1
            if(vid == ROOT){
                num_edges = 0;
                edge_source_id = new std::vector<int64_t>();
                edge_target_id = new std::vector<int64_t>();
                edge_weight = new std::vector<VertexWeit>();
                // process message
                for ( ; ! pmsgs->done(); pmsgs->next() ) {
                    VertexMess msg = pmsgs->getValue();
                    if (msg.type == 1){
                        num_edges += 1;
                        edge_source_id->push_back(msg.source_id);
                        edge_target_id->push_back(msg.target_id);
                        edge_weight->push_back(msg.edge_weight);
                    }
                }

                InitAliasTable();
            } else {

            }
        } else {
            //root vertex sample and send result to related vertex;
            if(vid == ROOT){

            }

            if(getSuperstep() != 1 || vid != ROOT){
                // process message
                for ( ; ! pmsgs->done(); pmsgs->next() ) {
                    VertexMess msg = pmsgs->getValue();
                    if (msg.type == 2){

                    } else if (msg.type == 3){

                    }
                }
            }
        }
    }
	/* The alias sampling algorithm, which is used to sample an edge in O(1) time. */
    void InitAliasTable() {
    	alias = (long long *)malloc(num_edges*sizeof(long long));
    	prob = (double *)malloc(num_edges*sizeof(double));
    	if (alias == NULL || prob == NULL) {
    		printf("Error: memory allocation failed!\n");
    		exit(1);
    	}
    	double *norm_prob = (double *)malloc(num_edges*sizeof(double));
    	long long *large_block = (long long *)malloc(num_edges*sizeof(long long));
    	long long *small_block = (long long *)malloc(num_edges*sizeof(long long));
    	if (norm_prob == NULL || large_block == NULL || small_block == NULL) {
    		printf("Error: memory allocatiion failed!\n");
    		exit(1);
    	}

    	double sum = 0;
    	long long cur_small_block, cur_large_block;
    	long long num_small_block = 0, num_large_block = 0;

    	/* need edge_weight[] */
    	for (long long k = 0; k != num_edges; k++) sum += edge_weight[k];
    	for (long long k = 0; k != num_edges; k++) norm_prob[k] = edge_weight[k] * num_edges / sum;

    	for (long long k = num_edges - 1; k >= 0; k--) {
    		if (norm_prob[k] < 1)
    			small_block[num_small_block++] = k;
    		else
    			large_block[num_large_block++] = k;
    	}

    	while (num_small_block && num_large_block) {
    		cur_small_block = small_block[--num_small_block];
    		cur_large_block = large_block[--num_large_block];
    		prob[cur_small_block] = norm_prob[cur_small_block];
    		alias[cur_small_block] = cur_large_block;
    		norm_prob[cur_large_block] = norm_prob[cur_large_block] + norm_prob[cur_small_block] - 1;
    		if (norm_prob[cur_large_block] < 1)
    			small_block[num_small_block++] = cur_large_block;
    		else
    			large_block[num_large_block++] = cur_large_block;
    	}

    	while (num_large_block) prob[large_block[--num_large_block]] = 1;
    	while (num_small_block) prob[small_block[--num_small_block]] = 1;

    	free(norm_prob);
    	free(small_block);
    	free(large_block);
    }
};

class VERTEX_CLASS_NAME(Graph): public Graph {
public:
    VERTEX_CLASS_NAME(Aggregator)* aggregator;

public:
    // argv[0]: PageRankVertex.so
    // argv[1]: <input path>
    // argv[2]: <output path>
    // argv[3]: k
    void init(int argc, char* argv[]) {

        setNumHosts(5);
        setHost(0, "localhost", 1411);
        setHost(1, "localhost", 1421);
        setHost(2, "localhost", 1431);
        setHost(3, "localhost", 1441);
        setHost(4, "localhost", 1451);

        if (argc < 3) {
           printf ("Usage: %s <input path> <output path>\n", argv[0]);
           exit(1);
        }

        m_pin_path = argv[1];
        m_pout_path = argv[2];

        aggregator = new VERTEX_CLASS_NAME(Aggregator)[1];
        regNumAggr(1);
        regAggr(0, &aggregator[0]);
        int k_val = atoi(argv[3]);
        aggregator->setKValue(k_val);
        printf("K value: %d\n", k_val);
    }

    void term() {
        delete[] aggregator;
    }
};

/* STOP: do not change the code below. */
extern "C" Graph* create_graph() {
    Graph* pgraph = new VERTEX_CLASS_NAME(Graph);

    pgraph->m_pin_formatter = new VERTEX_CLASS_NAME(InputFormatter);
    pgraph->m_pout_formatter = new VERTEX_CLASS_NAME(OutputFormatter);
    pgraph->m_pver_base = new VERTEX_CLASS_NAME();

    return pgraph;
}

extern "C" void destroy_graph(Graph* pobject) {
    delete ( VERTEX_CLASS_NAME()* )(pobject->m_pver_base);
    delete ( VERTEX_CLASS_NAME(OutputFormatter)* )(pobject->m_pout_formatter);
    delete ( VERTEX_CLASS_NAME(InputFormatter)* )(pobject->m_pin_formatter);
    delete ( VERTEX_CLASS_NAME(Graph)* )pobject;
}