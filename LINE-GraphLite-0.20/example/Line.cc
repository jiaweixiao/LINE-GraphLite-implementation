/* 1, 2018E8013261005, Zhang Yunfei */

#include <stdio.h>
#include <string.h>
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
struct VertexMsg{
    int type; // message type: 1,2,3
    int64_t start_vid; // start vertex of edge
    int64_t end_vid;   // end vertex of edge

    /* Type 1 Message: all vertexes send vertex and edge data to root vertex  */
    /* ---------------------------------1------------------------------------ */
    int w_edge; // weight of edge
    long degree_weight_sum; // all degree weight of vertex
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
    real emb_vertex[VERTEX_DIM];
    real emb_context[VERTEX_DIM];
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
        int weight = 0;
        
        int value = 0;
        int outdegree = 0;
        
        const char *line= getEdgeLine();

        // Note: modify this if an edge weight is to be read
        //       modify the 'weight' variable

        sscanf(line, "%lld %lld", &from, &to);
        addEdge(from, to, &weight);

        last_vertex = from;
        ++outdegree;
        for (int64_t i = 1; i < m_total_edge; ++i) {
            line= getEdgeLine();

            // Note: modify this if an edge weight is to be read
            //       modify the 'weight' variable

            sscanf(line, "%lld %lld", &from, &to);
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

class VERTEX_CLASS_NAME(): public Vertex <VertexVal, int, VertexMsg> {
private:
    int k; //input parameter K

    /* Parameters in ROOT */
    long long num_edges;
    // Parameters for edge sampling
    long long *alias;
    double *prob;

    int *vertex_hash_table, *neg_table;

public:
    void compute(MessageIterator* pmsgs) {

        int64_t vid = getVertexId();

        if (getSuperstep() == 0) {
            // all vertexes send edge information to root vertex
            OutEdgeIterator outEdgeIterator = getOutEdgeIterator()
			// int64_t out_degree = outEdgeIterator.size();
            int weight_sum = 0;
            VertexMess message;
            message.type = 1;
            message.start_vid = vid;
            for ( ; ! outEdgeIterator.done(); outEdgeIterator.next() ) {
                int64_t targetId = outEdgeIterator.target();
                int weight = outEdgeIterator.getValue();
                weight_sum += weight;
                message.end_vid = targetId;
                message.w_edge = weight;
                message.degree_weight_sum = weight_sum;
                sendMessageTo(ROOT, message);
            }
        } else {
            // root vertex initialize sample tables in superStep 1
            if(getSuperstep() == 1 && vid == ROOT){

            }
            //root vertex sample and send result to related vertex;
            if(vid == ROOT){

            }

            if(getSuperstep() != 1 || vid != ROOT){
                // process message
                for ( ; ! pmsgs->done(); pmsgs->next() ) {

                    VertexMess mess = pmsgs->getValue();
                    if (mess.type == 2){

                    } else if (mess.type == 3){

                    }
                }
            }
        }
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
