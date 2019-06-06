/**
 * @author  SUN Qian, XIAO Jiawei, ZHANG Yunfei
 * @version 0.1
 * 
 * @email sunqian18s@ict.ac.cn, zhangyunfei18s@ict.ac.cn, zhangyunfei18s@ict.ac.cn
 *
 */

#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <map>
#include <gsl/gsl_rng.h>

#include "GraphLite.h"

#define VERTEX_CLASS_NAME(name) LINE##name

#define EPS 1e-6
#define SAMPLE_NUM 16
#define NEG_NUM 10
#define VERTEX_DIM 200
#define ROOT 1

#define MAX_STRING 100
#define SIGMOID_BOUND 6
#define NEG_SAMPLING_POWER 0.75
#define NEG_TABLE_SIZE 1e8

#define SIGMOID_TABLE_SIZE 1000

#define total_samples 100000


int order = 1; //input parameter order
int parallel_num = 1;

typedef float real;
typedef enum {START, END, NEG} vertex_type;//start, end or negative vertex
typedef double VertexWeit;
struct VertexMsg{
    int type; // message type: 0,1,2,3
    int64_t source_id; // source vertex of edge
    int64_t target_id; // target vertex of edge
    real rho_m;

    /* Type 0 Message: all vertices send vertex and edge data to root vertex  */
    /* ---------------------------------1------------------------------------ */
    VertexWeit weight;
    /* ---------------------------------------------------------------------- */

    /* Type 1 Message: all vertices send it's degree to root vertex  */
    /* ---------------------------------1------------------------------------ */
    VertexWeit degree;
    /* ---------------------------------------------------------------------- */

    /* Type 2 Message: root sends sample edge and vertex data to related
     * vertexes */
    /* ---------------------------------2------------------------------------ */
    int order; // order of proximity
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

        sscanf(line, "%lld %lld %lf", &from, &to, &weight);
        addEdge(from, to, &weight);

        last_vertex = from;
        ++outdegree;
        for (int64_t i = 1; i < m_total_edge; ++i) {
            line= getEdgeLine();

            // Note: modify this if an edge weight is to be read
            //       modify the 'weight' variable

            sscanf(line, "%lld %lld %lf", &from, &to, &weight);
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
    /* Parameters in ROOT */
    long long num_edges;
    real init_rho = 0.025, rho;
    long long count = 0, last_count = 0, current_sample_count = 0;
    unsigned long long seed = 1;

    // Parameters for edge sampling
    long long *alias;
    double *prob;

    int num_vertices;
    std::vector<VertexWeit> *edge_weight, *vertex_degree;
    std::vector<int64_t> *edge_source_id, *edge_target_id;
    std::vector<pair<int64_t,VertexWeit>> *vid_map; //vid:degree
    int *neg_table;
    real *sigmoid_table;

public:
    void compute(MessageIterator* pmsgs) {

        int64_t vid = getVertexId();

        if (getSuperstep() == 0) {
            // all vertices send edge information to root vertex
            OutEdgeIterator outEdgeIterator = getOutEdgeIterator()
            VertexMsg msg;
            VertexWeit degree = 0;
            msg.type = 0;
            msg.source_vid = vid;
            for ( ; ! outEdgeIterator.done(); outEdgeIterator.next() ) {
                msg.target_vid = outEdgeIterator.target();
                msg.weight = outEdgeIterator.getValue();
                degree += msg.weight;    
                sendMessageTo(ROOT, msg);
            }
            msg.type = 1;
            msg.degree = degree;
            sendMessageTo(ROOT, msg);
        } else if (getSuperstep() == 1){
            // root receives messages and initiates
            if(vid == ROOT){
                num_edges = 0;
                num_vertices = 0;
                edge_source_id = new std::vector<int64_t>();
                edge_target_id = new std::vector<int64_t>();
                edge_weight = new std::vector<VertexWeit>();
                vid_map = new std::vector<pair<int64_t,VertexWeit>>();
                // process message
                for ( ; ! pmsgs->done(); pmsgs->next() ) {
                    VertexMsg msg = pmsgs->getValue();
                    if (msg.type == 0){
                        num_edges += 1;
                        edge_source_id->push_back(msg.source_id);
                        edge_target_id->push_back(msg.target_id);
                        edge_weight->push_back(msg.weight);
                    }
                    if (msg.type == 1) {
                        num_vertices += 1;
                        vid_map->push_back(pair<int64_t,VertexWeit>(msg.source_id,msg.degree));
                    }
                }
                num_vertices = vid_map->size();
                // let index of vertex_degree is same to its vid
                vertex_degree = new std::vector<VertexWeit>(num_vertices);
                for (int64_t i = 0; i < num_vertices; ++i) {
                    *(vertex_degree->begin()+vid_map->at(i).first) = vid_map->at(i).second;
                }


                InitAliasTable();
                InitNegTable();
                InitSigmoidTable();

                gsl_rng_env_setup();
                gsl_T = gsl_rng_rand48;
                gsl_r = gsl_rng_alloc(gsl_T);
                gsl_rng_set(gsl_r, 314159265);

                // root vertex samples and sends message
                for (int i = 0; i < parallel_num; ++i) {
                    SampleAndSendMsg();
                }

            } else {
                // all vertexes initialize vectors superStep 1
                // TODO
            }
        } else {
            //root vertex sample and send result to related vertex;
            if(vid == ROOT){
                //judge for exit
                if (count < total_samples / parallel_num + 2){
                    if (count - last_count > 10000) {
                        current_sample_count += count - last_count;
                        last_count = count;
                        printf("%cRho: %f  Progress: %.3lf%%", 13, rho, (real)current_sample_count / (real)(total_samples + 1) * 100);
                        fflush(stdout);
                        rho = init_rho * (1 - current_sample_count / (real)(total_samples + 1));
                        if (rho < init_rho * 0.0001) rho = init_rho * 0.0001;
                    }
                    // root vertex samples and sends message
                    for (int i = 0; i < parallel_num; ++i) {
                        SampleAndSendMsg();
                    }
                }
            }

            // process message
            for ( ; ! pmsgs->done(); pmsgs->next() ) {
                VertexMsg msg = pmsgs->getValue();
                if (msg.type == 2){
                    /* this is a source vertex and need to send it vector to target and neg vertexes */
                    VertexVal val = getValue();
                    VertexMsg message;
                    message.type = 3;
                    message.source_id = vid;
                    message.source_vertex_type = START;
                    message.rho_m = msg.rho_m;
                    for (int i = 0; i < VERTEX_DIM; ++i) {
                        message.vec_v[i] = val.emb_vertex[i];
                    }
                    message.target_id = msg.target_id;
                    message.target_vertex_type = END;
                    // send message to end vertex of the edge
                    sendMessageTo(msg.target_id, message);
                    // send message to negative vertexes of the edge
                    for (int j = 0; j < NEG_NUM; ++j) {
                        message.target_id = msg.neg_vid[j];
                        message.target_vertex_type = NEG;
                        sendMessageTo(msg.neg_vid[j], message);
                    }
                } else if (msg.type == 3){
                    real vec_error[VERTEX_DIM]={0};
                    if(msg.target_vertex_type == START){
                        // update start vertex
                        VertexVal val = getValue();
                        for (int i = 0; i < VERTEX_DIM; ++i) {
                            val.emb_vertex[i] += msg.vec_v[i];
                        }

                    } else if(msg.target_vertex_type == END){
                        // update end vertex
                        VertexVal val = getValue();
                        if (order == 1) Update(msg.vec_v, val.emb_vertex, vec_error, 1, msg.rho_m);
                        if (order == 2) Update(msg.vec_v, val.emb_context, vec_error, 1, msg.rho_m);

                        VertexMsg message;
                        message.type = 3;
                        message.source_id = vid;
                        message.source_vertex_type = END;
                        for (int i = 0; i < VERTEX_DIM; ++i) {
                            message.vec_v[i] = vec_error[i];
                        }
                        message.target_id = msg.source_id;
                        message.target_vertex_type = START;
                        // send message to start vertex of the edge
                        sendMessageTo(msg.source_id, message);

                    } else if(msg.target_vertex_type == NEG){
                        // update negative vertex
                        VertexVal val = getValue();
                        if (order == 1) Update(msg.vec_v, val.emb_vertex, vec_error, 0, msg.rho_m);
                        if (order == 2) Update(msg.vec_v, val.emb_context, vec_error, 0, msg.rho_m);
                        VertexMsg message;
                        message.type = 3;
                        message.source_id = vid;
                        message.source_vertex_type = NEG;
                        for (int i = 0; i < VERTEX_DIM; ++i) {
                            message.vec_v[i] = vec_error[i];
                        }
                        message.target_id = msg.source_id;
                        message.target_vertex_type = START;
                        // send message to start vertex of the edge
                        sendMessageTo(msg.source_id, message);
                    }
                }
            }
        }
    }
    /* root vertex samples and sends message to vertex*/
    void SampleAndSendMsg(){
        int64_t u, v;
        curedge = SampleAnEdge(gsl_rng_uniform(gsl_r), gsl_rng_uniform(gsl_r));
        u = edge_source_id->at(curedge);
        v = edge_target_id->at(curedge);
        VertexMsg msg;
        msg.type = 2;
        msg.source_vid = u;
        msg.target_id = v;
        msg.order = order;
        msg.rho_m = rho;
        // NEGATIVE SAMPLING
        for (int d = 0; d < NEG_NUM; d++) {
            msg.neg_vid[d] = neg_table[Rand(seed)];
        }
        sendMessageTo(u, msg);
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
    
    long long SampleAnEdge(double rand_value1, double rand_value2) {
        long long k = (long long)num_edges * rand_value1;
        return rand_value2 < prob[k] ? k : alias[k];
    }

    /* Sample negative vertex samples according to vertex degrees */
    void InitNegTable()
    {
        double sum = 0, cur_sum = 0, por = 0;
        int vid = 0;
        neg_table = (int *)malloc(NEG_TABLE_SIZE * sizeof(int));
        for (int k = 0; k != num_vertices; k++) sum += pow(vertex_degree.at(k), NEG_SAMPLING_POWER);
        for (int k = 0; k != NEG_TABLE_SIZE; k++)
        {
            if ((double)(k + 1) / NEG_TABLE_SIZE > por)
            {
                cur_sum += pow(vertex_degree.at(vid), NEG_SAMPLING_POWER);
                por = cur_sum / sum;
                vid++;
            }
            neg_table[k] = vid - 1;
        }
    }

    /* Fastly compute sigmoid function */
    void InitSigmoidTable(){
        real x;
        sigmoid_table = (real *)malloc((SIGMOID_TABLE_SIZE + 1) * sizeof(real));
        for (int k = 0; k != SIGMOID_TABLE_SIZE; k++)
        {
            x = 2.0 * SIGMOID_BOUND * k / SIGMOID_TABLE_SIZE - SIGMOID_BOUND;
            sigmoid_table[k] = 1 / (1 + exp(-x));
        }
    }
    real FastSigmoid(real x){
        if (x > SIGMOID_BOUND) return 1;
        else if (x < -SIGMOID_BOUND) return 0;
        int k = (x + SIGMOID_BOUND) * SIGMOID_TABLE_SIZE / SIGMOID_BOUND / 2;
        return sigmoid_table[k];
    }

    /* Update embeddings */
    void Update(real *vec_u, real *vec_v, real *vec_error, int label, real rho_m){
        real x = 0, g;
        for (int c = 0; c != VERTEX_DIM; c++) x += vec_u[c] * vec_v[c];
        g = (label - FastSigmoid(x)) * rho_m;
        for (int c = 0; c != VERTEX_DIM; c++) vec_error[c] += g * vec_v[c];
        for (int c = 0; c != VERTEX_DIM; c++) vec_v[c] += g * vec_u[c];
    }

    /* Fastly generate a random integer */
    int Rand(unsigned long long &seed) {
        seed = seed * 25214903917 + 11;
        return (seed >> 16) % NEG_TABLE_SIZE;
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