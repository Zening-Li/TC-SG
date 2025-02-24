#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>

#include <string.h>
#include <math.h>

#include "include/mt19937ar.h"
#include "include/stats.hpp"
#include "include/parallel_hashmap/phmap.h"


using namespace std;
using phmap::flat_hash_map;

// initialization of statslib
stats::rand_engine_t engine(12306);


struct edge{
    int neighbor;
    int8_t sign;
};

static inline bool edge_less_int(const edge &e, int value) {
    return e.neighbor < value;
}

static inline bool edge_comp(const edge &a, const edge &b){
    return a.neighbor < b.neighbor;
}


FILE *open_file(const string &file_path, const char *mode){
    FILE *fp;
    if((fp = fopen(file_path.c_str(), mode)) == NULL){
        cout << "cannot open " << file_path << endl;
        exit(-1);
    }
    return fp;
}


bool check_file_existence(const string &file_path){
    ifstream ifs(file_path);
    return ifs.is_open();
}


int total_nodes_count(const string &edge_path){
    FILE *fp;
    char str[1025];
    char *ptr, *temp;
    int total_num_nodes = 0;

    fp = open_file(edge_path, "r");
    // find the number of nodes
    temp = fgets(str, 1024, fp);
    ptr = strchr(str, '#');
    total_num_nodes = atoi(ptr + 1);
    fclose(fp);
    // cout << "total_num_nodes: " << total_num_nodes << endl;
    return total_num_nodes;
}


// randomly generate 0, 1, 2, ..., size-1, and store the first 'num' values into rndperm
void create_random_permutation(std::vector<int> &rndperm, int size, int num) {
    std::vector<int> ordperm(size);
    // create an ordered permutation
    std::iota(ordperm.begin(), ordperm.end(), 0);
    
    // Fisher-Yates shuffle
    for (int i = 0; i < num; ++i) {
        int rnd = i + (genrand_int32() % (size - i));
        std::swap(ordperm[i], ordperm[rnd]);
        rndperm[i] = ordperm[i];
    }
    ordperm.clear();
    ordperm.shrink_to_fit();
}


void create_node_order(const string &dataset, int &num_sample_nodes, const int &total_num_nodes, int &num_iters, vector<vector<int>> &node_order){
    if(num_sample_nodes == total_num_nodes){ // use all nodes to calculate the number of triangles
        num_iters = 1;
        node_order = vector<vector<int>>(num_iters, vector<int>(total_num_nodes));
        for(int i=0; i<num_sample_nodes; ++i){
            node_order[0][i] = i;
        }
    } else{ // randomly generate the order of nodes --> node_order
        FILE *fp;
        string outfile = "./data/"+dataset+"_node_order_iters_"+to_string(num_iters)+".csv";
        int i, j;
        node_order = vector<vector<int>>(num_iters, vector<int>(total_num_nodes));
        if(check_file_existence(outfile)){
            fp = open_file(outfile, "r");
            char str[1025];
            char *temp;
            for(i=0; i<total_num_nodes; ++i){
                temp = fgets(str, 1024, fp);
                strtok(str, ",");
                for(j=0; j<num_iters; ++j){
                    node_order[j][i] = atoi(strtok(NULL, ","));
                }
            }
            fclose(fp);
        } else{
            for(j=0; j<num_iters; ++j){
                create_random_permutation(node_order[j], total_num_nodes, total_num_nodes);
            }
            fp = open_file(outfile, "w");
            for(i=0; i<total_num_nodes; ++i){
                fprintf(fp, "%d,", i);
                for(j=0; j<num_iters; ++j){
                    fprintf(fp, "%d,", node_order[j][i]);
                }
                fprintf(fp, "\n");
            }
            fclose(fp);
        }
    }
}


// load adjacency list from the edge file
vector<int> load_adj_list(const string &edge_path, vector<int> &node_order, vector<vector<edge>> &adj_list){
    int node1, node2, sign;
    int map_node1, map_node2;
    char str[1025], *tok;
    int num_nodes = adj_list.size();

    vector<int> degree(num_nodes, 0);

    FILE *fp = open_file(edge_path, "r");
    // skip the first line of data
    char *temp = fgets(str, 1024, fp);
    while(fgets(str, 1024, fp) != NULL){
        // 1st node --> node1
        // ascii(string) to int
        tok = strtok(str, ",");
        node1 = atoi(tok);
        // 2nd node --> node2
        tok = strtok(NULL, ",");
        node2 = atoi(tok);
        if(node1 == node2) continue;

        map_node1 = node_order[node1];
        map_node2 = node_order[node2];
        if(map_node1 >= num_nodes || map_node2 >= num_nodes) continue;

        // sign
        tok = strtok(NULL, ",");
        sign = atoi(tok);

        adj_list[map_node1].push_back({ map_node2, static_cast<int8_t>(sign) });
        degree[map_node1] += 1;
        adj_list[map_node2].push_back({ map_node1, static_cast<int8_t>(sign) });
        degree[map_node2] += 1;
    }
    fclose(fp);

    for(int i=0; i<num_nodes; ++i){
        sort(adj_list[i].begin(), adj_list[i].end(), edge_comp);
    }

    return degree;
}


pair<int, int> calc_w_sum_abs(const vector<edge> &edges_1, const vector<edge> &edges_2){
    int i = 0, j = 0;
    int w_pos = 0, w_neg = 0;
    int edge1_len = edges_1.size(), edge2_len = edges_2.size();
    while(i < edge1_len && j < edge2_len){
        if(edges_1[i].neighbor == edges_2[j].neighbor){
            if(edges_1[i].sign == edges_2[j].sign){
                ++w_pos;
            } else{
                ++w_neg;
            }
            ++i;
            ++j;
        } else if(edges_1[i].neighbor < edges_2[j].neighbor){
            ++i;
        } else{
            ++j;
        }
    }
    int w_sum = w_pos + w_neg;
    int w_abs = abs(w_pos - w_neg);
    return {w_sum, w_abs};
}


pair<int, int> calc_max_wedge(const vector<vector<edge>> &adj_list, const vector<int> &degree){
    int num_nodes = adj_list.size();

    vector<int> sorted_nodes(num_nodes);
    iota(sorted_nodes.begin(), sorted_nodes.end(), 0);
    sort(sorted_nodes.begin(), sorted_nodes.end(), [&](int a, int b){
        return degree[a] > degree[b];
    });

    int w_sum, w_abs;
    int max_w_sum = 0, max_w_diff = 0;
    vector<int> S;
    S.reserve(num_nodes);

    int u;

    for(int rank=0; rank<num_nodes; ++rank){
        u = sorted_nodes[rank];
        if(degree[u] <= max_w_sum && degree[u] <= (max_w_diff / 2)){
            break;
        }
        for(int v: S){
            pair<int, int> result = calc_w_sum_abs(adj_list[u], adj_list[v]);
            w_sum = result.first;
            w_abs = result.second;
            max_w_sum = max(max_w_sum, w_sum);
            max_w_diff = max(max_w_diff, 2 * w_abs);
        }
        S.push_back(u);
    }
    return {max_w_sum, max_w_diff};
}


void count_triangles_intersection_mixed(
    int i, int j, const vector<edge> &adj_i, const vector<edge> &adj_j, 
    int8_t sign_ij, long long &balanced_count, long long &unbalanced_count){
    // 1) find the start position of neighbors >= j
    //    [adj_i.begin()..p_end) and [adj_j.begin()..q_end) are the neighbors < j
    auto p_end = std::lower_bound(adj_i.begin(), adj_i.end(), j, edge_less_int);
    auto q_end = std::lower_bound(adj_j.begin(), adj_j.end(), j, edge_less_int);

    // #neighbors < j
    int len_i = (int)std::distance(adj_i.begin(), p_end);
    int len_j = (int)std::distance(adj_j.begin(), q_end);

    // 2) decide the way to do the intersection
    //   if one side is much shorter than the other, use binary search
    //   otherwise, use two pointers to merge scan
    const double ratio = 8.0;
    int ni, nj;
    int8_t sign_ik, sign_jk;
    if ((double)len_i < (double)len_j / ratio) {
        // use the binary search
        for(int p = 0; p < len_i; ++p){
            ni = adj_i[p].neighbor;
            sign_ik = adj_i[p].sign;

            auto it = std::lower_bound(adj_j.begin(), adj_j.begin()+len_j, ni, edge_less_int);

            if(it != (adj_j.begin()+len_j) && it->neighbor == ni){
                sign_jk = it->sign;
                if (sign_ij * sign_ik * sign_jk == 1) {
                    ++balanced_count;
                } else {
                    ++unbalanced_count;
                }
            }
        }
    } else if ((double)len_j < (double)len_i / ratio) {
        // use the binary search
        for(int q = 0; q < len_j; ++q){
            nj = adj_j[q].neighbor;
            sign_jk = adj_j[q].sign;

            auto it2 = std::lower_bound(adj_i.begin(), adj_i.begin()+len_i, nj, edge_less_int);

            if(it2 != (adj_i.begin()+len_i) && it2->neighbor == nj){
                sign_ik = it2->sign;
                if (sign_ij * sign_ik * sign_jk == 1) {
                    ++balanced_count;
                } else {
                    ++unbalanced_count;
                }
            }
        }
    } else {
        // use two pointers to merge scan
        int p = 0, q = 0;

        while(p < len_i && q < len_j){
            ni = adj_i[p].neighbor;
            nj = adj_j[q].neighbor;
            if (ni == nj) {
                sign_ik = adj_i[p].sign;
                sign_jk = adj_j[q].sign;
                if (sign_ij * sign_ik * sign_jk == 1) {
                    ++balanced_count;
                } else {
                    ++unbalanced_count;
                }
                ++p;
                ++q;
            }
            else if (ni < nj) {
                ++p;
            }
            else {
                ++q;
            }
        }
    }
}


void calc_triangles(const vector<vector<edge>> &adj_list, long long &balan_count, long long &unbalan_count){
    balan_count = 0; // balanced triangles
    unbalan_count = 0; // unbalanced triangles

    int i, j;
    int8_t sign_ij;
    int num_nodes = adj_list.size();

    for(i=0; i<num_nodes; ++i){
        // make sure i > j
        const auto &nbr_i = adj_list[i];
        for(auto &edges : nbr_i){
            j = edges.neighbor;
            if(j >= i) break;
            sign_ij = edges.sign;

            count_triangles_intersection_mixed(i, j, nbr_i, adj_list[j], sign_ij, balan_count, unbalan_count);
        }
    }
}


// calculate the smooth upper bound of local sensitivity
double calculate_SU(int LS_1, int LS_2, int num_nodes, double beta){
    double SS = 0.0;
    double beta_inverse = 1.0 / beta;
    double s_optimal;
    int s_floor, s_ceil;
    // based on the LS_1
    s_optimal = beta_inverse - (double)LS_1;
    if(s_optimal <= 0.0){
        SS = std::max(SS, (double)LS_1);
    } else if(s_optimal <= (double)(num_nodes-2-LS_1)){
        s_floor = floor(s_optimal);
        SS = std::max(SS, exp(-1.0*beta*(double)s_floor) * (double)(LS_1 + s_floor));
        s_ceil = ceil(s_optimal);
        SS = std::max(SS, exp(-1.0*beta*(double)s_ceil) * (double)(LS_1 + s_ceil));
    } else{
        SS = std::max(SS, exp(-1.0*beta*(double)(num_nodes-2-LS_1)) * (double)(num_nodes-2));
    }
    // based on the LS_2
    s_optimal = beta_inverse - (double)LS_2 / 4.0;
    if(s_optimal <= 0.0){
        SS = std::max(SS, (double)LS_2);
    } else if(s_optimal <= (double)(2*num_nodes-4-LS_2)/4.0){
        s_floor = floor(s_optimal);
        SS = std::max(SS, exp(-1.0*beta*(double)s_floor) * (double)(LS_2 + 4 * s_floor));
        s_ceil = ceil(s_optimal);
        SS = std::max(SS, exp(-1.0*beta*(double)s_ceil) * (double)(LS_2 + 4 * s_ceil));
    } else{
        SS = std::max(SS, exp(-1.0*beta*(double)(2*num_nodes-4-LS_2)/4.0) * (double)(2*num_nodes-4));
    }
    return SS;
}


int main(int argc, char *argv[]){
    if(argc < 2){
        cout << "Usage: " << argv[0] << " [dataset] ([eps (default: 0.1)] [mech (default: 3)] [sample nodes (default: 1.0)] [iters (default: 1)] [repeats (default: 100)])" << endl;
        cout << "[dataset]: dataset name" << endl;
        cout << "[eps]: privacy budget" << endl;
        cout << "[mech]: mechanism (0: GS 1: SS (baseline) 2: SS (improved) 3: SS (upper bound)) but only work for mech=0 or 3" << endl;
        cout << "[sample nodes]: percentage of sample nodes (1.0: total nodes)" << endl;
        cout << "[iters]: number of iterations (works when the nodes are sampled)" << endl;
        cout << "[repeats]: number of repeats" << endl;
        exit(-1);
    }
    // initialization of Mersennne Twister
    unsigned long init[4] = { 0x123, 0x234, 0x345, 0x456 }, length = 4;
    init_by_array(init, length);

    string dataset = argv[1];
    // default parameters
    double eps = 0.1;
    int mech = 3;
    double sample_percentage = 1.0;
    int num_iters = 1;
    int repeats = 100;

    string edge_path = "./data/" + dataset + "_undirect.csv";
    int total_num_nodes;
    int num_sample_nodes;
    double delta;
    vector<vector<int>> node_order;
    vector<vector<edge>> adj_list;
    vector<int> degree;
    // number of balanced and unbalanced triangles
    vector<long long> num_triangles(2);
    vector<double> num_ns_triangles(2);
    double absolute_error_l1 = 0.0, absolute_error_l1_avg = 0.0;
    double relative_error_l1 = 0.0, relative_error_l1_avg = 0.0;

    // assign the values of command line arguments to variables
    if(argc > 2){
        eps = atof(argv[2]);
        if(argc > 3){
            mech = atoi(argv[3]);
            if(mech != 0 && mech != 3){
                cout << "mech should be 0 or 3" << endl;
                exit(-1);
            }
            if(argc > 4){
                sample_percentage = atof(argv[4]);
                if(argc > 5){
                    num_iters = atoi(argv[5]);
                    if(argc > 6) repeats = atoi(argv[6]);
                }
            }
        }
    }
    // get the total number of nodes --> num_nodes
    total_num_nodes = total_nodes_count(edge_path);
    num_sample_nodes = (int)(sample_percentage * (double)total_num_nodes);
    create_node_order(dataset, num_sample_nodes, total_num_nodes, num_iters, node_order);
    cout << "dataset: " << dataset << " eps: " << eps << " mech: " << mech << " total_num_nodes: " << total_num_nodes << " sample_percentage: " << sample_percentage << " num_iters: " << num_iters << " repeats: " << repeats << endl;
    delta = 1.0 / (10.0 * (double)(num_sample_nodes) * (double)(num_sample_nodes - 1) / 2.0);

    for(int i=0; i<num_iters; ++i){
        // initialization
        adj_list = vector<vector<edge>>(num_sample_nodes);
        // load edges from the edge file --> adj_list and degree
        degree = load_adj_list(edge_path, node_order[i], adj_list);
        // calculate the true number of wedges and triangles
        calc_triangles(adj_list, num_triangles[0], num_triangles[1]);
        
        // calculate the noisy number of balanced and unbalanced triangles
        if(mech == 0){
            double GS;
            GS = double(2 * num_sample_nodes - 4);
            // cout << "mech: " << mech << " GS: " << GS << endl;
            for(int repeat=0; repeat<repeats; ++repeat){
                if(genrand_real2() >= delta){
                    num_ns_triangles[0] = (double)num_triangles[0] + stats::rlaplace(0.0, GS/eps, engine);
                    num_ns_triangles[1] = (double)num_triangles[1] + stats::rlaplace(0.0, GS/eps, engine);
                } else{ // probability of privacy breach
                    num_ns_triangles[0] = (double)num_triangles[0];
                    num_ns_triangles[1] = (double)num_triangles[1];
                }
                // post-processing
                num_ns_triangles[0] = std::max(0.0, num_ns_triangles[0]);
                num_ns_triangles[1] = std::max(0.0, num_ns_triangles[1]);
                absolute_error_l1 = fabs(num_ns_triangles[0]-(double)num_triangles[0]) + fabs(num_ns_triangles[1]-(double)num_triangles[1]);
                relative_error_l1 = absolute_error_l1 / ((double)num_triangles[0] + (double)num_triangles[1]);
                absolute_error_l1_avg += absolute_error_l1;
                relative_error_l1_avg += relative_error_l1;
            }
        } else{
            int LS_1=0, LS_2=0;
            double SS = 0.0;
            double beta = eps / (4.0*(2.0+log(2.0/delta)));

            pair<int, int> result = calc_max_wedge(adj_list, degree);
            LS_1 = result.first;
            LS_2 = result.second;
            SS = calculate_SU(LS_1, LS_2, adj_list.size(), beta);
            // cout << "mech: " << mech << " LS_1: " << LS_1 << " LS_2: " << LS_2 << " SS: " << SS << endl;
            for(int repeat=0; repeat<repeats; ++repeat){
                // num_ns_triangles[0] = num_triangles[0] + stats::rcauchy(0.0, sqrt(2.0)*SS/eps, engine);
                num_ns_triangles[0] = num_triangles[0] + stats::rlaplace(0.0, 2.0*SS/eps, engine);
                // num_ns_triangles[1] = num_triangles[1] + stats::rcauchy(0.0, sqrt(2.0)*SS/eps, engine);
                num_ns_triangles[1] = num_triangles[1] + stats::rlaplace(0.0, 2.0*SS/eps, engine);
                num_ns_triangles[0] = std::max(0.0, num_ns_triangles[0]);
                num_ns_triangles[1] = std::max(0.0, num_ns_triangles[1]);
                absolute_error_l1 = fabs(num_ns_triangles[0]-(double)num_triangles[0]) + fabs(num_ns_triangles[1]-(double)num_triangles[1]);
                relative_error_l1 = absolute_error_l1 / ((double)num_triangles[0] + (double)num_triangles[1]);
                absolute_error_l1_avg += absolute_error_l1;
                relative_error_l1_avg += relative_error_l1;
            }
        }
        // release memory
        vector<int>().swap(degree);
        vector<vector<edge>>().swap(adj_list);
    }
    cout << "Average Absolute Error: " << absolute_error_l1_avg/(double)(num_iters*repeats) << " Average Relative Error: " << relative_error_l1_avg/(double)(num_iters*repeats) << endl;
    cout << "**********************************************************************" << endl;
    vector<vector<int>>().swap(node_order);
    return 0;
}
