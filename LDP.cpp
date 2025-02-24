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
stats::rand_engine_t engine(1776);


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


// randomly generate 0, 1, ..., size-1, and store the first 'num' values into rndperm
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
vector<int> load_adj_list(const string &edge_path, vector<int> &node_order, vector<flat_hash_map<int, int8_t>> &adj_list, vector<int> &degree){
    int node1, node2, sign;
    int map_node1, map_node2;
    char str[1025], *tok;
    int num_nodes = adj_list.size();
    vector<int> ldegree(num_nodes, 0);

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

        adj_list[map_node1][map_node2] = static_cast<int8_t>(sign);
        ++degree[map_node1];
        adj_list[map_node2][map_node1] = static_cast<int8_t>(sign);
        ++degree[map_node2];
        // calculate the degree from the lower triangular part of the adjacency matrix
        if(map_node1 > map_node2) ++ldegree[map_node1];
        else ++ldegree[map_node2];
    }
    fclose(fp);
    return ldegree;
}


void calc_triangles(vector<flat_hash_map<int, int8_t>> &adj_list, vector<long long> &num_triangles){
    int i, j, k;
    int8_t sign_ij, sign_ik, sign_jk;
    flat_hash_map<int, int8_t>::iterator iter_1, iter_2, iter_3;

    int num_nodes = adj_list.size();
    for(i=0; i<num_nodes; ++i){
        for(iter_1 = adj_list[i].begin(); iter_1 != adj_list[i].end(); ++iter_1){
            j = iter_1->first;
            // i > j > k
            if(i <= j) continue;
            sign_ij = iter_1->second;
            for(iter_2 = adj_list[i].begin(); iter_2 != adj_list[i].end(); ++iter_2){
                k = iter_2->first;
                if(j <= k) continue;
                sign_ik = iter_2->second;

                iter_3 = adj_list[j].find(k);
                if(iter_3 != adj_list[j].end()){
                    sign_jk = iter_3->second;
                    if((sign_ij * sign_ik * sign_jk) == 1) num_triangles[0] += 1;
                    else num_triangles[1] += 1;
                }
            }
        }
    }
}


// calculate the noisy max degree
void calc_noisy_max_degree(vector<int> &degree, double eps, int &max_degree, double &ns_max_degree){
    double ns_degree;
    max_degree = 0;
    ns_max_degree = 0.0;

    // calculate the max of noisy degree --> ns_max_degree (Laplace mechanism)
    int num_nodes = degree.size();
    for(int i=0; i<num_nodes; ++i){
        if(max_degree < degree[i]) max_degree = degree[i];
        ns_degree = (double)degree[i] + stats::rlaplace(0.0, 1.0/eps, engine);
        if(ns_max_degree < ns_degree) ns_max_degree = ns_degree;
    }

    // if ns_max_degree exceeds num_nodes - 1, then use num_nodes - 1
    if(ns_max_degree > (double)num_nodes - 1.0){
        ns_max_degree = (double)num_nodes - 1.0;
    }
}


// GRR mechanism
inline int8_t GRR(int8_t v, double q) {
    double rand_v = genrand_real1();
    if(rand_v >= q) return v;
    double half = q / 2.0;
    switch (v) {
        case 0:
            return rand_v < half ? 1 : -1;
        case 1:
            return rand_v < half ? 0 : -1;
        case -1:
            return rand_v < half ? 0 : 1;
        default:
            return 0;
    }
}


// calculate #balanced and #unbalanced triangles in the non-interactive local model
void calc_triangle_non_interactive(vector<flat_hash_map<int, int8_t>> &adj_list, double eps, bool estimate_flag, vector<double> &num_ns_triangles){
    int num_nodes = adj_list.size();
    // noisy adjacency list
    vector<flat_hash_map<int, int8_t>> ns_adj_list(num_nodes);
    // noisy degree
    vector<int> ns_degree(num_nodes, 0);
    long long num_ns_balanced_triangles = 0;
    long long num_ns_unbalanced_triangles = 0;
    // flip probability
    double filp_prob;
    int i, j, k;
    int8_t sign_ij, sign_ik, sign_jk;
    flat_hash_map<int, int8_t>::iterator iter_1, iter_2, iter_3;

    // the flip probability
    filp_prob = 2.0 / (exp(eps) + 2.0);

    // flip 0/1 in adj_list with probability q --> ns_adj_list
    for(i=0; i<num_nodes; ++i){
        auto &cur_adj = adj_list[i];
        for(j=0; j<i; ++j){
            iter_1 = cur_adj.find(j);
            if(iter_1 == cur_adj.end()) sign_ij = GRR(0, filp_prob);
            else sign_ij = GRR(iter_1->second, filp_prob);
            if(sign_ij != 0){
                ns_adj_list[i][j] = sign_ij;
                ++ns_degree[i];
                ++ns_degree[j];
            }
        }
    }

    // calculate #balanced and unbalanced triangles in the noisy graph
    for(i=0; i<num_nodes; ++i){
        for(iter_1 = ns_adj_list[i].begin(); iter_1 != ns_adj_list[i].end(); ++iter_1){
            j = iter_1->first;
            // i > j > k
            if(i<=j) continue;
            sign_ij = iter_1->second;
            for (iter_2 = ns_adj_list[i].begin(); iter_2 != ns_adj_list[i].end(); ++iter_2) {
                k = iter_2->first;
                if(j<=k) continue;
                sign_ik = iter_2->second;
                iter_3 = ns_adj_list[j].find(k);
                if(iter_3 != ns_adj_list[j].end()){
                    sign_jk = iter_3->second;
                    if((sign_ij * sign_ik * sign_jk) == 1) ++num_ns_balanced_triangles;
                    else ++num_ns_unbalanced_triangles;
                }
            }
        }
    }

    // empirical estimation
    if(estimate_flag){
        long long num_ns_2_stars = 0;
        long long num_ns_edges = 0;
        for(i=0; i<num_nodes; ++i){
            num_ns_edges += (long long)ns_degree[i];
            num_ns_2_stars += ((long long)ns_degree[i] * ((long long)ns_degree[i]-1)) / 2;
        }

        // calculate the number of 2-edges in the noisy graph --> num_ns_2_edges
        long long num_ns_2_edges = num_ns_2_stars - 3 * (num_ns_balanced_triangles+num_ns_unbalanced_triangles);
        // calculate the number of 1-edge in the noisy graph --> num_ns_1_edges
        long long num_ns_1_edges = num_ns_edges*((long long)num_nodes-2) - 2*num_ns_2_edges - 3 * (num_ns_balanced_triangles+num_ns_unbalanced_triangles);
        // calculate the number of 0-edge in the noisy graph --> num_ns_0_edges
        long long num_ns_0_edges = (long long)(num_nodes*(num_nodes-1)*(num_nodes-2)/6) - (num_ns_balanced_triangles+num_ns_unbalanced_triangles) - num_ns_2_edges - num_ns_1_edges;

        double e_eps = exp(eps);
        double e_eps_2 = e_eps * e_eps;
        double e_eps_plus_1_2 = (e_eps+1.0) * (e_eps+1.0);
        double e_eps_plus_1_3 = (e_eps+1.0) * (e_eps+1.0) * (e_eps+1.0);
        double e_eps_minus_1_3 = (e_eps-1.0) * (e_eps-1.0) * (e_eps-1.0);

        // balanced triangles
        num_ns_triangles[0] = (e_eps_plus_1_3+3.0*(e_eps+1.0))*(double)num_ns_balanced_triangles - (3.0*e_eps_plus_1_2+1.0)*(double)num_ns_unbalanced_triangles - e_eps_2*(double)num_ns_2_edges + 2.0*e_eps*(double)num_ns_1_edges - 4.0*(double)num_ns_0_edges;
        num_ns_triangles[0] /= e_eps_minus_1_3;
        // unbalanced triangles
        num_ns_triangles[1] = (e_eps_plus_1_3+3.0*(e_eps+1.0))*(double)num_ns_unbalanced_triangles - (3.0*e_eps_plus_1_2+1.0)*(double)num_ns_balanced_triangles - e_eps_2*(double)num_ns_2_edges + 2.0*e_eps*(double)num_ns_1_edges - 4.0*(double)num_ns_0_edges;
        num_ns_triangles[1] /= e_eps_minus_1_3;
    }
    // without empirical estimation
    else{
        num_ns_triangles[0] = (double)num_ns_balanced_triangles;
        num_ns_triangles[1] = (double)num_ns_unbalanced_triangles;
    }
    vector<int>().swap(ns_degree);
    vector<flat_hash_map<int, int8_t>>().swap(ns_adj_list);
}


// calculate #balanced and unbalanced triangles in the interactive local model
void calc_triangles_interactive(vector<flat_hash_map<int, int8_t>> &adj_list, vector<int> &degree, vector<int> &ldegree, bool SS_tag, double eps_degree, double eps_adj, double eps_triangle, double delta, vector<double> &num_ns_triangles){
    int num_nodes = adj_list.size();
    // noisy adjacency list
    vector<flat_hash_map<int, int8_t>> ns_adj_list(num_nodes);
    // #balanced and unbalanced triangles for each user
    vector<long long> num_ns_b_triangles_user(num_nodes, 0);
    vector<long long> num_ns_u_triangles_user(num_nodes, 0);
    vector<long long> num_2_stars_user(num_nodes, 0);
    double num_ns_b_triangles, num_ns_u_triangles;
    double filp_prob, q;
    int max_degree, ns_max_degree_floor;
    double ns_max_degree;
    int i, j, k, index;
    int8_t sign_ij, sign_ik, sign_jk;
    flat_hash_map<int, int8_t>::iterator iter_1, iter_2, iter_3, iter_4;

    vector<int> rndperm;
    double beta = eps_triangle / (4.0*(2.0+log(2.0/delta)));

    // flip probability
    q = 1.0 / (exp(eps_adj) + 2.0);
    filp_prob = 2.0 * q;

    for(i=0; i<num_nodes; ++i){
        for(iter_1=adj_list[i].begin(); iter_1!=adj_list[i].end(); ++iter_1){
            j = iter_1->first;
            // it must be i > j > k
            if(i<=j) continue;
            sign_ij = iter_1->second;
            for(iter_2=adj_list[i].begin(); iter_2!=adj_list[i].end(); ++iter_2){
                k = iter_2->first;
                if(j<=k) continue;
                sign_ik = iter_2->second;
                num_2_stars_user[i] += 1;

                // GRR --> ns_adj_list
                iter_3 = ns_adj_list[j].find(k);
                // if ns_adj_list[j][k] does not exist
                if(iter_3 == ns_adj_list[j].end()){
                    iter_4 = adj_list[j].find(k);
                    if(iter_4 == adj_list[j].end()) sign_jk = GRR(0, filp_prob);
                    else sign_jk = GRR(iter_4->second, filp_prob);
                    ns_adj_list[j][k] = sign_jk;
                } else{
                    sign_jk = iter_3->second;
                }
                if(sign_ij * sign_ik * sign_jk == 1){
                    num_ns_b_triangles_user[i] += 1;
                } else if (sign_ij * sign_ik * sign_jk == -1){
                    num_ns_u_triangles_user[i] += 1;
                }
            }
        }
    }

    if(SS_tag == false){
        // deleted adjacency list after projection
        vector<flat_hash_map<int, int8_t>> deleted_adj_list(num_nodes);
        // noisy max degree: degree + Lap
        // use the noisy max degree to calculate the global sensitivity
        calc_noisy_max_degree(degree, eps_degree, max_degree, ns_max_degree);

        // if max_degree exceeds ns_max_degree, then perform graph projection
        if((double)max_degree > ns_max_degree){
            ns_max_degree_floor = (int)floor(ns_max_degree);
            for(i=0; i<num_nodes; ++i){
                // if degree[i] exceeds ns_max_degree, then perform graph projection
                if(degree[i] > ns_max_degree){
                    // randomly generate 0, 1, ..., degree[i]-1 --> rndperm
                    rndperm = vector<int>(degree[i]);
                    create_random_permutation(rndperm, degree[i], degree[i]);

                    // randomly delete (degree[i] - max_deg_ns) edges from adj_list[i]
                    index = 0;
                    for(iter_1=adj_list[i].begin(); iter_1!=adj_list[i].end(); ++iter_1){
                        if(rndperm[index] >= ns_max_degree_floor){
                            j = iter_1->first;
                            // deleted edge --> deleted_adj_list[i][j]
                            // positive edge or negative edge
                            deleted_adj_list[i][j] = 1;
                        }
                        ++index;
                    }
                    vector<int>().swap(rndperm);

                    num_ns_b_triangles_user[i] = 0;
                    num_ns_u_triangles_user[i] = 0;
                    num_2_stars_user[i] = 0;
                    for(iter_1=adj_list[i].begin(); iter_1!=adj_list[i].end(); ++iter_1){
                        j = iter_1->first;
                        // it must be i > j > k
                        // continue if the edge is deleted
                        if(i <= j || deleted_adj_list[i].count(j) != 0) continue;
                        sign_ij = iter_1->second;
                        for(iter_2=adj_list[i].begin(); iter_2!=adj_list[i].end(); ++iter_2){
                            k = iter_2->first;
                            // continue if the edge is deleted
                            if(j <= k || deleted_adj_list[i].count(k) != 0) continue;
                            sign_ik = iter_2->second;
                            num_2_stars_user[i] += 1;
                            sign_jk = ns_adj_list[j][k];
                            if(sign_ij * sign_ik * sign_jk == 1){
                                num_ns_b_triangles_user[i] += 1;
                            } else if(sign_ij * sign_ik * sign_jk == -1){
                                num_ns_u_triangles_user[i] += 1;
                            }
                        }
                    }
                }
            }
        }

        vector<flat_hash_map<int, int8_t>>().swap(deleted_adj_list);

        double GS = std::max(ns_max_degree, 2.0*(ns_max_degree-1.0));
        for(i=0; i<num_nodes; ++i){
            // balanced triangles
            num_ns_b_triangles = (double)num_ns_b_triangles_user[i] - q*(double)num_2_stars_user[i];
            // unbalanced triangles
            num_ns_u_triangles = (double)num_ns_u_triangles_user[i] - q*(double)num_2_stars_user[i];
            if(genrand_real2() >= delta){
                num_ns_b_triangles += stats::rlaplace(0.0, GS/eps_triangle, engine);
                num_ns_u_triangles += stats::rlaplace(0.0, GS/eps_triangle, engine);
            }
            num_ns_triangles[0] += num_ns_b_triangles;
            num_ns_triangles[1] += num_ns_u_triangles;
        }
    } else{
        int ldegree_user, LS_s_upper_bound;
        double SS;
        double s_optimal;
        int s_floor, s_ceil;
        // calculate the smooth upper bound of local sensitivity
        for(i=0; i<num_nodes; ++i){
            ldegree_user = ldegree[i];
            // distance s = 0
            LS_s_upper_bound = std::max(2*(ldegree_user-1), ldegree_user);
            SS = (double)LS_s_upper_bound;
            // distance s = 1
            LS_s_upper_bound = std::max(std::min(2*ldegree_user, 2*(i-1)), std::min(ldegree_user+1, i));
            SS = std::max(SS, exp(-1.0*beta)*(double)LS_s_upper_bound);
            // distance s >= 2
            s_optimal = 1.0/beta - (double)ldegree_user + 1.0;
            if(s_optimal < double(i-ldegree_user)){
                if(s_optimal > 2.0){
                    s_floor = floor(s_optimal);
                    SS = std::max(SS, exp(-1.0*beta*(double)s_floor) * 2.0 * (double)(ldegree_user+s_floor-1));
                    s_ceil = ceil(s_optimal);
                    SS = std::max(SS, exp(-1.0*beta*(double)s_ceil) * 2.0 * (double)(ldegree_user+s_ceil-1));
                } else if(1.0 < s_optimal && s_optimal <= 2.0){
                    SS = std::max(SS, exp(-2.0*beta) * 2.0 * (double)(ldegree_user+1));
                }
            } else{
                SS = std::max(SS, exp(-1.0*beta*(double)(i-ldegree_user)) * 2.0 * (double)(i-1));
            }
            num_ns_b_triangles = (double)num_ns_b_triangles_user[i] - q*(double)num_2_stars_user[i];
            // num_ns_b_triangles = (double)num_ns_b_triangles_user[i];
            // add the noise for each user
            num_ns_b_triangles +=  stats::rlaplace(0.0, 2.0*SS/eps_triangle, engine);

            // unbalanced triangles
            num_ns_u_triangles = (double)num_ns_u_triangles_user[i] - q*(double)num_2_stars_user[i];
            // num_ns_u_triangles = (double)num_ns_u_triangles_user[i];
            num_ns_u_triangles += stats::rlaplace(0.0, 2.0*SS/eps_triangle, engine);

            num_ns_triangles[0] += num_ns_b_triangles;
            num_ns_triangles[1] += num_ns_u_triangles;
        }
    }
    // empirical estimate
    num_ns_triangles[0] /= (1.0 - 3.0 * q);
    num_ns_triangles[1] /= (1.0 - 3.0 * q);
    // post-processing
    num_ns_triangles[0] = std::max(0.0, num_ns_triangles[0]);
    num_ns_triangles[1] = std::max(0.0, num_ns_triangles[1]);

    vector<long long>().swap(num_2_stars_user);
    vector<long long>().swap(num_ns_u_triangles_user);
    vector<long long>().swap(num_ns_b_triangles_user);
    vector<flat_hash_map<int, int8_t>>().swap(ns_adj_list);
}


int main(int argc, char *argv[]){
    if(argc < 2){
        cout << "Usage: " << argv[0] << " [dataset] ([eps (default: 1.0)] [eps-alloc (default: 2-9-9)]) [mech (default: 3)] [sample nodes (default: 1.0)] [iters (default: 1)] [repeats (default: 100)])" << endl;
        cout << "[dataset]: dataset name" << endl;
        cout << "[eps]: privacy budget" << endl;
        cout << "[eps-alloc]: privacy budget allocation: degree-adj-triangle" << endl;
        cout << "[mech]: mechanism (0: non-interactive local model (w/o empirical estimation) 1: non-interactive local model (w/ empirical estimation) 2: interactive local model (GS projection), 3: interactive local model (SS))" << endl;
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
    double total_eps = 1.0;
    int mech = 3;
    double sample_percentage = 1.0;
    int num_iters = 1;
    int repeats = 100;
    double eps_alloc[3] = {2.0, 9.0, 9.0};

    string edge_path = "./data/" + dataset + "_undirect.csv";
    int total_num_nodes, num_sample_nodes;
    double delta;
    char *tok;
    vector<vector<int>> node_order;
    vector<flat_hash_map<int, int8_t>> adj_list;
    vector<int> degree;
    vector<int> ldegree;
    // #balanced and unbalanced triangles
    vector<long long> num_triangles(2);
    vector<double> num_ns_triangles(2);
    double absolute_error_l1 = 0.0, absolute_error_l1_avg = 0.0;
    double relative_error_l1 = 0.0, relative_error_l1_avg = 0.0;

    // assign the values of command line arguments to variables
    if(argc > 2){
        total_eps = atof(argv[2]);
        if(argc > 3){
            if((tok = strtok(argv[3], "-")) == NULL){
                cout << "Error: incorect [eps-alloc]" << endl;
                exit(-1);
            }
            eps_alloc[0] = atof(tok);
            if((tok = strtok(NULL, "-")) == NULL){
                cout << "Error: incorect [eps-alloc]" << endl;
                exit(-1);
            }
            eps_alloc[1] = atof(tok);
            if((tok = strtok(NULL, "-")) == NULL){
                cout << "Error: incorect [eps-alloc]" << endl;
                exit(-1);
            }
            eps_alloc[2] = atof(tok);
            if(argc > 4){
                mech = atoi(argv[4]);
                if(argc > 5){
                    sample_percentage = atof(argv[5]);
                    if(argc > 6){
                        num_iters = atoi(argv[6]);
                    if(argc > 7) repeats = atoi(argv[7]);
                    }
                }
            }
        }
    }
    // get #total nodes --> num_nodes
    total_num_nodes = total_nodes_count(edge_path);
    num_sample_nodes = (int)(sample_percentage * (double)total_num_nodes);
    create_node_order(dataset, num_sample_nodes, total_num_nodes, num_iters, node_order);
    double eps_degree, eps_adj, eps_triangle;
    if(mech == 0 || mech == 1){
        eps_adj = total_eps;
    } else if(mech == 2){
        eps_degree = total_eps / (eps_alloc[0]+eps_alloc[1]+eps_alloc[2]) * eps_alloc[0];
        eps_adj = total_eps / (eps_alloc[0]+eps_alloc[1]+eps_alloc[2]) * eps_alloc[1];
        eps_triangle = total_eps / (eps_alloc[0]+eps_alloc[1]+eps_alloc[2]) * eps_alloc[2];
    } else{
        eps_adj = total_eps / (eps_alloc[1]+eps_alloc[2]) * eps_alloc[1];
        eps_triangle = total_eps / (eps_alloc[1]+eps_alloc[2]) * eps_alloc[2];
    }
    delta = 1.0 / (10.0 * (double)(num_sample_nodes));

    cout << "dataset: " << dataset << " total eps: " << total_eps << " " << eps_alloc[0] << "-" << eps_alloc[1] << "-" << eps_alloc[2] << " mech: " << mech << " total_num_nodes: " << total_num_nodes << " sample_percentage: " << sample_percentage << " num_iters: " << num_iters << " repeats: " << repeats << endl;

    for(int i=0; i<num_iters; ++i){
        // Initialization
        adj_list = vector<flat_hash_map<int, int8_t>>(num_sample_nodes);
        degree = vector<int>(num_sample_nodes);
        num_triangles[0] = 0;
        num_triangles[1] = 0;
        // load edges from the edge file --> adj_list
        ldegree = load_adj_list(edge_path, node_order[i], adj_list, degree);
        calc_triangles(adj_list, num_triangles);

        for(int repeat=0; repeat<repeats; ++repeat){
            num_ns_triangles[0] = 0.0;
            num_ns_triangles[1] = 0.0;
            if(mech == 0){
                calc_triangle_non_interactive(adj_list, eps_adj, false, num_ns_triangles);
            } else if(mech == 1){
                calc_triangle_non_interactive(adj_list, eps_adj, true, num_ns_triangles);
            } else if(mech == 2){
                calc_triangles_interactive(adj_list, degree, ldegree, false, eps_degree, eps_adj, eps_triangle, delta, num_ns_triangles);
            } else{
                calc_triangles_interactive(adj_list, degree, ldegree, true, eps_degree, eps_adj, eps_triangle, delta, num_ns_triangles);
            }
            absolute_error_l1 = fabs(num_ns_triangles[0]-(double)num_triangles[0]) + fabs(num_ns_triangles[1]-(double)num_triangles[1]);
            relative_error_l1 = absolute_error_l1 / ((double)num_triangles[0] + (double)num_triangles[1]);
            absolute_error_l1_avg += absolute_error_l1;
            relative_error_l1_avg += relative_error_l1;
        }
        vector<int>().swap(degree);
        vector<int>().swap(ldegree);
        vector<flat_hash_map<int, int8_t>>().swap(adj_list);
    }
    cout << "Average Absolute Error: " << absolute_error_l1_avg/(double)(num_iters*repeats) << " Average Relative Error: " << relative_error_l1_avg/(double)(num_iters*repeats) << endl;
    cout << "**********************************************************************" << endl;
    vector<vector<int>>().swap(node_order);
    return 0;
}
