#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <unordered_map>
#include <algorithm>

#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#include "include/stats.hpp"
#include "include/mt19937ar.h"


using namespace std;

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


int get_total_num_nodes(const string &edge_path){
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
void generate_random_permutation(vector<int> &rndperm, int size, int num){
    int rnd;
    vector<int> ordperm(size);
    // 0, 1, 2, ..., size-1 --> ordperm
    for(int i=0; i<size; ++i){
        ordperm[i] = i;
    }
    for(int i=0; i<num; ++i){
        rnd = genrand_int32() % (size - i);
        rndperm[i] = ordperm[rnd];
        for(int j=rnd+1; j<size-i; ++j){
            ordperm[j - 1] = ordperm[j];
        }
    }
    vector<int>().swap(ordperm);
}


void generate_node_order(const string &dataset, int &num_sample_nodes, const int &total_num_nodes, int &num_iters, vector<vector<int>> &node_order){
    if(num_sample_nodes == total_num_nodes){ // use all nodes to calculate the number of triangles
        num_iters = 1;
        node_order = vector<vector<int>>(num_iters, vector<int>(total_num_nodes));
        for(int i=0; i<num_sample_nodes; ++i){
            node_order[0][i] = i;
        }
    } else{ // randomly generate the order of nodes --> node_order
        FILE *fp;
        string outfile = "./results/"+dataset+"/node_order_iters_"+to_string(num_iters)+".csv";
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
                generate_random_permutation(node_order[j], total_num_nodes, total_num_nodes);
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
void load_adj_list(const string &edge_path, vector<int> &node_order, vector<unordered_map<int,int>> &adj_list){
    int node1, node2, sign;
    FILE *fp;
    char str[1025];
    char *tok, *temp;

    int num_nodes = adj_list.size();
    fp = open_file(edge_path, "r");
    // skip the first line of data
    temp = fgets(str, 1024, fp);
    while(fgets(str, 1024, fp) != NULL){
        // 1st node --> node1
        // ascii(string) to int
        tok = strtok(str, ",");
        node1 = atoi(tok);
        // 2nd node --> node2
        tok = strtok(NULL, ",");
        node2 = atoi(tok);
        // sign
        tok = strtok(NULL, ",");
        sign = atoi(tok);
        if(node1 == node2) continue;
        // if both nodes exist, add the edge
        if(node_order[node1] < num_nodes && node_order[node2] < num_nodes){
            adj_list[node_order[node1]][node_order[node2]] = sign;
            adj_list[node_order[node2]][node_order[node1]] = sign;
        }
    }
    fclose(fp);
}


void calculate_degree(vector<unordered_map<int,int>> &adj_list, vector<int> &degree, vector<int> &ldegree){
    int i, j;
    unordered_map<int, int>::iterator iter_1;

    int num_nodes = adj_list.size();
    for(i=0; i<num_nodes; ++i){
        for(iter_1=adj_list[i].begin(); iter_1!=adj_list[i].end(); ++iter_1){
            j = iter_1->first;
            // degree = positive edges + negtive edges
            degree[i] += 1;
            // calculate the degree from the lower triangular part of the adjacency matrix
            if(i > j) ldegree[i] += 1;
        }
    }
}


void calculate_triangles(vector<unordered_map<int,int>> &adj_list, vector<long long> &num_triangles){
    int i, j, k, sign_ij, sign_ik, sign_jk;
    unordered_map<int, int>::iterator iter_1, iter_2;

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
                if(adj_list[j].count(k) > 0){
                    sign_jk = adj_list[j][k];
                    if((sign_ij * sign_ik * sign_jk) == 1) num_triangles[0] += 1;
                    else num_triangles[1] += 1;
                }
            }
        }
    }
}


// calculate the noisy max degree
void calculate_noisy_max_degree(vector<int> &degree, double eps, int &max_degree, double &ns_max_degree){
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


// generalized randomized response mechanism
inline int GRR(int v, double q, vector<int> &values){
    if(genrand_real1() >= q){ // not flip
        return v;
    } else{
        vector<int> candidate_values;
        for(int value : values){
            if(value != v){
                candidate_values.emplace_back(value);
            }
        }
        return candidate_values[genrand_int32() % 2];
    }
}


// calculate balanced and unbalanced triangle counts in the non-interactive local model
void calculate_triangle_non_interactive(vector<unordered_map<int, int>> &adj_list, double eps, int estimate_flag, vector<double> &num_ns_triangles){
    int num_nodes = adj_list.size();
    // noisy adjacency list
    vector<unordered_map<int, int>> ns_adj_list(num_nodes);
    // noisy degree
    vector<int> ns_degree(num_nodes, 0);
    long long num_ns_balanced_triangles = 0;
    long long num_ns_unbalanced_triangles = 0;
    long long num_ns_edges = 0;
    // flip probability
    double filp_prob;
    int i, j, k, sign_ij, sign_ik, sign_jk;
    unordered_map<int, int>::iterator iter_1, iter_2;

    // the flip probability
    filp_prob = 2.0 / (exp(eps) + 2.0);
    vector<int> values = {0, +1, -1};

    // flip 0/1 in adj_list with probability q --> ns_adj_list
    for(i=0; i<num_nodes; ++i){
        for(j=0; j<i; ++j){
            if(adj_list[i].count(j) == 0) sign_ij = GRR(0, filp_prob, values);
            else sign_ij = GRR(adj_list[i][j], filp_prob, values);
            if(sign_ij != 0){
                ns_adj_list[i][j] = sign_ij;
                ns_adj_list[j][i] = sign_ij;
            }
        }
    }
    // calculate the degree of each user in the noisy graph
    for(i=0; i<num_nodes; ++i){
        for(iter_1=ns_adj_list[i].begin(); iter_1!=ns_adj_list[i].end(); ++iter_1) ++ns_degree[i];
    }
    // calculate the total number of edges
    for(i=0; i<num_nodes; ++i) num_ns_edges += (long long)ns_degree[i];
    num_ns_edges /= 2;

    // calculate the number of balanced and unbalanced triangles in the noisy graph
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
                if(ns_adj_list[j].count(k) > 0){
                    sign_jk = ns_adj_list[j][k];
                    if((sign_ij * sign_ik * sign_jk) == 1) ++num_ns_balanced_triangles;
                    else ++num_ns_unbalanced_triangles;
                }
            }
        }
    }

    // empirical estimation
    if(estimate_flag == 1){
        long long num_ns_2_stars = 0;
        for(i=0; i<num_nodes; ++i){
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
    vector<unordered_map<int, int>>().swap(ns_adj_list);
}


// calculate balanced and unbalanced triangle counts in the interactive local model
void calculate_triangle_interactive(vector<unordered_map<int, int>> &adj_list, int SS_tag, double eps_degree, double eps_adj, double eps_triangle, double delta, vector<double> &num_ns_triangles){
    int num_nodes = adj_list.size();
    vector<int> degree(num_nodes);
    vector<int> ldegree(num_nodes);
    calculate_degree(adj_list, degree, ldegree);
    // noisy adjacency list
    vector<unordered_map<int, int>> ns_adj_list(num_nodes);
    // the number of balanced and unbalanced triangles for each user
    vector<long long> num_ns_b_triangles_user(num_nodes, 0);
    vector<long long> num_ns_u_triangles_user(num_nodes, 0);
    vector<long long> num_2_stars_user(num_nodes, 0);
    double num_ns_b_triangles, num_ns_u_triangles;
    double filp_prob, q;
    int max_degree, ns_max_degree_floor;
    double ns_max_degree;
    int i, j, k, index, sign_ij, sign_ik, sign_jk;
    unordered_map<int, int>::iterator iter_1, iter_2;
    vector<int> rndperm;
    double beta = eps_triangle / (4.0*(2.0+log(2.0/delta)));

    // flip probability
    q = 1.0 / (exp(eps_adj) + 2.0);
    filp_prob = 2.0 * q;
    vector<int> values = {0, +1, -1};

    // generalized randomized response --> ns_adj_list
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
                // if ns_adj_list[j][k] does not exist
                if(ns_adj_list[j].count(k) == 0){
                    if(adj_list[j].count(k) == 0) sign_jk = GRR(0, filp_prob, values);
                    else sign_jk = GRR(adj_list[j][k], filp_prob, values);
                    ns_adj_list[j][k] = sign_jk;
                }
                sign_jk = ns_adj_list[j][k];
                if((sign_ij * sign_ik * sign_jk) == 1){
                    num_ns_b_triangles_user[i] += 1;
                } else if((sign_ij * sign_ik * sign_jk) == -1){
                    num_ns_u_triangles_user[i] += 1;
                }
            }
        }
    }
    if(SS_tag == 0){
        // deleted adjacency list after projection
        vector<unordered_map<int, int>> deleted_adj_list(num_nodes);
        // noisy max degree: degree + Lap
        // use the noisy max degree to calculate the global sensitivity
        calculate_noisy_max_degree(degree, eps_degree, max_degree, ns_max_degree);

        // if max_degree exceeds ns_max_degree, then perform graph projection
        if((double)max_degree > ns_max_degree){
            ns_max_degree_floor = (int)floor(ns_max_degree);
            for(i=0; i<num_nodes; ++i){
                // if degree[i] exceeds ns_max_degree, then perform graph projection
                if(degree[i] > ns_max_degree){
                    // randomly generate 0, 1, ..., degree[i]-1 --> rndperm
                    rndperm = vector<int>(degree[i]);
                    generate_random_permutation(rndperm, degree[i], degree[i]);

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
        vector<unordered_map<int, int>>().swap(deleted_adj_list);

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
                } else if(1.0 < s_optimal <= 2.0){
                    SS = std::max(SS, exp(-2.0*beta) * 2.0 * (double)(ldegree_user+1));
                }
            } else{
                SS = std::max(SS, exp(-1.0*beta*(double)(i-ldegree_user)) * 2.0 * (double)(i-1));
            }
            num_ns_b_triangles = (double)num_ns_b_triangles_user[i] - q*(double)num_2_stars_user[i];
            // add the noise for each user
            num_ns_b_triangles +=  stats::rlaplace(0.0, 2.0*SS/eps_triangle, engine);
            num_ns_triangles[0] += num_ns_b_triangles;

            // unbalanced triangles
            num_ns_u_triangles = (double)num_ns_u_triangles_user[i] - q*(double)num_2_stars_user[i];
            num_ns_u_triangles += stats::rlaplace(0.0, 2.0*SS/eps_triangle, engine);
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
    vector<unordered_map<int, int>>().swap(ns_adj_list);
}


int main(int argc, char *argv[]){
    if(argc < 2){
        cout << "Usage: " << argv[0] << " [dataset] ([eps (default: 1.0)] [eps-alloc (default: 2-9-9)]) [mech (default: 3)] [sample nodes (default: -1)] [iters (default: 1)] [repeats (default: 100)])" << endl;
        cout << "[dataset]: dataset name" << endl;
        cout << "[eps]: privacy budget" << endl;
        cout << "[eps-alloc]: privacy budget allocation: degree-adj-triangle" << endl;
        cout << "[mech]: mechanism (0: non-interactive local model (w/o empirical estimation) 1: non-interactive local model (w/ empirical estimation) 2: interactive local model (GS projection), 3: interactive local model (SU))" << endl;
        cout << "[sample nodes]: number of sample nodes (1.0: total nodes)" << endl;
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

    string edge_path = "./data/" + dataset + "_edges.csv";
    int total_num_nodes;
    int num_sample_nodes;
    double delta;
    char *tok;
    vector<vector<int>> node_order;
    vector<unordered_map<int, int>> adj_list;
    // number of balanced and unbalanced triangles
    vector<long long> num_triangles(2);
    vector<double> num_ns_triangles(2);
    double absolute_error_l1 = 0.0, absolute_error_l1_avg = 0.0;
    double relative_error_l1 = 0.0, relative_error_l1_avg = 0.0;

    // Assign the values of command line arguments to variables
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
                    sample_percentage = atoi(argv[5]);
                    if(argc > 6){
                        num_iters = atoi(argv[6]);
                    if(argc > 7) repeats = atoi(argv[7]);
                    }
                }
            }
        }
    }
    // get the total number of nodes --> num_nodes
    total_num_nodes = get_total_num_nodes(edge_path);
    num_sample_nodes = (int)(sample_percentage * (double)total_num_nodes);
    generate_node_order(dataset, num_sample_nodes, total_num_nodes, num_iters, node_order);
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

    cout << "dataset: " << dataset << " total eps: " << total_eps << " " << eps_alloc[0] << "-" << eps_alloc[1] << "-" << eps_alloc[2] << " mech: " << mech << " total_num_nodes: " << total_num_nodes << " num_sample_nodes: " << num_sample_nodes << " num_iters: " << num_iters << " repeats: " << repeats << endl;
    delta = 1.0 / (10.0 * (double)(num_sample_nodes));

    for(int i=0; i<num_iters; ++i){
        // Initialization
        adj_list = vector<unordered_map<int, int>>(num_sample_nodes);
        num_triangles[0] = 0;
        num_triangles[1] = 0;
        // load edges from the edge file --> adj_list
        load_adj_list(edge_path, node_order[i], adj_list);
        calculate_triangles(adj_list, num_triangles);
        for(int repeat=0; repeat<repeats; ++repeat){
            num_ns_triangles[0] = 0.0;
            num_ns_triangles[1] = 0.0;
            if(mech == 0){
                calculate_triangle_non_interactive(adj_list, eps_adj, 0, num_ns_triangles);
            } else if(mech == 1){
                calculate_triangle_non_interactive(adj_list, eps_adj, 1, num_ns_triangles);
            } else if(mech == 2){
                calculate_triangle_interactive(adj_list, 0, eps_degree, eps_adj, eps_triangle, delta, num_ns_triangles);
            } else{
                calculate_triangle_interactive(adj_list, 1, eps_degree, eps_adj, eps_triangle, delta, num_ns_triangles);
            }
            absolute_error_l1 = fabs(num_ns_triangles[0]-(double)num_triangles[0]) + fabs(num_ns_triangles[1]-(double)num_triangles[1]);
            relative_error_l1 = absolute_error_l1 / ((double)num_triangles[0] + (double)num_triangles[1]);
            absolute_error_l1_avg += absolute_error_l1;
            relative_error_l1_avg += relative_error_l1;
        }
        vector<unordered_map<int, int>>().swap(adj_list);
    }
    cout << "Average Absolute Error: " << absolute_error_l1_avg/(double)(num_iters*repeats) << " Average Relative Error: " << relative_error_l1_avg/(double)(num_iters*repeats) << endl;
    cout << "*******************************************************************************" << endl;
    vector<vector<int>>().swap(node_order);
    return 0;
}
