#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <map>
#include <unordered_map>
#include <algorithm>

#include <string.h>
#include <math.h>
#include <time.h>
// #include <unistd.h>

#include "include/stats.hpp"
#include "include/mt19937ar.h"


using namespace std;

// initialization of statslib
stats::rand_engine_t engine(12306);


struct w_t_tuple{
    int w_sum;
    int w_abs;
    int w_min;
    int t;
};

struct w_t_pair{
    int w_sum;
    int t;
    int intersection;
};

struct w_t_triplet{
    int w_abs;
    int w_min;
    int t;
    int intersection;
};


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
    temp = fgets(str, 1024, fp);
    ptr = strchr(str, '#');
    total_num_nodes = atoi(ptr + 1);
    fclose(fp);
    // cout << "total_num_nodes: " << total_num_nodes << endl;
    return total_num_nodes;
}


// randomly generate 0, 1, 2, ..., size-1, and store the first 'num' values into rndperm
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
void load_adj_list(const string &edge_path, vector<int> &node_order, vector<unordered_map<int, int>> &adj_list){
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


void calculate_degree(vector<unordered_map<int,int>> &adj_list, vector<int> &degree){
    int num_nodes = adj_list.size();
    unordered_map<int, int>::iterator iter_1;
    for(int i=0; i<num_nodes; ++i){
        for(iter_1=adj_list[i].begin(); iter_1!=adj_list[i].end(); ++iter_1){
            // degree = positive edges + negative edges
            degree[i] += 1;
        }
    }
}


// calculate the true number of wedges, balanced triangles, and unbalanced triangles
void calculate_wedge_triangle(vector<unordered_map<int, int>> &adj_list, vector<unordered_map<int, pair<int,int>>> &num_wedges, long long &num_wedge_pairs, vector<long long> &num_triangles){
    int i, j, k, sign_ij, sign_ik, sign_jk;
    unordered_map<int, int>::iterator iter_1, iter_2;
    int num_nodes = adj_list.size();
    num_wedge_pairs = 0;
    for(i=0; i<num_nodes; ++i){
        for(iter_1=adj_list[i].begin(); iter_1!=adj_list[i].end(); ++iter_1){
            j = iter_1->first;
            sign_ij = iter_1->second;
            for(iter_2=adj_list[i].begin(); iter_2!=adj_list[i].end(); ++iter_2){
                k = iter_2->first;
                // wedge: j < k
                if(j>=k) continue;
                sign_ik = iter_2->second;
                // calculate the number of wedges
                if(num_wedges[j].count(k) == 0){
                    if((sign_ij * sign_ik == 1)) num_wedges[j][k] = make_pair(1, 0);
                    else num_wedges[j][k] = make_pair(0, 1);
                    ++num_wedge_pairs;
                } else{
                    if((sign_ij * sign_ik == 1)) num_wedges[j][k].first += 1;
                    else num_wedges[j][k].second += 1;
                }
                // triangle: i < j < k
                if(i>=j) continue;
                // calculate the number of triangles
                if(adj_list[j].count(k) > 0){
                    sign_jk = adj_list[j][k];
                    if((sign_ij * sign_ik * sign_jk) == 1) num_triangles[0] += 1;
                    else num_triangles[1] += 1;
                }
            }
        }
    }
}


// calculate local sensitivity at distance s
int calculate_LS_distance_s(vector<unordered_map<int, int>> &adj_list, vector<int> &degree_list, vector<unordered_map<int, pair<int,int>>> &num_wedges, int distance){
    int num_nodes = adj_list.size();
    int LS_s = 0;
    int i, j, w_ij_1, w_ij_2, w_ij_min, a_ij, t_ij, s, LS_s_1, LS_s_2;
    int LS_s_1_upper = num_nodes - 2;
    int LS_s_2_upper = 2 * num_nodes - 4;
    for(i=0; i<num_nodes; ++i){
        for(j=i+1; j<num_nodes; ++j){
            // make s equal distance
            s = distance;
            if(num_wedges[i].count(j) == 0){
                w_ij_1 = 0;
                w_ij_2 = 0;
            } else{
                w_ij_1 = num_wedges[i][j].first;
                w_ij_2 = num_wedges[i][j].second;
            }
            if(adj_list[i].count(j) == 0){
                a_ij = 0;
            } else{
                // a_ij = adj_list[i][j];
                a_ij = 1;
            }
            LS_s_1 = w_ij_1 + w_ij_2;
            LS_s_2 = 2 * abs(w_ij_1 - w_ij_2);
            w_ij_min = std::min(w_ij_1, w_ij_2);
            t_ij = degree_list[i] + degree_list[j] - 2 * (w_ij_1 + w_ij_2 + a_ij);
            // compute LS_s_1
            LS_s_1 += (s + std::min(s, t_ij)) / 2;
            LS_s_1 = std::min(LS_s_1, LS_s_1_upper);
            // compute LS_s_2
            if(a_ij == 0){
                s -= 1;
            }
            if(w_ij_min >= s){
                LS_s_2 += 4 * s;
            } else{
                if(t_ij >= (s-w_ij_min)){
                    LS_s_2 += (2 * w_ij_min + 2 * s);
                } else{
                    LS_s_2 += (2 * w_ij_min + 2 * ((s + t_ij + w_ij_min) / 2));
                    LS_s_2 = std::min(LS_s_2, LS_s_2_upper);
                }
            }
            LS_s = std::max(LS_s, std::max(LS_s_1, LS_s_2));
        }
    }
    // cout << "s = " << s << ", local sensitivity at distance s is " << LS_s << endl;
    return LS_s;
}


// calculate local sensitivity
int calculate_LS(vector<unordered_map<int,int>> &adj_list, vector<unordered_map<int, pair<int,int>>> &num_wedges){
    int num_nodes = adj_list.size();
    int LS = 0;
    int i, j, w_pos, w_neg;
    unordered_map<int,pair<int,int>>::iterator iter;
    for(i=0; i<num_nodes; ++i){
        for(iter=num_wedges[i].begin(); iter != num_wedges[i].end(); ++iter){
            j = iter->first;
            w_pos = iter->second.first;
            w_neg = iter->second.second;
            LS = std::max(LS, w_pos+w_neg);
            if(adj_list[i].count(j) != 0){
                LS = std::max(LS, 2 * abs(w_pos-w_neg));
            }
        }
    }
    // cout << "Local sensitivity: " << LS << endl;
    return LS;
}


// calculate smooth sensitivity: baseline method
double calculate_SS_baseline(vector<unordered_map<int,int>> &adj_list, vector<int> &degree_list, vector<unordered_map<int,pair<int,int>>> &num_wedges, double beta){
    int num_nodes = adj_list.size();
    // compute the local sensitivity
    int LS = calculate_LS(adj_list, num_wedges);
    // cout << "LS: " << LS << endl;
    // when s = 0, smooth sensitivity = exp(-beta * 0) * LS
    double SS = (double)LS;
    int LS_s;
    int s_max = 2 * (num_nodes - 2) + 1;
    for(int s=1; s<s_max; ++s){
        LS_s = calculate_LS_distance_s(adj_list, degree_list, num_wedges, s);
        SS = std::max(SS, exp(-1.0 * (double)s * beta) * (double)LS_s);
    }
    return SS;
}


void calculate_LS_upper_bound(vector<unordered_map<int,int>> &adj_list, vector<unordered_map<int, pair<int,int>>> &num_wedges, int &LS_1, int &LS_2){
    int num_nodes = adj_list.size();
    int w_pos, w_neg;
    unordered_map<int,pair<int,int>>::iterator iter;
    for(int i=0; i<num_nodes; ++i){
        for(iter=num_wedges[i].begin(); iter != num_wedges[i].end(); ++iter){
            w_pos = iter->second.first;
            w_neg = iter->second.second;
            LS_1 = std::max(LS_1, w_pos+w_neg);
            LS_2 = std::max(LS_2, 2 * abs(w_pos-w_neg));
        }
    }
}


// calculate smooth upper bound on local sensitivity
double calculate_SU(vector<unordered_map<int,int>> &adj_list, int LS_1, int LS_2, double beta){
    int num_nodes = adj_list.size();
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


bool customCompare_1(const w_t_tuple &a, const w_t_tuple &b){
    return a.w_sum > b.w_sum;
}

bool customCompare_2(const w_t_tuple &a, const w_t_tuple &b){
    return a.w_abs > b.w_abs;
}

bool customCompare_3(const w_t_tuple &a, const w_t_tuple &b){
    return (a.w_abs + a.w_min) > (b.w_abs + b.w_min);
}


// calculate smooth sensitivity
double calculate_SS(vector<unordered_map<int,int>> &adj_list, vector<int> &degree_list, vector<unordered_map<int,pair<int,int>>> &num_wedges, long long num_wedge_pairs, double beta){
    int num_nodes = adj_list.size();
    double SS = 0.0;
    vector<w_t_tuple> w_t_vec(num_wedge_pairs+1);
    w_t_tuple w_t_ij;
    int i, j, w_ij_1, w_ij_2, a_ij;
    long long index = 0;
    int t_max = 0;
    for(i=0; i<num_nodes; ++i){
        for(j=i+1; j<num_nodes; ++j){
            if(adj_list[i].count(j) == 0) a_ij = 0;
            else a_ij = 1;
            if(num_wedges[i].count(j) == 0){
                t_max = std::max(t_max, degree_list[i]+degree_list[j]-2*a_ij);
            } else{
                w_ij_1 = num_wedges[i][j].first;
                w_ij_2 = num_wedges[i][j].second;
                w_t_ij.w_sum = w_ij_1 + w_ij_2;
                w_t_ij.w_abs = abs(w_ij_1 - w_ij_2);
                w_t_ij.w_min = std::min(w_ij_1, w_ij_2);
                w_t_ij.t = degree_list[i] + degree_list[j] - 2 * (w_ij_1 + w_ij_2 + a_ij);
                w_t_vec[index] = w_t_ij;
                ++index;
            }
        }
    }

    // sort the tuple according to "w_sum" and calculate the LS_s_1
    sort(w_t_vec.begin(), w_t_vec.end(), customCompare_1);

    w_t_vec[index] = {.w_sum=0, .w_abs=0, .w_min=0, .t=t_max};

    int w_sum_0, t_0, w_sum_1, t_1, intersection;
    vector<w_t_pair> surs_pair;
    w_t_pair prev_sur_pair, sur_pair;
    w_sum_0 = w_t_vec[0].w_sum;
    t_0 = w_t_vec[0].t;
    index = 1;
    while(index <= num_wedge_pairs && w_t_vec[index].w_sum == w_sum_0){
        t_0 = std::max(t_0, w_t_vec[index].t);
        ++index;
    }
    prev_sur_pair = {.w_sum=w_sum_0, .t=t_0, .intersection=2*num_nodes-3};
    surs_pair.emplace_back(prev_sur_pair);

    while(index <= num_wedge_pairs){
        w_sum_1 = w_t_vec[index].w_sum;
        t_1 = w_t_vec[index].t;
        ++index;
        while(index <= num_wedge_pairs && w_t_vec[index].w_sum == w_sum_1){
            t_1 = std::max(t_1, w_t_vec[index].t);
            ++index;
        }
        w_sum_0 = prev_sur_pair.w_sum;
        t_0 = prev_sur_pair.t;
        intersection = 2 * (w_sum_0 - w_sum_1) + t_0;
        if(intersection < t_1){
            surs_pair.back().intersection = intersection;
            sur_pair = {.w_sum=w_sum_1, .t=t_1, .intersection=2*num_nodes-3};
            surs_pair.emplace_back(sur_pair);
            prev_sur_pair = sur_pair;
        }
    }

    int LS_s_1;
    int s_max = 2 * (num_nodes - 2 - surs_pair.back().w_sum) - surs_pair.back().t + 1;
    auto iter = surs_pair.begin();
    for(int s=0; s<s_max; ++s){
        if(s >= iter->intersection){
            ++iter;
        }
        LS_s_1 = iter->w_sum + (s + std::min(s, iter->t)) / 2;
        LS_s_1 = std::min(LS_s_1, num_nodes - 2);
        SS = std::max(SS, exp(-1.0 * beta * (double)s) * (double)LS_s_1);
    }
    vector<w_t_pair>().swap(surs_pair);
    // sort the tuple according to "w_abs" and calculate the LS_s_2
    sort(w_t_vec.begin(), w_t_vec.end(), customCompare_2);
    vector<w_t_triplet> surs_triplet_1, surs_triplet_2, surs_triplet_3;
    w_t_triplet prev_sur_triplet, sur_triplet;
    int w_abs_0, w_abs_1, w_min_0, w_min_1;
    w_abs_0 = w_t_vec[0].w_abs;
    w_min_0 = w_t_vec[0].w_min;
    t_0 = w_t_vec[0].t;
    index = 1;
    while(index <= num_wedge_pairs && w_t_vec[index].w_abs == w_abs_0){
        if(w_t_vec[index].w_min > w_min_0){
            w_min_0 = w_t_vec[index].w_min;
            t_0 = w_t_vec[index].t;
        }
        ++index;
    }
    prev_sur_triplet = {.w_abs=w_abs_0, .w_min=w_min_0, .t=t_0, .intersection=2*num_nodes-3};
    surs_triplet_1.emplace_back(prev_sur_triplet);

    while(index <= num_wedge_pairs){
        w_abs_1 = w_t_vec[index].w_abs;
        w_min_1 = w_t_vec[index].w_min;
        t_1 = w_t_vec[index].t;
        ++index;
        while(index <= num_wedge_pairs && w_t_vec[index].w_abs == w_abs_1){
            if(w_t_vec[index].w_min > w_min_1){
                w_min_1 = w_t_vec[index].w_min;
                t_1 = w_t_vec[index].t;
            }
            ++index;
        }
        w_abs_0 = prev_sur_triplet.w_abs;
        w_min_0 = prev_sur_triplet.w_min;
        t_0 = prev_sur_triplet.t;
        intersection = w_abs_0 - w_abs_1 + w_min_0;
        if(intersection < w_min_1){
            surs_triplet_1.back().intersection = intersection;
            sur_triplet = {.w_abs=w_abs_1, .w_min=w_min_1, .t=t_1, .intersection=2*num_nodes-3};
            surs_triplet_1.emplace_back(sur_triplet);
            prev_sur_triplet = sur_triplet;
        }
    }

    w_abs_0 = w_t_vec[0].w_abs;
    w_min_0 = w_t_vec[0].w_min;
    t_0 = w_t_vec[0].t;
    index = 1;
    while(index <= num_wedge_pairs && w_t_vec[index].w_abs == w_abs_0){
        if(3*w_t_vec[index].w_min+w_t_vec[index].t > 3*w_min_0+t_0){
            w_min_0 = w_t_vec[index].w_min;
            t_0 = w_t_vec[index].t;
        }
        ++index;
    }
    prev_sur_triplet = {.w_abs=w_abs_0, .w_min=w_min_0, .t=t_0, .intersection=6*num_nodes-9};
    surs_triplet_2.emplace_back(prev_sur_triplet);

    while(index <= num_wedge_pairs){
        w_abs_1 = w_t_vec[index].w_abs;
        w_min_1 = w_t_vec[index].w_min;
        t_1 = w_t_vec[index].t;
        ++index;
        while(index <= num_wedge_pairs && w_t_vec[index].w_abs == w_abs_1){
            if(3*w_t_vec[index].w_min+w_t_vec[index].t > 3*w_min_1+t_1){
                w_min_1 = w_t_vec[index].w_min;
                t_1 = w_t_vec[index].t;
            }
            ++index;
        }
        w_abs_0 = prev_sur_triplet.w_abs;
        w_min_0 = prev_sur_triplet.w_min;
        t_0 = prev_sur_triplet.t;
        intersection = 2*(w_abs_0-w_abs_1) + 3*w_min_0 + t_0;
        if(intersection < 3*w_min_1){
            surs_triplet_2.back().intersection = intersection;
            sur_triplet = {.w_abs=w_abs_1, .w_min=w_min_1, .t=t_1, .intersection=6*num_nodes-9};
            surs_triplet_2.emplace_back(sur_triplet);
            prev_sur_triplet = sur_triplet;
        }
    }

    // sort the tuple according to "w_abs + w_min" and calculate the LS_s_2
    sort(w_t_vec.begin(), w_t_vec.end(), customCompare_3);
    w_abs_0 = w_t_vec[0].w_abs;
    w_min_0 = w_t_vec[0].w_min;
    t_0 = w_t_vec[0].t;
    index = 1;
    while(index<=num_wedge_pairs && (w_t_vec[index].w_abs+w_t_vec[index].w_min)==(w_abs_0+w_min_0)){
        if(w_t_vec[index].w_min+w_t_vec[index].t > w_min_0+t_0){
            w_abs_0 = w_t_vec[index].w_abs;
            w_min_0 = w_t_vec[index].w_min;
            t_0 = w_t_vec[index].t;
        }
        ++index;
    }
    prev_sur_triplet = {.w_abs=w_abs_0, .w_min=w_min_0, .t=t_0, .intersection=2*num_nodes-3};
    surs_triplet_3.emplace_back(prev_sur_triplet);

    while(index <= num_wedge_pairs){
        w_abs_1 = w_t_vec[index].w_abs;
        w_min_1 = w_t_vec[index].w_min;
        t_1 = w_t_vec[index].t;
        ++index;
        while(index<=num_wedge_pairs && (w_t_vec[index].w_abs+w_t_vec[index].w_min)==(w_abs_1+w_min_1)){
            if(w_t_vec[index].w_min+w_t_vec[index].t > w_min_1+t_1){
                w_abs_1 = w_t_vec[index].w_abs;
                w_min_1 = w_t_vec[index].w_min;
                t_1 = w_t_vec[index].t;
            }
            ++index;
        }
        w_abs_0 = prev_sur_triplet.w_abs;
        w_min_0 = prev_sur_triplet.w_min;
        t_0 = prev_sur_triplet.t;
        intersection = 2*(w_abs_0-w_abs_1) + 3*w_min_0 - 2*w_min_1 + t_0;
        if(intersection < w_min_1+t_1){
            surs_triplet_3.back().intersection = intersection;
            sur_triplet = {.w_abs=w_abs_1, .w_min=w_min_1, .t=t_1, .intersection=2*num_nodes-3};
            surs_triplet_3.emplace_back(sur_triplet);
            prev_sur_triplet =sur_triplet;
        }
    }

    int LS_s_2_1, LS_s_2_2, LS_s_2_3;
    s_max = 2*(num_nodes-2-surs_triplet_1.back().w_abs-2*surs_triplet_1.back().w_min)+surs_triplet_1.back().w_min-surs_triplet_1.back().t+1;
    s_max = std::max(s_max, 2*(num_nodes-2-surs_triplet_2.back().w_abs-2*surs_triplet_2.back().w_min)+surs_triplet_2.back().w_min-surs_triplet_2.back().t+1);
    s_max = std::max(s_max, 2*(num_nodes-2-surs_triplet_3.back().w_abs-2*surs_triplet_3.back().w_min)+surs_triplet_3.back().w_min-surs_triplet_3.back().t+1);

    auto iter_1 = surs_triplet_1.begin();
    auto iter_2 = surs_triplet_2.begin();
    auto iter_3 = surs_triplet_3.begin();
    for(int s=0; s<s_max; ++s){
        if(s >= iter_1->intersection){
            ++iter_1;
        }
        LS_s_2_1 = 2*iter_1->w_abs + 2*std::min(iter_1->w_min,s) + 2*((s+std::min(s,iter_1->w_min+iter_1->t))/2);

        if(3*s >= iter_2->intersection){
            ++iter_2;
        }
        LS_s_2_2 = 2*iter_2->w_abs + 2*std::min(iter_2->w_min,s) + 2*((s+std::min(s,iter_2->w_min+iter_2->t))/2);

        if(s >= iter_3->intersection){
            ++iter_3;
        }
        LS_s_2_3 = 2*iter_3->w_abs + 2*std::min(iter_3->w_min,s) + 2*((s+std::min(s,iter_3->w_min+iter_3->t))/2);

        LS_s_2_3 = std::max(std::max(LS_s_2_1, LS_s_2_2), LS_s_2_3);
        LS_s_2_3 = std::min(LS_s_2_3, 2*num_nodes-4);
        SS = std::max(SS, exp(-1.0 * beta * (double)s) * (double)LS_s_2_3);
    }
    vector<w_t_triplet>().swap(surs_triplet_3);
    vector<w_t_triplet>().swap(surs_triplet_2);
    vector<w_t_triplet>().swap(surs_triplet_1);
    vector<w_t_tuple>().swap(w_t_vec);
    return SS;
}


int main(int argc, char *argv[]){
    if(argc < 2){
        cout << "Usage: " << argv[0] << " [dataset] ([eps (default: 0.1)] [mech (default: 2)] [sample nodes (default: 1.0)] [iters (default: 1)] [repeats (default: 100)])" << endl;
        cout << "[dataset]: dataset name" << endl;
        cout << "[eps]: privacy budget" << endl;
        cout << "[mech]: mechanism (0: GS 1: SS (baseline) 2: SS 3: SU" << endl;
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
    int mech = 2;
    double sample_percentage = 1.0;
    int num_iters = 1;
    int repeats = 100;

    string edge_path = "./data/" + dataset + "_edges.csv";
    int total_num_nodes;
    int num_sample_nodes;
    double delta;
    vector<vector<int>> node_order;
    vector<unordered_map<int, int>> adj_list;
    vector<int> degree;
    //number of wedges: num_wedges[i][j].first: ++ and --, num_wedges[i][j].second: +-
    vector<unordered_map<int, pair<int,int>>> num_wedges;
    long long num_wedge_pairs;
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
    total_num_nodes = get_total_num_nodes(edge_path);
    num_sample_nodes = (int)(sample_percentage * (double)total_num_nodes);
    generate_node_order(dataset, num_sample_nodes, total_num_nodes, num_iters, node_order);
    cout << "dataset: " << dataset << " eps: " << eps << " mech: " << mech << " total_num_nodes: " << total_num_nodes << " sample_percentage: " << sample_percentage << " num_iters: " << num_iters << " repeats: " << repeats << endl;
    delta = 1.0 / (10.0 * (double)(num_sample_nodes) * (double)(num_sample_nodes - 1) / 2.0);

    for(int i=0; i<num_iters; ++i){
        // initialization
        adj_list = vector<unordered_map<int, int>>(num_sample_nodes);
        degree = vector<int>(num_sample_nodes, 0);
        num_wedges = vector<unordered_map<int, pair<int, int>>>(num_sample_nodes);
        num_triangles[0] = 0;
        num_triangles[1] = 0;
        // load edges from the edge file --> adj_list and degree
        load_adj_list(edge_path, node_order[i], adj_list);
        calculate_degree(adj_list, degree);
        // calculate the true number of wedges and triangles
        calculate_wedge_triangle(adj_list, num_wedges, num_wedge_pairs, num_triangles);
        // calculate the noisy number of balanced and unbalanced triangles
        if(mech == 0){
            double GS;
            GS = double(2 * num_sample_nodes - 4);
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
            double SS = 0.0;
            double beta = eps / (4.0*(2.0+log(2.0/delta)));
            // judge: SS = LS_uppper_bound or SS != LS_upper_bound
            int LS_1=0, LS_2=0;
            calculate_LS_upper_bound(adj_list, num_wedges, LS_1, LS_2);
            double beta_inverse = 1.0 / beta;
            // cout << "LS_1: " << LS_1 << " LS_2: " << LS_2 << " beta_inverse: " << beta_inverse << endl;
            if(beta_inverse <= double(LS_1) && beta_inverse <= double(LS_2)/4.0){
                SS = std::max((double)LS_1, (double)LS_2);
            }
            if(mech == 1){
                SS = calculate_SS_baseline(adj_list, degree, num_wedges, beta);
            } else if(mech == 2){
                SS = calculate_SS(adj_list, degree, num_wedges, num_wedge_pairs, beta);
            } else if(mech == 3){
                calculate_LS_upper_bound(adj_list, num_wedges, LS_1, LS_2);
                SS = calculate_SU(adj_list, LS_1, LS_2, beta);
            }
            for(int repeat=0; repeat<repeats; ++repeat){
                num_ns_triangles[0] = num_triangles[0] + stats::rlaplace(0.0, 2.0*SS/eps, engine);
                num_ns_triangles[1] = num_triangles[1] + stats::rlaplace(0.0, 2.0*SS/eps, engine);
                num_ns_triangles[0] = std::max(0.0, num_ns_triangles[0]);
                num_ns_triangles[1] = std::max(0.0, num_ns_triangles[1]);
                absolute_error_l1 = fabs(num_ns_triangles[0]-(double)num_triangles[0]) + fabs(num_ns_triangles[1]-(double)num_triangles[1]);
                relative_error_l1 = absolute_error_l1 / ((double)num_triangles[0] + (double)num_triangles[1]);
                absolute_error_l1_avg += absolute_error_l1;
                relative_error_l1_avg += relative_error_l1;
            }
        }
        vector<unordered_map<int, pair<int, int>>>().swap(num_wedges);
        vector<int>().swap(degree);
        vector<unordered_map<int, int>>().swap(adj_list);
    }
    cout << "Average Absolute Error: " << absolute_error_l1_avg/(double)(num_iters*repeats) << " Average Relative Error: " << relative_error_l1_avg/(double)(num_iters*repeats) << endl;
    cout << "*******************************************************************************" << endl;
    vector<vector<int>>().swap(node_order);
    return 0;
}
