#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <queue>
#include <utility>
#include <iostream>

using namespace std;

vector<pair<int,int>> dijkstra(
    unordered_map<int,vector<pair<int,int>>> graph, 
    int source_vertex, 
    const unordered_set<int>& final_vertices
) {
    // Returns a vector of (v,d) :: pair<int,int> where v is 
    // a vertex in final vertices and d is the length of the
    // shortest path from source_vertex to v in the given
    // multigraph. The multigraph is assumed to be represented 
    // by the map graph, a mapping from vertices (int)
    // to vectors of (neigbour, d) pairs.


    // initialize the result vector of travel costs:
    vector<pair<int,int>> travel_costs;

    // initialize the cost dictionary, where we initially
    // set all path lengths to other nodes to "infinity" = 10000
    vector<int> cost(graph.size(), 10000);

    // Create the priority queue
    priority_queue<pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> Q;

    // Pute the source node in the queue:
    Q.emplace(pair<int,int>{0,source_vertex});

    // create a hash set of finished nodes:
    unordered_set<int> finished;

    // Canonical Dijkstra loop
    while (!Q.empty()) {
        int total_cost = Q.top().first;
        int top = Q.top().second;
        Q.pop();
        cost[top] = total_cost;
        
        if (final_vertices.find(top) != final_vertices.end() && finished.find(top) == finished.end())
            // if the current top node is in final_vertices, we should emplace it in the result
            // adjacency vector.
            travel_costs.push_back(pair<int,int>{top, total_cost});

        finished.emplace(top); 
        for(pair<int,int> nbr_dt : graph[top]) {
            int dt = nbr_dt.second;
            int neighbour = nbr_dt.first;
            if (finished.find(neighbour) == finished.end() && cost[neighbour] > dt + total_cost) {
                Q.emplace(pair<int,int>{dt + total_cost, neighbour});
                cost[neighbour] = dt + total_cost;
            }
        }
    }
    return travel_costs;
}

unordered_map<int,vector<pair<int,int>>> reduce (
    const unordered_map<int,vector<pair<int,int>>>& graph, 
    const vector<int>& final_vertices_list
) {

    // Reduce a graph (V,E) to (V',E') where V' are the vertices in
    // final_vertices_list and (v,v', d) in E' if and only if there
    // is a shortest path from v to v' in (V,E) of length d.

    unordered_set<int> final_vertices_set;
    unordered_map<int,vector<pair<int,int>>> reduced_graph;
    for (int vertex : final_vertices_list)
        final_vertices_set.emplace(vertex);
    for (int vertex : final_vertices_list) {
        vector<pair<int,int>> reduced_neigbhours = dijkstra(graph, vertex, final_vertices_set);
        reduced_graph[vertex] = reduced_neigbhours;
    }
    return reduced_graph;
}

int pth_binary_digit(int n, int p) {
    // computes the pth binary digit of positive integer n >= 0
    return (n % (1<<(p+1))) / (1<<p);
}
    
class WaterProblem {
    // Stores |V|, the cardinality of V.
    int v;                                  

    // Stores V as a vector of integers.
    // We number the intersections n = 0, 1, ... , v-1 .
    vector<int> intersections;

    // Stores |W|, the cardinality of W
    int w;

    // Stores W, a list of v in V that are waterpumps
    vector<int> waterpumps_list;

    // For w = waterpumps_list[i], we define 
    // index_in_waterpumps_list[w] = i
    // This means that we can enumerate 
    // the waterpumps 1 to w-1 and also
    // quickly look up the enum value of any waterpump w.
    unordered_map<int,int> index_in_waterpumps_list;

    // This is t from the problem statement (V,W,E,t), 
    // but we capitalize it to distinguish from local use
    // of t in subproblems.
    int T;

    // This stores for every v in V, a list of (v',d) 
    // such that (v,v',d) in E or (v',v,d) in E
    // This means we treat the undirected multigraph as a
    // directed multigraph, but 
    // with v -> v' iff v' -> v exists. This allows 
    // us to quickly look up neighbours in both directions
    // along an edge (v,v',d), which will be needed for 
    // Dijkstra's Algorithm
    unordered_map<int,vector<pair<int,int>>> graph;

public:
    WaterProblem() {
        // The natural number e (number of edges) 
        // as in the assignment:
        int e;                  

        // Read the first line of input into 
        // the described constants:
        cin >> v >> w >> e >> T;

        // read the consecutive w lines "w" and 
        // push these waterpumps into the vector.
        // also store the index in the vector i in 
        // the map index_in_waterpumps_list[w].
        for (int i = 0; i < w; i++) {
            int waterpump;
            cin >> waterpump;
            waterpump -= 1;                          // 0 based numbering
            waterpumps_list.push_back(waterpump);
            index_in_waterpumps_list[waterpump] = i;
        }
        
        // for each v in V, initialize graph[vertex] 
        // as an empty list for
        // elements of type (int, int), where the 
        // first integer will be v' (a neighbour) and the
        // second d, for (v,v',d) or (v',v,d) in E.
        for (int vertex = 0; vertex < v; vertex++) {
            intersections.push_back(vertex);
            graph[vertex] = vector<pair<int,int>>{};
        }


        // read the consecutive v lines "v  v' d" and 
        // store two edges: v gets neighbour (v',d) 
        // and v' gets (v,d)
        for (int _ = 0; _ < e; _++) {
            int v1, v2, d;
            cin >> v1 >> v2 >> d;
            graph[v1-1].push_back(pair<int,int>{v2-1,d}); // 0 based numbering
            graph[v2-1].push_back(pair<int,int>{v1-1,d});
        }
    }

    bool is_elem(int w, int set_number) {
        // returns whether w in W', where set_number = d(W')
        return pth_binary_digit(set_number, index_in_waterpumps_list[w]) == 1;
    }

    bool is_water_pump(int w) {
        // returns whether w is a water pump
        return index_in_waterpumps_list.find(w) != index_in_waterpumps_list.end();
    }

    int add_to_set(int w, int set_number) {
        // returns d(W'+w) given set_number = d(W') and w
        return set_number + (1 << index_in_waterpumps_list[w]);
    }

    bool is_full_set(int set_number) {
        // returns whether W' == W where set_number = d(W')
        return set_number + 1 == (1<<w);
    }
        
    
    void reduce_graph() {
        copy(waterpumps_list.begin(), waterpumps_list.end(), intersections.begin());
        intersections.resize(waterpumps_list.size());
        if (!is_water_pump(0))
            intersections.push_back(0);
        graph = reduce(graph, intersections);    
    }
        
        

    int simplify_solve() {
        // Reduce the graph:
        reduce_graph();

        // Solve the reduced problem"
        return solve();
    }
        

    int solve() {
        // Create the memoization map outside the function max_water_pumed:
        unordered_map<int,unordered_map<int,unordered_map<int,int>>> memo;

        // Initialize a lower bound on the best possible value,
        // which is 0.
        int lower_bound = 0;

        // Dive into the recursion
        return max_water_pumped(T, 0, 0, w, memo, 0, lower_bound);
    }

    int upper_bound(int t, int nr_available_pumps) {
        // returns an upper bound for how much total water we can
        // still pump away given that we are in a state with 
        // timelimit t and nr_available_pumps pumps that are
        // not yet turned on.
    
        const int m = nr_available_pumps; 
            // this is how many pumps can still be turned on in the 
            // remaining time, if we could travel between 
            // intersections in 0 minutes
        return 200 * m * t - 1000 * m * (m + 1); 
    }

        int max_water_pumped(
        int t, 
        int current_node, 
        int set_number,  
            // The three variables (t, current_node, set_number) 
            // completely describe the state.

        int nr_available_pumps,  
            // This number can in principle be calculated 
            // from set_number, but since that adds O(|W|) complexity 
            // to each call to this function, we keep it 
            // as an extra argument.

        
        unordered_map<int,unordered_map<int,unordered_map<int,int>>>& memo, 
            // the memoization map 
            
        int current_water_pumped, 
            // an extra variable to keep track of how much water has 
            // been pumped by parent calls along the currently explored
            // branch.
            
        int& attained_lower_bound 
            // a reference to the value of the best solution 
            // that has been found yet. 
            // Initially, this reference holds 0. 
            // If we are in a leaf call and a better path is found, 
            // the reference is updated to current_water_pumped
        ) {
        
        if (current_water_pumped + upper_bound(t, nr_available_pumps) 
                    <= attained_lower_bound) {
            // We cannot hope to do better.
            // Might as well return 0 because we don't
            // care about what we do next
            // But DON'T store 0 in the memo map! It breaks the
            // invariance that the memo map contains correct values.        
            return 0;
        }

        else if (t <= 10 || is_full_set(set_number)) {
            // If here, we know that we have found a better solution,
            // otherwise control would have gone to the first if-branch.
            
            // And that we also have reached a base case:
            // Therefore, it is a good idea to update the
            // attained_lower_bound to a new optimal value.
            
            attained_lower_bound = current_water_pumped;
            return 0;
        }
        
        else if (memo.find(t) != memo.end() 
        && memo[t].find(current_node) != memo[t].end() 
        && memo[t][current_node].find(set_number) != memo[t][current_node].end()) {
        
            // This means that we have already done the calculation 
            // for the state (t, current_node, set_number)
            
            return memo[t][current_node][set_number];
        }
         
        
        else if (is_water_pump(current_node) 
            && !is_elem(current_node, set_number)) {

            // The value for the state (t, current_node, set_number) is not
            // present in the memoization table

            // If the current_node is a water pump, but not yet turned on,
            // we should turn it on immediately

            // Hence add current_node to W'

            int set_number_plus_current_node = add_to_set(current_node, set_number);
            nr_available_pumps -= 1;
            current_water_pumped += 200 * (t - 10); 
                // add the grand total immediately; 
                // this is fine if we do it consistently.


            memo[t][current_node][set_number] = max_water_pumped(
                t - 10, 
                current_node, 
                set_number_plus_current_node, 
                nr_available_pumps, 
                memo, 
                current_water_pumped, 
                attained_lower_bound) + 200 * (t - 10);
            return memo[t][current_node][set_number];
        }
        
            
        
        else {
            int max_water = 0;

            // for each neighbour w that is a water pump that is 
            // not-yet turned on, such that it can be reached in dt
            // time such that t - dt > 10 (i.e. there is time left to 
            // turn the pump on), recursively compute 
            // max_water_pumped(t - dt, w, set_number)
            // and return the maximum over all such w.

            
            for(pair<int,int> w_dt : graph[current_node]) {
                int w = w_dt.first;
                int dt= w_dt.second;
                if (is_water_pump(w)
                    && !is_elem(w, set_number) 
                    && t - dt > 10)
                    max_water = max(max_water, 
                        max_water_pumped(
                            t - dt, 
                            w, 
                            set_number, 
                            nr_available_pumps, 
                            memo, 
                            current_water_pumped, 
                            attained_lower_bound));
            }
            memo[t][current_node][set_number] = max_water;
            return memo[t][current_node][set_number];
        }
    }
};

int main(void) {
    WaterProblem water_problem;
    int water_pumped = water_problem.simplify_solve();
    cout << water_pumped;
    return 0;
}
