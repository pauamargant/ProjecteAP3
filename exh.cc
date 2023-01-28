#include <fstream>
#include <iostream>
#include <vector>
#include <climits>
#include <algorithm>
#include <random>
#include <ctime>
#include <queue>
using namespace std;

/*
Projecte AP3 Curs 2023

Autors: Violeta Sanchez, Pau Amargant
*/

/*
---------
Important
---------
The following code is used in all three algorithms (exh, greedy, mh)
*/

// Global variables
string input, output;
double t1, t2;

// Struct which stores car class information
struct Car_Class
{
    int id;
    int amount;
    vector<int> req_upgrades; // The ith element is 0 if the Class contains the ith upgrade, 0 otherwise
    int n_req = 0;
};

// Struct which stores upgrade (millora) information
struct upgrade
{
    int c;
    int n;
};

/*
Struct used to store information about each input case.
Contains the information which is provided in the input in the adequate structure.
It stays constant throughout the execution of the program.
*/
struct problem_data
{
    int C, M, K;
    vector<upgrade> upgrades;
    vector<Car_Class> classes;
};

/*
Class used to simulate the production line during the execution of the program.
Contains information which changes during the search for an optimal solution.

For each station the capacities vector k contains the number of cars which needed and are inside
its window at each moment of the execution. By calling the update_penalization or remove_penalization
functions the ugprade stations capacities are updated according to the added or removed car.

This class is intended to be used together with a vector which contains the ordered cars.
Cars must be added to the first available position in said vector.

*/
class ProductionLinea
{
public:
    vector<int> k;
    vector<int> classes_left;
    vector<int> upgrades_left;
    problem_data *E;

    ProductionLinea(problem_data &Ent)
        : k(Ent.M, 0), classes_left(Ent.K, 0), upgrades_left(Ent.M, 0), E(&Ent)
    {
        // We calculate the number of classes and upgades required
        for (Car_Class c : E->classes)
        {
            classes_left[c.id] = c.amount;
            for (int i = 0; i < E->M; ++i)
            {
                upgrades_left[i] += c.req_upgrades[i] * c.amount;
            }
        }
    }
    // Adds car with of class Class to position pos in vector sol and updates the penalization
    // It returns the added penalization.
    int add_car(int pos, int Class, vector<int> &sol)
    {
        sol[pos] = Class;
        classes_left[Class]--;
        return update_penalization(pos, sol);
    }
    // Removes car with of class Class from position pos in vector sol and updated the penalizations
    void remove_car(int pos, vector<int> &sol)
    {
        undo_penalization(pos, sol);
        classes_left[sol[pos]]++;
        sol[pos] = -1;
    }

    // Updated the production line when a car is added at position pos
    // Capacities of stations (k) are updated and the value for the penalitzation at position
    // pos i returned.
    int update_penalization(int pos, const vector<int> &sol)
    {
        int added_class = sol[pos];
        vector<int> added_upgrades = E->classes[added_class].req_upgrades;
        int pen = 0;
        for (int id_upgrade = 0; id_upgrade < E->upgrades.size(); ++id_upgrade)
        {
            int n = E->upgrades[id_upgrade].n;
            if (pos >= n)
            {
                int removed_class = sol[pos - n];
                k[id_upgrade] -= E->classes[removed_class].req_upgrades[id_upgrade];
            }
            k[id_upgrade] += added_upgrades[id_upgrade];
            pen += max(k[id_upgrade] - E->upgrades[id_upgrade].c, 0);
        }
        return pen;
    }
    // Reverts changes which are made when a car is inserted at position pos and
    // function actualitzar_pen is called
    void undo_penalization(int pos, const vector<int> &sol)
    {
        int removed_class = sol[pos];
        vector<int> removed_upgrades = E->classes[removed_class].req_upgrades;
        for (int id_upgrade = 0; id_upgrade < E->upgrades.size(); ++id_upgrade)
        {
            k[id_upgrade] -= removed_upgrades[id_upgrade];
            int n_e = E->upgrades[id_upgrade].n;
            // We moved the upgrade window if necessary
            if (pos - n_e >= 0)
            {
                int previous_removed = sol[pos - n_e];
                k[id_upgrade] += E->classes[previous_removed].req_upgrades[id_upgrade];
            }
        }
    }
    // Returns the values of the penalitzations at the end of the production line
    // No changes are made to variables k.
    int final_penalization(int pos, const vector<int> &cotxes)
    {
        int pen = 0;
        for (int id_upgrade = 0; id_upgrade < E->upgrades.size(); ++id_upgrade)
        {
            int left = pos - E->upgrades[id_upgrade].n + 1;
            int k_i = k[id_upgrade];
            while (left < pos)
            {
                int classe_treta = cotxes[left];

                left++;
                k_i -= E->classes[classe_treta].req_upgrades[id_upgrade];
                pen += max(k_i - E->upgrades[id_upgrade].c, 0);
            }
        }
        return pen;
    }
};
// Structure used to contain information related to nods in the search tree. A node is created
// at each position in the solution vector for any car_class which can be added (with penalization under the minimum)
struct node
{
    int id;             // Id of the class added at position pos
    int new_pen;        // Penalization corresponding to the last car added
    ProductionLinea LP; // LiniaProduccio information
};

// Given two solutions (nodes) returns true if n1 is better than n2 according to the criterion used.
struct compNode
{
    bool operator()(const node &n1, const node &n2)
    {
        return not(n1.new_pen < n2.new_pen or (n1.new_pen == n2.new_pen and
                                               ((n1.LP.E->classes[n1.id].n_req * n1.LP.classes_left[n1.id] > n2.LP.E->classes[n2.id].n_req * n2.LP.classes_left[n2.id] or
                                                 n1.LP.classes_left[n1.id] > n2.LP.classes_left[n2.id]) or
                                                (n1.LP.classes_left[n1.id] == n2.LP.classes_left[n2.id] and n1.LP.E->classes[n1.id].n_req > n2.LP.E->classes[n2.id].n_req))));
    }
};
// Used to calculate the execution time
double now()
{
    return clock() / double(CLOCKS_PER_SEC);
}

// Writes to the output file a solution to the problem together with its penalitzacii
void write_solution(const vector<int> &solucio, int penalitzacio)
{

    ofstream out;
    out.open(output);
    out.setf(ios::fixed);
    out.precision(1);

    t2 = now();
    out << penalitzacio << " " << t2 - t1 << endl; // Time spent to find the solution.
    out << solucio[0];
    for (int i = 1; i < solucio.size(); ++i)
    {
        out << " " << solucio[i];
    }
    out << endl;
    out.close();
}
// Reads the millores input and returns it in a vector
vector<upgrade> read_upgrades(ifstream &in, int M)
{

    // For each millora we read capacitat c
    vector<upgrade> upgrades(M);
    for (int i = 0; i < M; ++i)
    {
        in >> upgrades[i].c;
    }
    // For each millora we read the finestra lenth n
    for (int i = 0; i < M; ++i)
    {
        in >> upgrades[i].n;
    }
    return upgrades;
}

// Reads class input and returns a vector of car classes
vector<Car_Class> read_classes(ifstream &in, int K, int M)
{
    vector<Car_Class> classes(K);
    // For each class we read its id and the booleans which determinate
    // which millores are requierd.
    for (int i = 0; i < K; ++i)
    {
        int id;
        in >> id;
        classes[id].id = id;
        in >> classes[id].amount;
        vector<int> m(M);
        for (int j = 0; j < M; ++j)
        {
            in >> m[j];
            classes[id].n_req += m[j];
        }

        classes[id].req_upgrades = m;
    }
    return classes;
}

/*
-----------------
Exhaustive Search
-----------------
The following code is specific to the exhaustive search function
*/

// Returns a priority queue which returns the Nodes in the search tree which must be visited in the following "iteration".
// It only containts nodes in which the penalization is lower than the minimum found.
priority_queue<node, vector<node>, compNode>
get_order(ProductionLinea &LP, int pos, vector<int> sol)
{
    priority_queue<node, vector<node>, compNode> pq;
    for (int id = 0; id < LP.E->K; ++id)
    {
        if (LP.classes_left[id] > 0)
        {
            node n{id, 0, LP};
            n.new_pen = n.LP.add_car(pos, id, sol);
            pq.push(n);
        }
    }
    return pq;
}

// Exhaustive search function
// It looks for an optimal solution by, at each call, oreding the following search nodes
// to be visited based on the added penalization.
void gen(int pos, problem_data &E, ProductionLinea &LP, vector<int> &sol_p, int pen, int &pen_min)
{

    if (pen_min == 0)
        return;

    auto &[C, M, K, mi, classes] = E;

    if (pos == C)
    {

        pen += LP.final_penalization(pos - 1, sol_p);
        if (pen < pen_min)
        {
            pen_min = pen;
            write_solution(sol_p, pen);
        }
    }
    else
    {
        priority_queue<node, vector<node>, compNode> order = get_order(LP, pos, sol_p);
        while (not order.empty())
        {

            node n = order.top();
            order.pop();
            if (n.new_pen + pen < pen_min)
            {
                sol_p[pos] = n.id;
                gen(pos + 1, E, n.LP, sol_p, pen + n.new_pen, pen_min);
            }
        }
    }
}

// Given the problem input through the standard input writes solutions into the output file
void find_solution(ifstream &in)
{
    problem_data E;

    in >> E.C >> E.M >> E.K;
    E.upgrades = read_upgrades(in, E.M);
    E.classes = read_classes(in, E.K, E.M);

    vector<int> sol_p(E.C, -1);
    ProductionLinea LP(E);
    int pen_min = INT_MAX;
    gen(0, E, LP, sol_p, 0, pen_min);
}

int main(int argc, char **argv)
{
    input = argv[1];
    output = argv[2];

    t1 = now();
    ifstream in(input);
    find_solution(in);
}
