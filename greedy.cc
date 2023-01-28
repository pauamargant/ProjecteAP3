#include <fstream>
#include <iostream>
#include <vector>
#include <climits>
#include <algorithm>
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
----------------
Greedy Heuristic
----------------
The following code is specific to the greedy algorithm.
*/

// Given the current situation of the production line and two car classes returns true if it's better to add
// this_class than best_class
// The criterion used is based on minimizing the added penalization and, in case of draw taking into account the number
// of upgrades/cars left of that kind and also the complexity of the class.
bool comp_sol(const ProductionLinea &LP, const Car_Class &this_class, const Car_Class &best_class, int pen, int min_pen)
{
    return (pen < min_pen or (pen <= min_pen and
                              ((LP.classes_left[this_class.id] * LP.E->classes[this_class.id].n_req > LP.classes_left[best_class.id] * LP.E->classes[best_class.id].n_req or LP.classes_left[this_class.id] > LP.classes_left[best_class.id]) or
                               (LP.classes_left[this_class.id] == LP.classes_left[best_class.id] and this_class.n_req > best_class.n_req))));
}

// Calculates a solution to the problem using a Greedy approach.
// The following criterium is used:
//   - At each position we calculate the car that, if added, adds the less penalization
//     or its penalization is equal to the minimum and other conditions are satisfied
//     (explained in the comp_sol function)
void greedy(problem_data E, ProductionLinea LP)
{
    vector<int> sol(E.C, -1);
    int pos = 0;
    int penalitzacio = 0;
    while (pos < E.C)
    {

        int min_pen = INT_MAX;
        int best_class = 0;
        for (auto classe : E.classes)
        {
            if (LP.classes_left[classe.id] > 0)
            {
                sol[pos] = classe.id;
                int this_pen = LP.update_penalization(pos, sol);
                if (comp_sol(LP, classe, LP.E->classes[best_class], this_pen, min_pen))
                {
                    min_pen = this_pen;
                    best_class = classe.id;
                }
                LP.undo_penalization(pos, sol);
                sol[pos] = -1;
            }
        }
        LP.add_car(pos, best_class, sol);
        penalitzacio += min_pen;
        pos++;
    }
    penalitzacio += LP.final_penalization(pos - 1, sol);
    write_solution(sol, penalitzacio);
}

void find_solution(ifstream &in)
{
    problem_data E;

    in >> E.C >> E.M >> E.K;
    E.upgrades = read_upgrades(in, E.M);
    E.classes = read_classes(in, E.K, E.M);

    vector<int> sol_p(E.C, -1);
    ProductionLinea LP(E);
    greedy(E, LP);
}

int main(int argc, char **argv)
{

    input = argv[1];
    output = argv[2];
    t1 = now();
    ifstream in(input);
    find_solution(in);
}
