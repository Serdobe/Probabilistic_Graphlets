/* -------------------------------------------------
      _       _     ___                            
 __ _| |_ _ _(_)___/ __| __ __ _ _ _  _ _  ___ _ _ 
/ _` |  _| '_| / -_)__ \/ _/ _` | ' \| ' \/ -_) '_|
\__, |\__|_| |_\___|___/\__\__,_|_||_|_||_\___|_|  
|___/                                          
    
gtrieScanner: quick discovery of network motifs
Released under Artistic License 2.0
(see README and LICENSE)

Pedro Ribeiro - CRACS & INESC-TEC, DCC/FCUP

----------------------------------------------------
G-Trie Implementation and associated methods

Last Update: 11/02/2012
---------------------------------------------------- */

#ifndef _GTRIE_
#define _GTRIE_

#include "Common.h"
#include "Graph.h"

#define BASE_FORMAT      95
#define BASE_FIRST       ' '
#define BASE_BITS        6

class GraphTree; // forward declaration

class GTrieNode {
 private:
  bool _is_intset_included(list<int> a, list<int> b);
  bool _is_pairset_included(list<iPair> a, list<iPair> b);

 public:

  static int *mymap;
  static bool *used;
  static bool **adjM;
  static int **fastnei;
  static int *numnei;
  static int glk;
  static int numNodes;
  static bool isdir;
  static double *prob;

  //For probabilities:

  static int Subgraph_Degree;

  static vector<vector<float>> Probability_Matrix_GTries;
  static vector<float>Graphlet_zero;
  static vector<float> Graphlet_one;
  static vector<float> Graphlet_two;
  static vector<vector<float>> Expected_Value_2_Nodes_Orbit;
  static vector<vector<float>> Expected_Value_3_Nodes_Orbit;
  static float expected_Graphlet_zero;
  static float expected_Graphlet_one;
  static float expected_Graphlet_two;
  static int Count_Prob_Graph_Zero;

  static float Graphlet_three;
  static float Graphlet_four;
  static float Graphlet_five;
  static float Graphlet_six;
  static float Graphlet_seven;
  static float Graphlet_eight;
  static vector<vector<float>> Expected_Value_4_Nodes_Orbit;
  static float expected_Graphlet_three; //For the printing function.
  static float expected_Graphlet_four;
  static float expected_Graphlet_five;
  static float expected_Graphlet_six;
  static float expected_Graphlet_seven;
  static float expected_Graphlet_eight;


  ///////////////////////////////////////////////


  list< list<int> >   this_node_cond; // This node must be bigger than all these nodes
  list< list<iPair> > cond;           // List of symmetry breaking conditions
  list<GTrieNode *> child;            // List of child g-trie nodes

  bool cond_ok;                       // no need to check for conditions
  bool cond_this_ok;                  // no need to check for this node conditions

  int depth;          // Depth of g-trie node

  bool is_graph;       // Is this node the end of a subGraph?
  int frequency;      // Frequency of this particular subGraph

  bool *in;           // Outward edges
  bool *out;          // Inward edges
  int total_in;       // Number of inward  edges
  int total_out;      // Number of outward edges
  int total_edges;    // Number of inward + outward edges

  int nconn;          // Number of connected nodes
  int *conn;          // Connected nodes


  GTrieNode(int d);   // Create g-trie node with depth 'd'
  ~GTrieNode();

  void insert(Graph *g);
  void insertCond(Graph *g, list<iPair> *cond);
  void showAsText(FILE *f);

  int countNodes();
  int countGraphs();
  int countGraphPaths();
  int maxDepth();

  void zeroFrequency();
  void showFrequency();

  int frequencyGraph(Graph *g);

  void goCondDir();
  void goCondUndir();
  void goCondSample();

  //For probabilities//

  void Function_To_Probabilities(vector<string> &Binary, vector<vector<float>> &Probability_Matrix, int degree);

  ////////////////////

  void insertConditionsFiltered(list<iPair> *cond);

  double countOccurrences();
  int countGraphsApp();

  void populateGraphTree(GraphTree *tree, char *s, int size);
  void populateMap(mapStringInt *m, char *s, int size);

  void writeToFile(FILE *f);
  void readFromFile(FILE *f);

  void Function_Probabilities_Calculaton(vector<string> line_G_tries, Graph *g, int degree_graphlets);

  void cleanConditions();
  void clean(int a, int b);
  void makeConditionsArray();
};

class GTrie {
 private:
  GTrieNode *_root;

 public:
  GTrie();
  ~GTrie();

  void insertGraph(Graph *g);
  void insertGraphCond(Graph *g, list<iPair> *cond);
  void insertGraphString(int size, const char* s);
  void insertGraphNautyString(int size, const char* s, bool dir, int label);

  int frequencyGraphString(int size, const char *s);

  void showAsText(FILE *f);
  double compressionRate();
  int maxDepth();

  void census(Graph *g);
  void censusSample(Graph *g, double *p);

  void showFrequency();

  void unlist();

  int countGraphsApp();
  int countGraphs();
  double countOccurrences();

  void cleanConditions();

  void writeToFile(char *s);
  void readFromFile(char *s);
  void readSubgraphs(int size, bool dir, char *s);

  void populateGraphTree(GraphTree *tree, int size);
  void populateMap(mapStringInt *m, int size);
};


#endif
