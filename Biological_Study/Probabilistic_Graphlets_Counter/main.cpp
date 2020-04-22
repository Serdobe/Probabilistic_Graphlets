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
Main File

Last Update: 11/02/2012
---------------------------------------------------- */

#include "CmdLine.h"

// "Global" Variables (accessible on every src file)
bool  Global::show_occ;
FILE *Global::occ_file;
bool Global::show_probabilities;
FILE *Global::probabilities_file_Total_Nodes;
FILE *Global::probabilities_file_Orbit_Node;
bool Global::show_graphlet_zero;

//To paralelize:

bool Global::show_paral;
int Global::num_from;
int Global::num_to;

// Main Function
int main(int argc, char **argv) {

  CmdLine::init(argc, argv);
  CmdLine::decide_action();  
  CmdLine::finish();

  return 0;
 }

