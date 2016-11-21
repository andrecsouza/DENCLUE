


/* 
 *  Copyright 2006 Andre Cardoso de Souza
 *  
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  */




#ifndef DENCLUE_H
#define DENCLUE_H


/** INCLUSIONS **/
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <getopt.h>
#include "hypercube.h"
#include "hyperspace.h"
#include "dataset.h"
#include "denclue_functions.h"
using namespace std;



#define MAX_FILENAME 64
#define MAXSIZE_LINE 1024

/** STRUCTS **/

/** Arguments of denclue algorithm.
 * */
typedef struct arguments_struct {

    unsigned int dimension;

    double sigma;  // Influence of an entity in its neighborhood
    double xi;     // Minimum density level for a density-attractor to be significant

    FILE *input_file;  // Stream to the output file
    FILE *output_file; // Stream to the input file

    char input_filename[MAX_FILENAME];  // Name of the output file
    char output_filename[MAX_FILENAME]; // Name of the input file

} arguments_t;




/** METHODS **/


// Main function of DENCLUE
int main( int argc, char **argv);


/** Parse command line arguments and store them in a struct
 *
 *  @param argc argc from main function
 *  @param argv argv from main function
 *  @param arguments struct to store the arguments
 *
 * @return True if all args succesfully parsed. False, otherwise.
 *
 * */
bool parse_args( int argc, char **argv, arguments_t& arguments_t );


/** Print usage of the program.
 *
 * */
void usage();


/** Print clusters to a file.
 *
 *  @param clusters Map of density attractors to entities contained in the corresponding cluster.
 *  @param output_file File to write the clusters.
 *  @param xi Minimum density threshold
 *
 * */
void printOutput( const map< string, vector<DatasetEntity> >& clusters, FILE *output_file , double xi);



#endif

