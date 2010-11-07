
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







/** INCLUSIONS **/
#include "denclue.h"


/** METHODS **/


// Main function of DENCLUE
int main( int argc, char **argv){


    arguments_t args;

    /* Parse arguments */
    bool parsed_ok = parse_args(argc, argv, args);
    if( !parsed_ok ){
        usage();
        exit(1);
    }


    const unsigned int dimension = args.dimension;
    Dataset dataset(dimension);


    /* Read entities from input and store them */
    while( !feof(args.input_file) ){


        char input_line[MAXSIZE_LINE];

        if( fgets( input_line, MAXSIZE_LINE, args.input_file ) == NULL ) {

            if( feof(args.input_file) ) break;
            perror("Error reading input file");
        }

        const string entity_str( input_line );

        // Create an entity from the read line
        DatasetEntity entity(dimension);
        entity.buildEntityFromString( entity_str );

        // Add the created entity to the dataset
        dataset.addEntity(entity);

    }
    fclose(args.input_file);  // Finish file read



    /* Get lower and upper bounds of the dataset. Shouldn't be emptied. */
    const vector<double>& upper_bounds = dataset.retrieveUpperBound();
    const vector<double>& lower_bounds = dataset.retrieveLowerBound();



    /* Determine hypercubes in the dataset and associate each entity to one of
       them */
    HyperSpace spatial_region( upper_bounds, lower_bounds, args.sigma, args.xi, dimension);
    /*HyperSpace::hypercube_container const *hcubes =*/
    spatial_region.determineSpatialRegions();


    cout << "HyperSpace defined, inserting entities" << endl;

    // Insert entities in the appropriate hypercubes
    Dataset::iterator iter(dataset);
    for( iter.begin() ; !iter.end() ; iter++){

        DatasetEntity ent = dataset.getEntity(*iter);
        spatial_region.insertEntity( ent );
    }


    cout << "Removing low populated hypercubes" << endl;

    /* Determine high populated cubes and remove empty hypercubes or hypercubes
     that are not neighbors of a high populated hypercube */
    spatial_region.removeLowPopulatedHypercubes();

    //DEBUG
    //cout << "Printing hypercubes" << endl;
    /*HyperSpace::hypercube_iterator h_iter = hcubes->begin();
    for( ; h_iter != hcubes->end() ; h_iter++){

        cout << h_iter->second << endl;
    }
    cout << "--------------------" << endl; // */




    cout << "Entities inserted, calculating density functions at each entity"
        << endl;


    /* Calculate density of each entity */
    HyperSpace::EntityIterator hs_iter(spatial_region);

    for( hs_iter.begin() ; !hs_iter.end() ; hs_iter++){

        HyperSpace::EntityIterator calculation_iter(spatial_region);
        calculation_iter.begin();

        double curr_density = DenclueFunctions::calculateDensity( *hs_iter ,
                calculation_iter, args.sigma );

        hs_iter->setDensity( curr_density );
    }


    cout << "Densities calculated, determining density-attractors" << endl;

    /* Determine density attractors and entities attracted by each of them */
    map< string, vector<DatasetEntity> > clusters;  // Map density-attractors to entities

    HyperSpace::EntityIterator iter_entities(spatial_region);
    iter_entities.begin();
    while( !iter_entities.end() ){


        HyperSpace::EntityIterator attractor_entity_iter(spatial_region);
        attractor_entity_iter.begin();

        DatasetEntity curr_attractor =
            DenclueFunctions::getDensityAttractor(*iter_entities,
                    spatial_region, attractor_entity_iter , args.sigma);


        // Ignores density-attractors that don't satisfy minimum density
        // restriction
        if( curr_attractor.getDensity() < args.xi ){

            iter_entities++;
            continue;
        }


        // Create a new cluster if necessary
        if( clusters.count( curr_attractor.getStringRepresentation() ) <= 0 ){

            vector<DatasetEntity> *vec = new vector<DatasetEntity>();
            clusters.insert( make_pair(curr_attractor.getStringRepresentation(), *vec) );
        }


        // Assign current entity to the cluster represented by its density-attractor
        clusters[curr_attractor.getStringRepresentation()].push_back(*iter_entities);

        // Move cursosr to next entity
        iter_entities++;
    }

    cout << "Density attractors determined, determining clusters" << endl;


    /* Merge clusters with a path between them */
    map< string, vector<DatasetEntity> >::iterator outer_iter = clusters.begin();
    while( outer_iter != clusters.end() ){


        // Try to merge a pair of clusters
        map< string, vector<DatasetEntity> >::iterator inner_iter = outer_iter;
        inner_iter++;


        while( inner_iter != clusters.end() ){


            // Build entities that represent each density-attractor
            ostringstream outer_str;
            ostringstream inner_str;
            outer_str << outer_iter->first << Constants::EOL;
            inner_str << inner_iter->first << Constants::EOL;
            DatasetEntity outer(args.dimension);
            DatasetEntity inner(args.dimension);
            outer.buildEntityFromString( outer_str.str() );
            inner.buildEntityFromString( inner_str.str() );


            // Mark ends of desired path as used in path's sequence
            map<string, bool> usedEntities;
            usedEntities[inner.getStringRepresentation()] = true;
            usedEntities[outer.getStringRepresentation()] = true;

            bool canMerge = DenclueFunctions::pathBetweenExists( outer, inner ,
                    spatial_region, args.xi, args.sigma, usedEntities);


            // Merge clusters if there's an appropriate path between their
            // density-attractors
            if( canMerge ){


                DenclueFunctions::AppendVector( outer_iter->second , inner_iter->second );

                clusters.erase( inner_iter++ );  // Erase appended vector and go to next cluster
                continue;
            }


            inner_iter++;

        }


        outer_iter++;
    }





    /* Print clusters representation to output file */

    printOutput( clusters, args.output_file , args.xi);

    cout << "Clusters written to output file " << args.output_filename << endl;


    return 0;

}




/** Parse command line arguments and store them in a struct
 *
 *  @param argc argc from main function
 *  @param argv argv from main function
 *  @param arguments struct to store the arguments
 *
 * @return True if all args succesfully parsed. False, otherwise.
 *
 * */
bool parse_args( int argc, char **argv, arguments_t& arguments ){


    bool parsed_ok = true;
    int curr_flag = 0;


    // Zeroes arguments
    memset((void *)&arguments, 0, sizeof(arguments_t));


    while( (curr_flag = getopt(argc, argv, "hd:s:x:i:o:")) != -1 ){

        switch(curr_flag){

            case 'd':  // number of dimensions in the dataset
                arguments.dimension = (unsigned) atoi(optarg);
                break;


            case 's':  // sigma
                arguments.sigma = atof(optarg);
                break;

            case 'x':  // xi
                arguments.xi = atof(optarg);
                break;

            case 'i': // input file
                memcpy((void *)arguments.input_filename, optarg, strlen(optarg));
                break;

            case 'o': // output file
                memcpy((void *)arguments.output_filename, optarg, strlen(optarg));
                break;

            default:
                parsed_ok = false;

        }


    }



    /* Verify validity of received values */
    if( arguments.dimension == 0 ){
        cerr << "Number of dimensions must be grater than zero" << endl;
        parsed_ok = false;
    }

    if( arguments.sigma == 0 ){
        cerr << "Sigma must be grater than zero" << endl;
        parsed_ok = false;
    }

    if( arguments.xi == 0 ){
        cerr << "Xi must be grater than zero" << endl;
        parsed_ok = false;
    }

    if( strlen(arguments.input_filename) <= 0 ){
        cerr << "Input file name must be defined and must exist" << endl;
        parsed_ok = false;
    }

    if( strlen(arguments.output_filename) <= 0 ){
        cerr << "Output file name must be defined and must exist" << endl;
        parsed_ok = false;
    }


    /* Open files */
    if( parsed_ok ){

        if( (arguments.input_file = fopen( arguments.input_filename, "r" )) == NULL ){
            perror("Error opening input file");
        }


        if( (arguments.output_file = fopen( arguments.output_filename, "w" )) == NULL ){
            perror("Error opening output file");
        }

    }


    return parsed_ok;
}


/** Print usage of the program.
 *
 * */
void usage(){


    cout << "-------------------------------------------" << endl;
    cout << "DENCLUE: density-based clustering algorithm" << endl;
    cout << "Parameters:" << endl;
    cout << "-d\t(number of dimensions of the dataset)" << endl;
    cout << "-s\t(sigma: inlfuence of an entity in its neighborhood)" << endl;
    cout << "-x\t(xi: minimum density level)" << endl;
    cout << "-i\t(input file name)" << endl;
    cout << "-o\t(output file name)" << endl;
    cout << "-h\t(print this help)" << endl;
    cout << "-------------------------------------------" << endl;

    return;
}


/** Print clusters to a file.
 *
 * -> Parameters:
 *  - clusters: Map of density attractors to entities contained in the
 *  corresponding cluster.
 *  - output_file: File to write the clusters.
 *  - xi: Minimum density threshold
 *
 * -> Return: void
 * */
void printOutput( const map< string, vector<DatasetEntity> >& clusters, FILE *output_file, double xi ){


    /* Iterate over the list of clusters, printing the entities that belong to
     * each cluster */
    map< string, vector<DatasetEntity> >::const_iterator iter = clusters.begin();
    unsigned ind_cluster = 0;
    while( iter != clusters.end() ){


        // Avoid printing clusters corresponding to density-attractors that
        // don't satisfy the minimum density restriction
        if( iter->second.empty() ){

            iter++;
            continue;
        }


        ostringstream output_str;
        output_str << "Cluster " << ++ind_cluster << "\tAttractor " << iter->first << endl;

        /* Print each entity of current cluster */
        vector<DatasetEntity>::const_iterator ent_iter = iter->second.begin();
        for( ; ent_iter != iter->second.end() ; ent_iter++){

            output_str << '\t' << *ent_iter << endl;
        }


        // Actually write on output file
        fprintf(output_file, output_str.str().c_str());


        iter++;
    }


    return;
}


