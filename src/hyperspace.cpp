


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





/* INCLUSIONS */
#include "hyperspace.h"


/* METHODS */


// Constructor
HyperSpace::HyperSpace( const vector<double>& up_bound , const vector<double>& lw_bound, double sigma , double xi, const unsigned num_dimensions) : dimension(num_dimensions), sigma(sigma), xi(xi) {


    /* Allocate storage for spatial bounds */
    this->lower_bounds = new double[num_dimensions];
    this->upper_bounds = new double[num_dimensions];


    double edge_length = this->hypercubeEdgeLenght();
    for(unsigned i=0 ; i < num_dimensions ; i++){

        this->lower_bounds[i] = lw_bound[i];
        this->upper_bounds[i] = edge_length * ceill(up_bound[i] / edge_length ); // Ensure that all regions will have the same size
    }

}


// Copy-constructor
HyperSpace::HyperSpace( const HyperSpace& other ) : dimension(other.dimension), sigma(other.sigma), xi(other.xi) {


    /* Allocate storage for spatial bounds */
    this->lower_bounds = new double[this->dimension];
    this->upper_bounds = new double[this->dimension];


    // Copy bounds
    for(unsigned i=0 ; i < this->dimension ; i++){

        this->lower_bounds[i] = other.lower_bounds[i];
        this->upper_bounds[i] = other.upper_bounds[i];
    }


    // Copy hypercubes
    this->hypercubes.clear();
    HyperSpace::hypercube_iterator it = other.hypercubes.begin();
    for( ; it != other.hypercubes.end() ; it++){
        this->hypercubes.insert(*it);
    }


}

/** Determine the regions of the space based on the parameter sigma.
 * The index string is formed using the appropriate function.
 *
 * @return a mapping of the multi-dimensional index string and the
 *  hypercubes.
 *
 * */
const map< string, HyperCube >* HyperSpace::determineSpatialRegions(){


    unsigned num_partitions[this->dimension]; // Number of partitions of each spatial component
    unsigned long num_hypercubes = 1;   // Number of hypercubes to create


    /* Calculate the number of partitions in each spatial component and
     * the total number of hypercubes to create. */
    double edge_length = this->hypercubeEdgeLenght();
    for( unsigned i=0; i < this->dimension ; i++){

        num_partitions[i] = (unsigned) ceil((this->upper_bounds[i] -
                    this->lower_bounds[i]) / (edge_length));

        num_hypercubes = num_hypercubes * num_partitions[i];
    }



    /* Instantiate all hypercubes */

    // Calculate upper bounds of hypercube at spatial origin
    double initial_upp_bounds[this->dimension];
    for(unsigned i=0 ; i < this->dimension ; i++){

        initial_upp_bounds[i] = this->lower_bounds[i] + this->hypercubeEdgeLenght();
    }
    string initial_upp_bounds_str = HyperCube::getKeyFromArray( initial_upp_bounds , this->dimension, this->hypercubeEdgeLenght() );

    // Execute the instantiation of hypercubes
    this->createHyperCubes( num_partitions , initial_upp_bounds_str, this->hypercubeEdgeLenght() );


    return &(this->hypercubes);

}


/** Instantiate all hypercubes in the HyperSpace. A hypercube and all
 * its neighbors are instantiated in a recursive manner.
 *  Neighborhood is determined by generating all hypercubes that have at most
 * (d-1) bound components with a difference of 1 sigma between their bound and
 * current bound.
 *
 *
 *  @param num_partitions Array with the number of partitions of each spatial component
 *  @param upp_bounds_str string representation of upper bounds of the initial hypercube.
 *  @param edge_length lower bounds of the initial hypercube.
 *
 *
 * */
void HyperSpace::createHyperCubes( const unsigned *num_partitions , const string& upp_bounds_str  , double edge_length){



    /** Instantiate a hypercube **/
    const string& curr_cube_key = upp_bounds_str;

    // Get bounds of hypercube
    double *upp_bounds = HyperCube::getArrayFromKey( curr_cube_key , this->dimension);
    double low_bounds[this->dimension];
    for(unsigned i=0 ; i < this->dimension ; i++){
        low_bounds[i] = upp_bounds[i] - this->hypercubeEdgeLenght();
    }


    // Build a new hypercube object
    HyperCube* cube = new HyperCube( this->dimension, upp_bounds, this->hypercubeEdgeLenght() );
    HyperCube& curr_cube = *cube;


    /** Determine neighbors of instantiated hypercube **/
    vector<string> neighbors_keys;  // Keys of all neighbors of the cube created above
    unsigned num_neighbors = (unsigned) pow(3*1.0 , this->dimension * 1.0);  // For simplicity, include curr cube

    /* Set first neighbor */
    int neighbor[this->dimension];
    for(unsigned i=0 ; i < this->dimension; i++){

        neighbor[i] = -1;
    }


    /* Generate neighbors. First neighbor is already generated */
    for(unsigned i=0 ; i < num_neighbors; i++){

        bool outside_space = false; // Indicate that a neighbor is out of bounds

        /* Calculate neighbor bounds */
        double next_upp_bounds[this->dimension];
        double next_low_bounds[this->dimension];

        for(unsigned dimension=0; dimension < this->dimension ; dimension++){

            double increment = this->hypercubeEdgeLenght();


            // Update upper bound, verifying its validity
            if( (upp_bounds[dimension] == this->upper_bounds[dimension]) && (neighbor[dimension] > 0) ){
                outside_space = true;
                break;
            }
            else if( (upp_bounds[dimension] == this->lower_bounds[dimension]) && (neighbor[dimension] < 0) ){
                outside_space = true;
                break;
            }
            next_upp_bounds[dimension] = upp_bounds[dimension] + ( neighbor[dimension] * increment );
            //outside_space = outside_space || (next_upp_bounds[dimension] <
            //        (this->lower_bounds[dimension] + increment)) ||
            //    (next_upp_bounds[dimension] > this->upper_bounds[dimension]);


            // Update lower bound, verifying its validity
            if( (low_bounds[dimension] == this->lower_bounds[dimension]) && (neighbor[dimension] < 0) ){
                outside_space = true;
                break;
            }
            else if( (low_bounds[dimension] == this->upper_bounds[dimension]) && (neighbor[dimension] > 0) ){
                outside_space = true;
                break;
            }
            next_low_bounds[dimension] = low_bounds[dimension] + ( neighbor[dimension] * increment );
            //outside_space = outside_space || (next_low_bounds[dimension] >
            //        (this->upper_bounds[dimension] - increment)) ||
            //    (next_low_bounds[dimension] < this->lower_bounds[dimension]);


        }



        // Store neighbor key if neighbor is a valid hypercube
        if( !outside_space ){

            string neighbor_str_key = HyperCube::getKeyFromArray(
                    next_upp_bounds, this->dimension,
                    this->hypercubeEdgeLenght() );

            if( neighbor_str_key != curr_cube_key ){  // Avoid inserting itself as neighbor

                neighbors_keys.push_back(neighbor_str_key);

            }
        }


        /* Generate next neighbor */
        for( unsigned dimension = this->dimension ; dimension > 0 ; dimension--){

            int index = dimension - 1;
            if( neighbor[index] == -1 ){  neighbor[index] =  0;    break;  }
            if( neighbor[index] ==  0 ){  neighbor[index] =  1;    break;  }
            if( neighbor[index] ==  1 ){  neighbor[index] = -1;    }
        }


    }  // End of neighborhood generation



    /* Assign neighbors to current hypercube and add it to the container */
    curr_cube.setNeighbors( neighbors_keys );

    pair<string, HyperCube> *element = new pair<string,HyperCube>(curr_cube_key, curr_cube);
    this->hypercubes.insert( *element ); // Add current cube to the list of hypercubes
    delete cube;
    delete element;


    /* For each neighbor that doesn't exist yet, invokes this creation method
     * recursively */
    vector<string>::const_iterator it = neighbors_keys.begin();
    for( ; it != neighbors_keys.end() ; it++){


        // Avoid instantiating hypercubes already instantiated
        if( this->hypercubes.count(*it) > 0 )  continue;
        else{


            /* Recursive invocation of this method */
            /*double *neighbor_upp_bounds = HyperCube::getArrayFromKey( *it ,
              this->dimension );

            // Calculate lower bounds of current neighbor
            double neighbor_low_bounds[this->dimension];
            bool outside_space = false;
            for(unsigned i=0 ; i < this->dimension ; i++){

            neighbor_low_bounds[i] = neighbor_upp_bounds[i] - this->hypercubeEdgeLenght();

            outside_space = neighbor_low_bounds[i] < this->lower_bounds[i];
            }
            if( outside_space ) continue;*/


            // Create current neighbor and its neighbors
            HyperSpace::createHyperCubes( num_partitions , *it , this->hypercubeEdgeLenght() );

            //delete[] neighbor_upp_bounds;

        }


    }


    delete[] upp_bounds;

}



/** Insert a dataset entity in the space.
 *
 *  @param entity The entity to insert.
 *
 * */
void HyperSpace::insertEntity( const DatasetEntity& entity ){



    /* Determine the hypercube that should contain the entity */
    double attr_values[this->dimension];
    for(unsigned i=0 ; i < this->dimension ; i++){

        double sigma = this->hypercubeEdgeLenght();
        attr_values[i] = sigma * (floor( entity.getComponentValue(i) / sigma) + 1);

    }


    /* Insert the entity in the corresponding hypercube */

    // Transform key in string
    string key = HyperCube::getKeyFromArray( attr_values, this->dimension , this->hypercubeEdgeLenght());


    if( this->hypercubes.count(key) > 0 ){  // Verify whether wanted hypercube exists
        this->hypercubes.find(key)->second.addObject(entity);
    }
    else {

        cerr << "[insertEntity] HyperCube for key (" << key << ") not found"<< endl;
    }


    return;
}


/** Remove low populated hypercubes, except those who are connected to
 * a high populated hypercube.
 *
 * */
void HyperSpace::removeLowPopulatedHypercubes(){


    /** Remove empty hypercubes **/

    vector<string> deleted_keys;
    this->high_populated_keys.clear();

    map< string, HyperCube>::iterator it = this->hypercubes.begin();
    while( it != this->hypercubes.end() ){


        // Mark high populated cube and store it internally
        if( it->second.numObjects() >= this->minimumObjectsInHypercubes() ){

            this->high_populated_keys.push_back( it->first ); // Store cube key
        }


        // Mark empty cubes
        if( it->second.isEmpty() ){

            deleted_keys.push_back(it->first); // Mark key as deleted
            this->hypercubes.erase(it++);  // Delete current neighbor and forward iterator
        }
        else ++it; // Forward iterator


    }


    /** Remove keys of empty hypercubes from other hypercubes **/
    it = this->hypercubes.begin();
    for( ; it != this->hypercubes.end() ; it++){

        it->second.removeEmptyNeighbors( deleted_keys );
    }


    /** Remove hypercubes that aren't connected to high populated hypercubes
     * */
    it = this->hypercubes.begin();
    while( it != this->hypercubes.end() ){


        if( !(it->second.isNeighbor( high_populated_keys, this->hypercubes )) ){  // Cube isn't connected to high populated cube
            this->hypercubes.erase(it++);
        }
        else ++it;
    }


    return;

}


/** Retrieve the number of entities in the spatial region.
 *
 * @return The number of entities associated with the HyperSpace
 *  object.
 * */
unsigned HyperSpace::getNumEntities(void) const{


    unsigned num_entities = 0;


    // Sum the number of entities in all hypercubes
    hypercube_iterator iter = this->hypercubes.begin();
    while( iter != this->hypercubes.end() ){

        num_entities += iter->second.numObjects();
        iter++;
    }


    return num_entities;
}



/** Methods of class EntityIterator **/



// Constructor
HyperSpace::EntityIterator::EntityIterator( HyperSpace& hs) : space(&hs) {}

// Destructor
HyperSpace::EntityIterator::~EntityIterator(){}


// Copy-constructor
HyperSpace::EntityIterator::EntityIterator( const EntityIterator& other) : space(other.space) {

    this->cube_keys_iterator = other.cube_keys_iterator;
    this->entities_iterator = other.entities_iterator;
}


/** Move the cursor to the beginning of entities.
 *
 * */
void HyperSpace::EntityIterator::begin(){


    this->cube_keys_iterator = this->space->high_populated_keys.begin();

    this->entities_iterator = this->space->hypercubes.find( *(this->cube_keys_iterator) )->second.retrieveObjects().begin();

}


/** Move the cursor to the next entity.
 *
 * */
void HyperSpace::EntityIterator::operator++(int){



    // Verify whether iteration over current hypercube ended
    this->entities_iterator++;
    if( this->entities_iterator == this->space->hypercubes.find( *(this->cube_keys_iterator) )->second.retrieveObjects().end() ){


        // Restart iteration in the next hypercube, unless the end of iteration was reached

        this->cube_keys_iterator++;

        if( this->cube_keys_iterator != this->space->high_populated_keys.end() ){

            this->entities_iterator = this->space->hypercubes.find( *(this->cube_keys_iterator) )->second.retrieveObjects().begin();
        }

    }


}


/** Retrieve the entity that the cursor is pointing to.
 *
 * */
DatasetEntity& HyperSpace::EntityIterator::operator*(){

    return *entities_iterator;

}


/** Retrieve a pointer to the entity that the cursor is pointing to.
 *
 * @return a pointer to the entity pointed to by the cursor.
 * */
DatasetEntity* HyperSpace::EntityIterator::operator->(){


    return &(*entities_iterator);
}


/** Verify whether the cursor is at the end of the list of
 * entities.
 *
 * @return True, if the cursor is at the end of the list of
 *  entities. False, otherwise.
 * */
bool HyperSpace::EntityIterator::end(){


    return (this->cube_keys_iterator == this->space->high_populated_keys.end());
}


