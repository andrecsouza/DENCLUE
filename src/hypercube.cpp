


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
#include "hypercube.h"


// Constructor
HyperCube::HyperCube( unsigned dimensions, const double *upper_bounds, double edge_length) : dimensions(dimensions), edge_length(edge_length) {


    this->hypercube_key = HyperCube::getKeyFromArray( upper_bounds, dimensions, edge_length );

    // Zeroes sum of entities components
    for(unsigned i=0 ; i < this->dimensions ; i++){
        this->entities_sum.push_back(0);
    }

}


// Copy-constructor
HyperCube::HyperCube( const HyperCube& other ) : dimensions(other.dimensions), edge_length(other.edge_length) {



    // Copy hypercube key
    this->hypercube_key = other.hypercube_key;


    // Copy objects
    this->objects.clear();
    for(unsigned i=0 ; i < other.objects.size() ; i++){
        this->objects.push_back(other.objects[i]);
    }


    // Copy objects sum
    this->entities_sum.clear();
    for(unsigned i=0; i < other.entities_sum.size() ; i++){

        this->entities_sum.push_back( other.entities_sum[i] );
    }


    // Copy neighbors
    this->neighbors.clear();
    for(unsigned i=0 ; i < other.neighbors.size() ; i++){
        this->neighbors.push_back(other.neighbors[i]);
    }

}



/** Insert an object in the hypercube.
 *
 *  @param object The object to insert
 *
 * */
void HyperCube::addObject( const DatasetEntity& object ){



    double *upper_bounds = HyperCube::getArrayFromKey(this->hypercube_key, this->dimensions);
    double lower_bounds[this->dimensions];


    /* Verify whether this object is inside the region represented by tehis
     * hypercube */
    bool outside_hypercube = false;
    for(unsigned i=0 ; i < this->dimensions ; i++){


        // Calculate lower bound of current component
        lower_bounds[i] = upper_bounds[i] - this->edge_length;


        double curr_comp_value = object.getComponentValue(i); // Value of i-th component
        if( (curr_comp_value < lower_bounds[i]) || (curr_comp_value >= upper_bounds[i]) ){

            // Object outside this spatial region
            cerr << "Entity " << object << " isn't inside this HyperCube" << endl;
            cerr << "Component " << i << " should be in [" << lower_bounds[i];
            cerr << "," << upper_bounds[i] << ")" << endl;

            outside_hypercube = true;
            break;
        }


    }


    // Free some storage
    delete[] upper_bounds;

    if( !outside_hypercube ){

        this->objects.push_back(object);  // Add object to hypercube

        // Update sum of entities components
        for(unsigned i=0 ; i < this->dimensions; i++){

            this->entities_sum[i] += object.getComponentValue(i);
        }

    }

    return;

}



/** Retrieve all objects in the hypercube.
 *
 * @return a vector with all objects in the hypercube
 * */
vector<DatasetEntity>& HyperCube::retrieveObjects(){

    return this->objects;

}


/** Assign a set of neighbors to this HyperCube. A representation of
 * each neighbor is stored, not the objects themselves. This
 * representation can be used as a key to retrieve a neighbor later.
 *
 *  @param neighbors set of neighbors of this spatial region.
 *
 * */
void HyperCube::setNeighbors( const vector<string>& neighbors ){


    // Copy each value in the vector to the internal vector
    vector<string>::const_iterator it = neighbors.begin();

    while( it != neighbors.end() ){


        bool isntSameHyperCube = ( *it != this->hypercube_key );


        if( isntSameHyperCube ){

            this->neighbors.push_back(*it);
        }
        it++;

    }

    return;
}


/** Create a string representation of a hypercube identifier from an
 * array.
 *
 *  @param array_key Array containing the upper bounds of a hypercube.
 *  @param dimension The number of dimensions in the hypercube.
 *  @param edge_length size of each uni-dimensional edge
 *
 * @return the string representation of the key
 *
 * */
string HyperCube::getKeyFromArray( const double *upper_bounds, unsigned dimension , double edge_length) {


    ostringstream str_buf;

    for(unsigned i=0 ; i < dimension ; i++){


        double curr_index =  upper_bounds[i] ;

        if( i != 0 )  str_buf << ',';
        str_buf << curr_index;
    }


    return str_buf.str();

}



/** Create an array representation of a hypercube identifier from a
 * string.
 *
 *  @param str_key String containing the key of a hypercube.
 *  @param dimension The number of dimensions in the hypercube.
 *
 * @return the upper bounds of the hypercube represented by the key.
 *  This array have 'dimension' elements
 *
 * */
double * HyperCube::getArrayFromKey( const string& str_key, unsigned dimension ){


    double *upp_bounds = new double[dimension];


    istringstream string_input(str_key);   // Buffer used to read array values more easily
    stringstream string_output;

    unsigned i = 0;


    /* Read each character in the input stream, building a value each time the
     * csv separator is seen and when the stream ends. This phase obtains a
     * multi-dimensional index that should be converted in bound values. */
    while( !string_input.eof() && (i < dimension) ){


        char curr_char = (char) string_input.get();


        // If the current character is the csv separator or end of line, obtain the component value
        if( (curr_char == Constants::CSV_SEPARATOR) || string_input.eof() ){

            string str;
            string_output >> str;
            upp_bounds[i] = atof(str.c_str());  // Store the component value


            i++;  // Go to next component
            string_output.clear();   // Restart character collecting

        }

        else{  // Append the character to the string representing the component value

            string_output << curr_char;
        }

    } // end of while


    return upp_bounds;

}


/** Retrieve a vector with all neighboring hypercubes' keys.
 *
 * @return a vector with all neighboring hypercubes' keys.
 * */
const vector<string>& HyperCube::getNeighbors(){


    /*vector<string> *neighbors = new vector<string>();

    // Copy each key to the output vector
    vector<string>::const_iterator it = this->neighbors.begin();
    for( ; it != this->neighbors.end() ; ++it){

    neighbors->push_back(*it);
    }*/

    return this->neighbors;
}


/** Remove keys of neighbors that are empty neighbors.
 *
 *  @param empty_neighbors: Vector with keys of empty neighbors
 *
 * */
void HyperCube::removeEmptyNeighbors( const vector<string>& empty_neighbors ){


    // Search for any of this hypercube's neighbors in the empty neighbors' list
    vector<string>::iterator it = this->neighbors.begin();
    while( it != this->neighbors.end() ){

        // Search for current neighbor in empties' list
        vector<string>::const_iterator wanted_key = find(
                empty_neighbors.begin() , empty_neighbors.end() , *it );


        if( wanted_key != empty_neighbors.end() ){  // Curr neighbor is and empty neighbor

            it = this->neighbors.erase(it);  // delete empty neihgbor
            continue;
        }

        it++;

    }


}


/** Verify whether this hypercube is neighbor of any of a list of hypercubes.
 *
 *  @param hypercube_keys Keys of hypercubes supposed to be neighbors.
 *  @param cubes Container of hypercubes.
 *
 * @return True, if any of the received list is a neighbor. False, otherwise.
 *
 * */
bool HyperCube::isNeighbor( const vector<string>& hypercube_keys , const HyperSpace::hypercube_container& cubes) const {

    bool is_neighbor = false;


    vector<string>::const_iterator it = hypercube_keys.begin();
    for( ; it != hypercube_keys.end() ; it++){


        // Search for current received key in the list of neighbors
        vector<string>::const_iterator wanted_key = find( this->neighbors.begin(), this->neighbors.end(), *it);
        bool found_it = (wanted_key != this->neighbors.end() );


        // More restrict neighborhood criterion: distance between means MUST be
        // less or equal to 2*edge_length
        if( found_it ){

            // Calculate distance between means of hypercubes
            DatasetEntity this_mean = this->getMeanElement();

            const HyperCube& neighbor_cube = cubes.find(*wanted_key)->second;
            DatasetEntity neighbor_mean = neighbor_cube.getMeanElement();

            DatasetEntity difference = this_mean - neighbor_mean;
            double distance = difference.getEuclideanNorm();


            // Verify whether distance satisfies minimum bound
            is_neighbor = ( distance <= (2*this->edge_length) );


            if( is_neighbor ) break;

        }

    }


    return is_neighbor;

}


/** Get the mean element of the hypercube.
 *
 * @return a DatasetEntity representing the mean of the hypercube.
 * */
DatasetEntity HyperCube::getMeanElement() const{


    DatasetEntity mean(this->dimensions);
    const unsigned num_entities = this->numObjects();

    // Calculate the mean values of each component, storing them in a string
    ostringstream mean_str;
    vector<double>::const_iterator it = this->entities_sum.begin();

    for( ; it != this->entities_sum.end() ; it++){

        double curr_component_mean = *it / (num_entities*1.0);

        if( it != this->entities_sum.begin() )  mean_str << ',';
        mean_str << curr_component_mean;
    }
    //mean_str << Constants::EOL;

    mean.buildEntityFromString( mean_str.str() );

    return mean;
}


