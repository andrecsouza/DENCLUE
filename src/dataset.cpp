
/* 
 * Copyright 2006 Andre Cardoso de Souza
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
 *
 */










/* INCLUSIONS */
#include "dataset.h"


/* METHODS */



/** Methods from class DatasetEntity **/


/** Build an entity from its representation in a character
 * sequence.
 *
 * @param input string used to build a dataset entity.
 *
 * @return void
 *
 * */
void DatasetEntity::buildEntityFromString( const string& input ){


    istringstream string_input(input);   // Buffer used to read entity values more easily
    stringstream string_output( ios::in | ios::out );

    unsigned int curr_component_index = 0;

    /* Read each character in the input stream, building a component value each time
     * the csv separator is seen and when the stream ends. */

    char curr_char = (char) string_input.get();
    if( string_output.eof() ) return;


    do{

        if( curr_component_index >= this->num_dimensions ) break;


        // If the current character is the csv separator or end of line, obtain the component value
        if( (curr_char == Constants::CSV_SEPARATOR) || (curr_char == Constants::EOL) ){

            string str;
            string_output >> str;
            this->attributes[curr_component_index] = atof(str.c_str());  // Store the component value
            curr_component_index++;  // Go to next component
            string_output.clear();   // Restart character collecting

            if( string_input.eof() )  break;
        }

        else{  // Append the character to the string representing the component value

            string_output << curr_char;
        }


        curr_char = (char) string_input.get();

    }while( !string_input.eof() ); // end of while


    return;
}


/** Retrieve a string representation for an entity
 * sequence.
 *
 * @return The string representation of the entity
 *
 * */
string DatasetEntity::getStringRepresentation( void ) const{


    ostringstream out_str;

    for(unsigned i=0 ; i < this->num_dimensions ; i++){

        if(i != 0) out_str << ',';
        out_str << this->attributes[i];
    }


    return out_str.str();
}



/** Retrieve the value of the i-th component of the entity.
 *
 *  @param component_index The index of the component whose value is required
 *
 * @return The value of the component specified by the index in the entity
 *
 * */
inline double DatasetEntity::getComponentValue( unsigned component_index ) const {


    // Verify whether the received index is compatible with the number of
    // dimensions of the entity
    if( component_index >= this->num_dimensions ){

        cerr << "[DatasetEntity::getComponentValue] Incompatible index (" <<
            component_index << "); entity has dimension " << this->num_dimensions << endl;
    }


    return this->attributes[component_index];


}

/** Calculate the difference between two DatasetEntity.
 *
 *  @param operand 2nd operand of difference
 *
 * @return a DatasetEntity representing the difference between the 2
 *  operands
 * */
DatasetEntity DatasetEntity::operator-( const DatasetEntity& operand ) const {


    unsigned dimension_output = (this->num_dimensions < operand.num_dimensions) ? this->num_dimensions : operand.num_dimensions ;
    DatasetEntity difference(dimension_output);


    ostringstream entity_str;

    // Calculate the difference of each component
    for(unsigned i=0 ; i < dimension_output ; i++){

        double curr_difference = this->getComponentValue(i) - operand.getComponentValue(i);

        if( i!= 0 ) entity_str << ',';
        entity_str << curr_difference;


    }
    entity_str << Constants::EOL;


    difference.buildEntityFromString( entity_str.str() );

    return difference;
}


/** Calculate the sum of two DatasetEntity.
 *
 *  @param operand 2nd operand of sum
 *
 * @return a DatasetEntity representing the sum between the 2
 *  operands
 * */
DatasetEntity DatasetEntity::operator+( const DatasetEntity& operand ) const{


    unsigned dimension_output = (this->num_dimensions < operand.num_dimensions) ? this->num_dimensions : operand.num_dimensions ;
    DatasetEntity sum(dimension_output);


    ostringstream entity_str;

    // Calculate the difference of each component
    for(unsigned i=0 ; i < dimension_output ; i++){

        double curr_sum = this->getComponentValue(i) + operand.getComponentValue(i);

        if( i!= 0 ) entity_str << ',';
        entity_str << curr_sum;
    }
    entity_str << Constants::EOL;

    sum.buildEntityFromString( entity_str.str() );

    return sum;
}


/** Calculate the Euclidean norm of a DatasetEntity.
 *
 * @return value of the Euclidean norm of the entity
 * */
double DatasetEntity::getEuclideanNorm() const {


    double component_squares_sum = 0;

    for(unsigned i=0 ; i < this->num_dimensions ; i++){

        component_squares_sum += pow(this->getComponentValue(i), 2.0);
    }


    double norm = sqrt(component_squares_sum);  /*/ (this->num_dimensions*1.0);*/

    return norm;
}


/** Calculate the distance between two dataset entities.
 *
 *  @param entity1 First operand.
 *  @param entity2 Second operand.
 *
 * @return The Euclidean distance between entities.
 * */
double DatasetEntity::distanceBetween( const DatasetEntity& entity1, const DatasetEntity& entity2 ){


    const DatasetEntity difference = entity1 - entity2;
    double distance = difference.getEuclideanNorm();

    return distance;
}



/** Methods from class Dataset **/



/** Instantiates a Dataset object. Initializes some attributes.
 *
 * */
Dataset::Dataset(unsigned dimensions){


    this->dimensions = dimensions;

    // Allocate all arrays in the object
    this->sum = new double[dimensions];
    this->upper_bound = new double[dimensions];
    this->lower_bound = new double[dimensions];



    // Zeroes all arrays in the object
    memset((void *) this->sum, 0, dimensions * sizeof(double));
    memset((void *) this->upper_bound, (char)(((unsigned)0xffff)* -1.0), dimensions * sizeof(double));
    memset((void *) this->lower_bound, (char)(((unsigned)0xffff)*1.0), dimensions * sizeof(double));

    return;
}


// Copy-constructor
Dataset::Dataset(const Dataset& other) : dimensions(other.dimensions){


    // Allocate all arrays in the object
    this->sum = new double[this->dimensions];
    this->upper_bound = new double[this->dimensions];
    this->lower_bound = new double[this->dimensions];

    // Copy arrays
    for(unsigned i=0 ; i < this->dimensions ; i++){

        this->sum[i] = other.sum[i];
        this->lower_bound[i] = other.lower_bound[i];
        this->upper_bound[i] = other.upper_bound[i];
    }


    // Copy entities
    Dataset::iterator it(other);
    for( it.begin() ; !it.end() ; it++){
        DatasetEntity entity = other.getEntity(*it);
        this->addEntity( entity );
    }


    // Copy hypercubes
    Dataset::hypercube_container::const_iterator iter = other.hypercubes.begin();
    for( ; iter != other.hypercubes.end() ; iter++){
        this->hypercubes.insert(*iter);
    } // */

}


/** Destroy a Dataset object.
 *
 * */
Dataset::~Dataset(){

    this->hypercubes.clear();
    this->entities.clear();

    // Free storage
    delete[] this->sum;
    delete[] this->upper_bound;
    delete[] this->lower_bound;

}



/** Insert an entity to this dataset.
 *
 *  @param entity The entity to insert.
 *
 * */
void Dataset::addEntity( const DatasetEntity& entity){

    // Push the entity back to the entities vector
    this->entities.push_back(entity);


    // Add each component value to the array of component sums
    unsigned int i = 0;
    for( i=0 ; i < entity.getNumOfDimensions() ; i++){

        this->sum[i] += entity.getComponentValue(i);

        // Updates the upper and lower bounds of the dataset
        this->upper_bound[i] = entity.getComponentValue(i) > this->upper_bound[i] ? entity.getComponentValue(i) : this->upper_bound[i];
        this->lower_bound[i] = entity.getComponentValue(i) < this->lower_bound[i] ? entity.getComponentValue(i) : this->lower_bound[i];


        // To avoid problems with precision loss, round bounds' values
        this->upper_bound[i] = ceill( this->upper_bound[i] );
        this->lower_bound[i] = floorl( this->lower_bound[i] );

    }

}



/** Retrieve the upper bound of each dataset component.
 *
 * @return the vector of upper bounds the dataset
 *  */
const vector<double>& Dataset::retrieveUpperBound() const {


    vector<double> *upper_bounds = new vector<double>();


    // Copy values of upper bounds stored in the dataset
    for(unsigned i=0 ; i < this->dimensions ; i++){

        upper_bounds->push_back(this->upper_bound[i]);
    }


    return (*upper_bounds);

}


/** Retrieve the lower bound of each dataset component.
 *
 * @return the vector of lower bounds the dataset
 *  */
const vector<double>& Dataset::retrieveLowerBound() const {


    vector<double> *lower_bounds = new vector<double>();


    // Copy values of lower bounds stored in the dataset
    for(unsigned i=0 ; i < this->dimensions ; i++){

        lower_bounds->push_back(this->lower_bound[i]);
    }


    return (*lower_bounds);

}


/** Calculate the product of one DatasetEntity by one scalar value.
 *
 *  @param scalar Number used to make the product of the entity
 *
 * @return a DatasetEntity representing the product
 *  */
DatasetEntity DatasetEntity::operator*( long double scalar ) const {



    DatasetEntity scaled(*this);


    for(unsigned i=0; i < scaled.getNumOfDimensions() ; i++){

        scaled.attributes[i] *= scalar;
    }

    return scaled;
}



