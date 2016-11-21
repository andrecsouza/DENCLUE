


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





#ifndef DATASET_H
#define DATASET_H

/* INCLUSIONS */
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cstring>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include "hypercube.h"
#include "constants.h"
using namespace std;

class HyperCube;


/* CLASSES */



/** @Class DatasetEntity
 * @brief This class represents a single dataset entity. It must have the same
 * dimension of the enclosing dataset. Each attribute of the point have numeric
 * values.
 *
 * */
class DatasetEntity {

    private:

        /*** Object attributes ***/
        double *attributes;  // Values in columns of this data entity
        unsigned num_dimensions;  // Dimension of this dataset

        double density;

    public:


        /*** Instance methods ***/

        //Constructor
        DatasetEntity(unsigned dimension) : num_dimensions(dimension) {

            this->attributes = new double[dimension];
            this->density = 0;
        }


        // Copy-constructor
        DatasetEntity(const DatasetEntity& other) : num_dimensions(other.num_dimensions) {


            this->attributes = new double[other.num_dimensions];
            this->density = other.density;

            // Copy each element of the array ttributes'
            for(unsigned i=0 ; i < other.num_dimensions ; i++){

                this->attributes[i] = other.attributes[i];
            }

        }


        //Destructor
        ~DatasetEntity(){   delete[] this->attributes;    }


        /** Copy an exiting entity to another.
         *
         *  @param clone Entity that will receive the copy
         *  @param source Entity that will be copied.
         *
         * @return a reference to the entity copy.
         * */
        DatasetEntity& operator=( const DatasetEntity& copy ){


            this->num_dimensions = copy.num_dimensions;

            delete[] this->attributes;
            this->attributes = new double[copy.num_dimensions];
            this->density = copy.density;

            // Copy each element of the array attributes'
            for(unsigned i=0 ; i < copy.num_dimensions ; i++){

                this->attributes[i] = copy.attributes[i];
            }

            return *this;

        }




        /** Build an entity from its representation in a character
         * sequence.
         *
         * @param input string used to build a dataset entity.
         *
         * @return void
         *
         * */
        void buildEntityFromString( const string& input );


        /** Retrieve a string representation for an entity
         * sequence.
         *
         * @return The string representation of the entity
         *
         * */
        string getStringRepresentation( void ) const;


        /** Set the value of density for the entity.
         *
         *  @param density Value of density.
         *
         * */
        void setDensity( double density ){  this->density = density;  }


        /** Get the value of density for the entity.
         *
         * @return value of density for the entity.
         *
         * */
        double getDensity( void ) const {  return this->density;  }


        /** Retrieve the value of the i-th component of the entity.
         *
         *  @param component_index The index of the component whose value is required
         *
         * @return The value of the component specified by the index in the entity
         *
         * */
        double getComponentValue( unsigned int component_index ) const ;


        /** Retrieve the number of dimensions of the entity.
         *
         * @return the number of dimensions of the entity
         *
         * */
        unsigned int getNumOfDimensions() const { return this->num_dimensions; }


        /** Print a representation of the entity in the output stream.
         *
         *  @param os The stream to write
         *  @param entity The representation of the entity
         *
         * @return: a reference to the received output stream
         * */
        friend ostream& operator<<(ostream& os, const DatasetEntity& entity) {

            // Print each component value
            os << '(';
            for(unsigned int i=0; i < entity.getNumOfDimensions() ; i++){

                if( i!= 0 ) os << ',';

                os << entity.getComponentValue(i);

            }
            os << ")" << " DENSITY ["<< entity.getDensity() <<"]";

            return os;

        }


        /** Calculate the difference between two DatasetEntity.
         *
         *  @param operand 2nd operand of difference
         *
         * @return a DatasetEntity representing the difference between the 2
         *  operands
         * */
        DatasetEntity operator-( const DatasetEntity& operand ) const;


        /** Calculate the sum of two DatasetEntity.
         *
         *  @param operand 2nd operand of sum
         *
         * @return a DatasetEntity representing the sum between the 2
         *  operands
         * */
        DatasetEntity operator+( const DatasetEntity& operand ) const;


        /** Calculate the product of one DatasetEntity by one scalar value.
         *
         *  @param scalar Number used to make the product of the entity
         *
         * @return a DatasetEntity representing the product
         *  */
        DatasetEntity operator*( long double scalar ) const;


        /** Verify whether an entity is equal to another
         *
         * @param other Entity to compare
         * @return True, if entities are equal. False, otherwise.
         * */
        bool operator==( const DatasetEntity& other ) const {

            return (DatasetEntity::distanceBetween(*this, other) == 0);
        }



        /** Verify whether an entity is different of another
         *
         * @param other Entity to compare
         * @return True, if entities are different. False, otherwise.
         * */
        bool operator!=( const DatasetEntity& other ) const {   return !( *this == other );   }


        /** Calculate the product of one DatasetEntity by one scalar value.
         *
         *  @param scalar Number used to make the product of the entity
         *  @param Entity to scale
         *
         * @return a DatasetEntity representing the product
         *  */
        friend DatasetEntity operator*( double scalar, const DatasetEntity& entity){
            return entity.operator*(scalar);
        }


        /** Compare two entities with the relational operator 'less than'.
         *
         *  @param other 2nd operand of comparation
         *
         * @return True, if this entity is less than 'other'.
         *  */
        bool operator<( const DatasetEntity& other ) const {


            for(unsigned i=0 ; i < this->getNumOfDimensions() ; i++){

                if( this->getComponentValue(i) != other.getComponentValue(i) ){
                    return (this->getComponentValue(i) < other.getComponentValue(i) );
                }
            }

            return true;
        }


        /** Calculate the Euclidean norm of a DatasetEntity.
         *
         * @return value of the Euclidean norm of the entity
         * */
        double getEuclideanNorm() const;


        /** Calculate the distance between two dataset entities.
         *
         *  @param entity1 First operand.
         *  @param entity2 Second operand.
         *
         * @return The Euclidean distance between entities.
         * */
        static double distanceBetween( const DatasetEntity& entity1, const DatasetEntity& entity2 );



}; // end of class DatasetEntity




/** @class Dataset
 * @brief This class represents a dataset. This is the data structure used for
 * storing dataset's contents and perform operations over dataset points.
 *
 * */
class Dataset {


    private:

        typedef  DatasetEntity dataset_entity;
        typedef HyperCube dataset_hypercube;

        /* Object attributes */
        unsigned dimensions;
        vector< dataset_entity > entities;   // entities of this dataset. The real data.

        map< int , dataset_hypercube > hypercubes;  // Mapping from entities into spatial regions

        double *sum;         // Sum of all values in each entity component.

        // Attributes that help to create a hypercube
        double *upper_bound;
        double *lower_bound;


    public:

        typedef map<int, dataset_hypercube> hypercube_container;

        // Constructor
        Dataset(unsigned dimension);

        // Copy-constructor
        Dataset(const Dataset&);

        // Destructor
        ~Dataset();


        /** Instance methods **/


        /** Insert an entity to this dataset.
         *
         *  @param entity The entity to insert.
         *
         * */
        void addEntity( const dataset_entity& entity);



        /** Retrieve an entity of this dataset.
         *
         *  @param entity The index of the entity.
         *
         * */
        dataset_entity getEntity( unsigned int index ) const {


            // Verify whether the index is valid
            if( index < this->entities.size() ){

                return this->entities[index];
            }
            else{

                // Print an error message
                cerr << " [Dataset::getEntity] Invalid index received (" << index << " of ";
                cerr << this->entities.size() << ")" << endl;

                return *(new DatasetEntity(dimensions));
            }


        }



        /** Retrieve the number of dimensions of this dataset.
         *
         * @return the number of dimensions of this dataset.
         *
         * */
        unsigned int getNumOfDimensions() const {  return dimensions; }



        /** Retrieve the number of entities in this dataset.
         *
         * @return the number of entities in this dataset.
         *
         * */
        unsigned int getNumOfEntities() const {  return this->entities.size(); }



        /** Retrieve the upper bound of each dataset component.
         *
         * @return the vector of upper bounds the dataset
         *  */
        const vector<double>& retrieveUpperBound() const;


        /** Retrieve the lower bound of each dataset component.
         *
         * @return the vector of lower bounds the dataset
         *  */
        const vector<double>& retrieveLowerBound() const ;



        /** @class Dataset::iterator
         *
         * @brief This class represents an iterator over the list of dataset's
         * entities.
         * */
        class iterator {


            private:
                const Dataset& dataset;
                unsigned int element_index;

            public:

                // Constructor
                iterator( const Dataset& ds ) : dataset(ds), element_index(0) {}


                // Destructor
                ~iterator(){}


                /** Instance methods and operators **/


                /** Initiate an iteration over the elements of the dataset
                 *
                 *  */
                void begin(){

                    this->element_index = 0;
                    return;
                }



                /** Move the cursor to the next element of the dataset. It's a
                 * postfix operator.
                 *
                 * @param int (indicates that the operator is postfix)
                 *  */
                void operator++(int){

                    this->element_index++;
                    return;
                }


                /** Retrieve the object pointed to by this iterator.
                 *
                 *  @param it Iterator pointing to the desired object.
                 *
                 * @return Index of current entity.
                 *  */
                friend unsigned int operator*( iterator& it ){

                    if( it.element_index < it.dataset.getNumOfEntities() ){
                        return it.element_index;
                    }
                    else return it.dataset.getNumOfEntities();
                }



                /** Verify whether the end of the list of elements was reached.
                 *
                 * @return True, if the end was reached. False, otherwise.
                 *  */
                bool end() const {

                    return (this->element_index >= this->dataset.entities.size() );
                }



        }; // end of class 'iterator'


        typedef iterator entities_iterator;


};  // end of class Dataset



#endif


