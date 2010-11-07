



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



#ifndef HYPERSPACE_H
#define HYPERSPACE_H


/* INCLUSIONS */
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <utility>
#include "hypercube.h"
#include "dataset.h"
using namespace std;


/* CLASSES */

class DatasetEntity;
class HyperCube;

/** @class HyperSpace
 *
 * @brief This class represents a multi-dimensional space. This space is divided into
 * hypercubes.
 *
 *
 * */
class HyperSpace {


    private:

        const unsigned dimension;

        /* Influence of a point in its neighborhood */
        const double sigma;

        /* Lower bound of density to consider */
        const double xi;

        /* Bounds of the multi-dimensional space */
        double *upper_bounds;
        double *lower_bounds;



    public:

        typedef DatasetEntity dataset_entity;
        typedef HyperCube space_hypercube;
        typedef map< string, space_hypercube> hypercube_container;
        typedef map< string, space_hypercube>::const_iterator hypercube_iterator;

    private:
        map< string, space_hypercube >   hypercubes;  // Regions in the space
        vector<string>   high_populated_keys;  // Regions in the space that satisfy entities minimum bound



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
        void createHyperCubes( const unsigned *num_partitions , const string&
                upp_bounds_str  , double edge_length);


        /** Retrieve the length of a partition (an edge of a hypercube).
         *
         * @return the length of a hypercube edge
         * */
        double hypercubeEdgeLenght() const {  return (2 * this->sigma);  }


        /** Retrieve the minimum number of entities of a high populated hypercube.
         *
         * @return the minimum number of entities of a high populated hypercube.
         * */
        double minimumObjectsInHypercubes() const {  return (this->xi / (2 *
                    this->dimension) );  }

    public:

        // Constructor
        HyperSpace( const vector<double>& up_bound, const vector<double>&
                lw_bound, double sigma , double xi, const unsigned num_dimensions);


        // Copy-constructor
        HyperSpace(const HyperSpace&);


        // Destructor
        ~HyperSpace(){

            /* Free storage */
            delete[] this->upper_bounds;
            delete[] this->lower_bounds;

            this->hypercubes.clear();
        }


        /** Determine the regions of the space based on the parameter sigma.
         * The index string is formed using the appropriate function.
         *
         * @return a mapping of the multi-dimensional index string and the
         *  hypercubes.
         *
         * */
        const map< string, space_hypercube >* determineSpatialRegions();



        /** Insert a dataset entity in the space.
         *
         *  @param entity The entity to insert.
         *
         * */
        void insertEntity( const dataset_entity& entity );


        /** Remove low populated hypercubes, except those who are connected to
         * a high populated hypercube.
         *
         * */
        void removeLowPopulatedHypercubes();


        /** Retrieve the number of entities in the spatial region.
         *
         * @return The number of entities associated with the HyperSpace
         *  object.
         * */
        unsigned getNumEntities(void) const;


        /** @class HyperSpace::EntityIterator
         *
         * @brief This class represents an iterator over all entities of all
         * hypercubes in a spatial region.
         *
         * */
        class EntityIterator {


            private:
                HyperSpace* space;
                vector<string>::const_iterator cube_keys_iterator;
                vector<DatasetEntity>::iterator entities_iterator;

            public:

                // Constructor
                EntityIterator( HyperSpace& );

                // Destructor
                ~EntityIterator();

                // Copy-constructor
                EntityIterator( const EntityIterator& );


                /** Move the cursor to the beginning of entities.
                 *
                 * */
                void begin();


                /** Move the cursor to the next entity.
                 *
                 * */
                void operator++(int);


                /** Retrieve the entity that the cursor is pointing to.
                 *
                 * @return a reference to the entity pointed to by the cursor.
                 * */
                DatasetEntity& operator*();


                /** Retrieve a pointer to the entity that the cursor is pointing to.
                 *
                 * @return a pointer to the entity pointed to by the cursor.
                 * */
                DatasetEntity* operator->();


                /** Verify whether the cursor is at the end of the list of
                 * entities.
                 *
                 * @return True, if the cursor is at the end of the list of
                 *  entities. False, otherwise.
                 * */
                bool end();


        };  // End of class entity_iterator

        friend class EntityIterator;

};  // End of class HyperSpace



#endif


