


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





#ifndef HYPERCUBE_H
#define HYPERCUBE_H


/* INCLUSIONS */
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include "dataset.h"
#include "constants.h"
#include "hyperspace.h"
using namespace std;


/* CLASSES */

class DatasetEntity;
class HyperSpace;

/**  @class HyperCube
 *
 * @brief This class represents a spatial region determined by a hypercube. The
 * hypercube is used to store objects associated to it. The hypercube is
 * determined by ranges of values of spatial dimensions.
 *
 * */
class HyperCube{


    private:

        /*** Attributes ***/
        const unsigned dimensions;
        const double edge_length;

        string hypercube_key;   // String representation of the upper bounds of the cube

        vector< DatasetEntity > objects;  // Objects associated to the hypercube

        vector< string > neighbors;   // HyperCubes adjacent to this spatial region

        vector< double > entities_sum;  // Sum of each entity component. It speeds hypercube mean calculation


    public:

        /*** Instance methods ***/

        // Constructor
        HyperCube( unsigned dimensions, const double *upper_bounds, double edge_length);

        // Copy-constructor
        HyperCube( const HyperCube& );


        // Destructor
        ~HyperCube(){

            this->objects.clear();
            this->neighbors.clear();

        }



        /** Retrieve the number of objects in the hypercube.
         *
         * @return the number of objects in the hypercube
         * */
        unsigned numObjects() const {  return this->objects.size();  }


        /** Print a representation of the hypercube.
         * */
        friend ostream& operator<<(ostream& os, const HyperCube& hc){

            cout << "-------- HyperCube ---------" << endl;
            cout << "Dimensions: " << hc.dimensions << endl;
            cout << "Edge lenght: " << hc.edge_length << endl;
            cout << "Key: (" << hc.hypercube_key << ')' << endl;


            cout << "# of entities: " << hc.numObjects() << endl;

            cout << "Entities sum: (";
            for(unsigned i=0 ; i < hc.entities_sum.size() ; i++){

                if( i != 0 )  cout << ',';
                cout << hc.entities_sum[i];
            }
            cout << ")" << endl;


            cout << "Neighbors: ";
            for(vector<string>::const_iterator it = hc.neighbors.begin() ; it != hc.neighbors.end() ; it++){
                if( it != hc.neighbors.begin() )  cout << " , ";
                cout << *it;
            }
            cout << endl;

            cout << "-------- end ---------" << endl;

            return os;
        }


        /** Insert an object in the hypercube.
         *
         *  @param object The object to insert
         *
         * */
        void addObject( const DatasetEntity& object );


        /** Retrieve all objects in the hypercube.
         *
         * @return a vector with all objects in the hypercube
         * */
        vector<DatasetEntity>& retrieveObjects();


        /** Assign a set of neighbors to this HyperCube. A representation of
         * each neighbor is stored, not the objects themselves. This
         * representation can be used as a key to retrieve a neighbor later.
         *
         *  @param neighbors set of neighbors of this spatial region.
         *
         * */
        void setNeighbors( const vector<string>& neighbors );


        /** Retrieve a vector with all neighboring hypercubes' keys.
         *
         * @return a vector with all neighboring hypercubes' keys.
         * */
        const vector<string>& getNeighbors();


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
        static string getKeyFromArray( const double *upper_bounds, unsigned dimension , double edge_length);


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
        static double * getArrayFromKey( const string& str_key, unsigned dimension );


        /** Verify whether the hypercube has no objects.
         *
         * @return true, if hypercube has no objects. False, otherwise.
         * */
        bool isEmpty() const {  return (this->numObjects() == 0 );   }


        /** Remove keys of neighbors that are empty neighbors.
         *
         *  @param empty_neighbors: Vector with keys of empty neighbors
         *
         * */
        void removeEmptyNeighbors( const vector<string>& empty_neighbors );


        /** Verify whether this hypercube is neighbor of any of a list of hypercubes.
         *
         *  @param hypercube_keys Keys of hypercubes supposed to be neighbors.
         *  @param cubes Container of hypercubes.
         *
         * @return True, if any of the received list is a neighbor. False, otherwise.
         *
         * */
        bool isNeighbor( const vector<string>& hypercube_keys , const map<string, HyperCube>& cubes) const ;


        /** Get the mean element of the hypercube.
         *
         * @return a DatasetEntity representing the mean of the hypercube.
         * */
        DatasetEntity getMeanElement() const;



};  // End of class Hypercube


#endif

