


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
#include "denclue_functions.h"


/* METHODS */




/** Calculate the influence of an entity in another. The chosen
 * influence function was the Gaussian Influence Function, defined by:
 * I(x,y) = exp { - [distance(x,y)**2] / [2*(sigma**2)] }
 *
 *  @param entity1 First entity of the operation
 *  @param entity2 Second entity of the operation
 *  @param sigma Parameter that ponderates the influence of an entity into another
 *
 * @return The value of influence.
 * */
long double DenclueFunctions::calculateInfluence( const DatasetEntity&
        entity_one, const DatasetEntity& entity_two , double sigma){



    long double distance = DatasetEntity::distanceBetween(entity_one, entity_two);

    // Verify whether the entities are the same (indirectly)
    if( distance == 0 ){

        return 0;  // Influence is zero if entities are the same
    }



    long double exponent = - powl( distance , 2) / (2.0 * powl(sigma,2) );
    long double influence = expl(exponent);


    return influence;
}



/** Calculate the density in an entity. It's defined as the sum of the
 * influence of each another entity of dataset.
 *
 *  @param entity The entity to calculate density
 *  @param iter Iterator over entities of dataset.
 *  @param sigma Parameter that ponderates the influence of an entity into another
 *
 * @return The value of density in entity.
 * */
long double DenclueFunctions::calculateDensity( const DatasetEntity& entity, HyperSpace::EntityIterator iter , double sigma){


    long double density = 0;

    while( !iter.end() ){

        density += DenclueFunctions::calculateInfluence( entity, *iter, sigma );
        iter++;
    }


    return density;
}



/** Calculate gradient of density functions in a given spatial point.
 *
 *  @param entity The spatial point used to calculate the gradient.
 *  @param iter Iterator over all dataset entities.
 *  @param sigma Parameter that ponderates the influence of an entity into another
 *
 * @return The vector that represents the gradient of the influence
 *  function in a given spatial point.
 * */
vector<double> DenclueFunctions::calculateGradient( const DatasetEntity& entity, HyperSpace::EntityIterator iter, double sigma ){


    vector<double> gradient;
    for( unsigned i=0 ; i < entity.getNumOfDimensions(); i++){

        gradient.push_back(0);
    }


    // Iterate over all entities and calculate the factors of gradient
    for( ; !iter.end() ; iter++){

        const DatasetEntity& other_entity = *iter;
        double curr_influence = DenclueFunctions::calculateInfluence(entity, other_entity, sigma);

        // Calculate the gradient function for each dimension of data
        for(unsigned i=0 ; i < entity.getNumOfDimensions() ; i++){

            double curr_difference = (iter->getComponentValue(i) - entity.getComponentValue(i));
            gradient[i] += curr_difference * curr_influence;

        }

    }


    return gradient;

}


/** Find density-attractor for an entity. The density-attractor is
 * obtained executing a hill climbing algorithm.
 *
 *  @param entity The spatial point used to calculate the gradient.
 *  @param iter Iterator over all dataset entities.
 *  @param sigma Parameter that ponderates the influence of an entity into
 *  another
 *
 * @return An entity that is the density-attractor for the given
 *  entity.
 * */
const DatasetEntity DenclueFunctions::getDensityAttractor( const DatasetEntity& entity, const HyperSpace& spatial_region, HyperSpace::EntityIterator iter, double sigma ){


    const double delta = 1;

    HyperSpace::EntityIterator initial_iter = iter;

    //DatasetEntity last_attractor(entity); // Set the initial density-attractor to the received entity
    DatasetEntity curr_attractor(entity);
    DatasetEntity *found_attractor = NULL;


    // Execute the hill climbing algorithm until it finds the local maxima of density function
    unsigned MAX_ITERATIONS = 1000;
    bool reachedTop = false;
    do{

        // Avoid infinite loops
        if( --MAX_ITERATIONS <= 0 )  break;


        // Store last calculated values for further comparison
        DatasetEntity last_attractor(curr_attractor);


        // Calculate the gradient of density function at current candidate to attractor
        HyperSpace::EntityIterator gradient_iter = initial_iter;
        vector<double> curr_gradient =
            DenclueFunctions::calculateGradient(last_attractor, initial_iter,
                    sigma);


        // Build an entity to represent the gradient
        ostringstream grad_entity_str;
        for(unsigned i=0 ; i < curr_gradient.size() ; i++){

            if( i != 0 ) grad_entity_str << Constants::CSV_SEPARATOR;
            grad_entity_str << curr_gradient[i];

        }
        grad_entity_str << Constants::EOL;

        DatasetEntity grad_entity(entity.getNumOfDimensions());
        grad_entity.buildEntityFromString(grad_entity_str.str());


        // Calculate next candidate to attractor
        long double grad_entity_norm = grad_entity.getEuclideanNorm();


        assert( grad_entity_norm > 0 );

        curr_attractor = last_attractor + ( ( (long double)(delta/grad_entity_norm)) * grad_entity );


        // Calculate density in current attractor
        HyperSpace::EntityIterator density_iter = initial_iter;
        density_iter.begin();
        double curr_density = calculateDensity( curr_attractor, density_iter, sigma );

        curr_attractor.setDensity(curr_density);

        // Verify whether local maxima was found
        reachedTop = ( curr_attractor.getDensity() < last_attractor.getDensity() );
        if( reachedTop ) found_attractor = new DatasetEntity(last_attractor);


    }while( !reachedTop );

    if( MAX_ITERATIONS <= 0 )  found_attractor = new DatasetEntity(curr_attractor);


    return *found_attractor;
}


/** Verify whether an appropriate path between two density-attractor
 * exists.  Each entity in the path MUST satisfy the minimum density
 * restriction.
 *
 *  @param attractor1 Attractor where the path must start
 *  @param attractor2 Attractor where the path must end
 *  @param xi Minimum density threshold
 *
 * @return True, if is possible to establish a path between the given
 * density-attractors.  False, otherwise.
 * */
bool DenclueFunctions::pathBetweenExists( const DatasetEntity& attractor1,
        const DatasetEntity& attractor2, HyperSpace& hs,  double xi, double
        sigma, map<string, bool>& usedEntities ){


    usedEntities[ attractor1.getStringRepresentation() ] = true;
    usedEntities[ attractor2.getStringRepresentation() ] = true;


    /* If the distance between entities is less or equal to sigma, a path can
     * be established between them */
    if( DatasetEntity::distanceBetween(attractor1, attractor2) <= sigma ){

        //DEBUG
        //cout << "Path found between " << attractor1 << " and " << attractor2 << endl;

        return true;
    }

    //DEBUG
    //cout << "Trying to find path between " << attractor1 << " and " << attractor2 << endl;


    vector<HyperSpace::EntityIterator>  curr_path;  // Path being constructed


    /* Try to create a path between received density-attractors using all
     * possibilities (backtracking). This is ugly, but it works.
     * */
    HyperSpace::EntityIterator iter(hs);
    iter.begin();
    //while( !iter.end() )
    while( usedEntities.size() <= hs.getNumEntities() ){



        /* If we iterated over all entities ... */
        if( iter.end() ){


            if( curr_path.empty() )  return false;  // Couldn't estalish even a path of size 1

            // Back one step if there's an entity to add to the path
            //if( usedEntities.size() != hs.getNumEntities() )
            if( curr_path.size() < (hs.getNumEntities() - 2) ){


                // Mark current path end as unused, remove it from the path and
                // move forward the cursor of the new path end
                usedEntities[ curr_path.back()->getStringRepresentation() ] = false;
                curr_path.pop_back();


                // If there isn't more possibilities, return fasse
                if( curr_path.empty() )  return false;


                // Find an unused entity
                bool isUsed = false;
                bool reachedEnd = false;
                do{
                    curr_path.back()++;  // Move cursor forward

                    reachedEnd = ( curr_path.back().end() );
                    if( reachedEnd ) break;  // Will back one more step

                    isUsed = usedEntities[
                        curr_path.back()->getStringRepresentation() ];

                }while( isUsed );


                iter = curr_path.back();  // Update loop iterator
                if( reachedEnd ) continue;


                if( !curr_path.empty() ){
                    curr_path.pop_back();
                }
                else  return false;

                //DEBUG
                /*cout << "Curr temp path: ";
                  cout << attractor1;
                  for(unsigned i=0; i < curr_path.size() ; i++){
                  cout << " -> ";
                  cout << *(curr_path[i]);
                  }
                  cout << endl; // */


            }

            else{  // If all possibilities have been tried, return false

                //DEBUG
                //cout << "No path found" << endl;

                return false;
            }


        }



        /* Avoid using the same entity twice */
        const DatasetEntity& curr_entity = *iter;
        if( usedEntities.count(curr_entity.getStringRepresentation()) > 0 ){

            if( usedEntities[curr_entity.getStringRepresentation()] ){
                iter++;
                continue;
            }
        }


        // Verify whether next entity can be part of the path
        const DatasetEntity& curr_path_end = curr_path.empty() ? attractor1 : *curr_path.back() ;
        if( (curr_entity.getDensity() >= xi) && ( DatasetEntity::distanceBetween(
                        curr_path_end, curr_entity ) < sigma ) ){


            // Add current entity to path and mark it as used
            curr_path.push_back( iter );
            usedEntities[curr_entity.getStringRepresentation()] = true;

            //DEBUG
            /*cout << "Curr path: ";
              cout << attractor1;
              for(unsigned i=0; i < curr_path.size() ; i++){
              cout << " -> ";
              cout << *(curr_path[i]);
              }
              cout << endl; // */

            // Verify whether 'attractor2' can be reached by recently inserted
            // entity
            if( DatasetEntity::distanceBetween( curr_entity , attractor2 ) < sigma ){


                //DEBUG
                /*cout << "Path found: ";
                  cout << attractor1;
                  for(unsigned i=0; i < curr_path.size() ; i++){
                  cout << " -> ";
                  cout << *(curr_path[i]);
                  }
                  cout << " -> " << attractor2 << endl; // */



                return true;
            }

            else{  // Restart iterator and try to add one more entity to the path

                iter.begin();
                continue;
            }


            /* Try to find a path between the recently inserted entity and
             * 'attractor2' */
            /*if( DenclueFunctions::pathBetweenExists( curr_entity,
              attractor2, hs, xi, sigma, usedEntities ) ){

            //DEBUG
            //cout << curr_entity << "->";

            return true;
            }
            else{

            // Remove current entity from path and mark it as unused
            usedEntities[curr_entity.getStringRepresentation()] = false;

            }*/

        }

        iter++;
    }


    return false;
}


/** Append one vector to the end of another.
 *
 *  @param dest Vector that will receive new elements.
 *  @param src Vector with elements to be appended to vector 'dest'
 *
 * */
void DenclueFunctions::AppendVector( vector<DatasetEntity>& dest, const
        vector<DatasetEntity>& src){


    for( unsigned i=0 ; i < src.size() ; i++){

        dest.push_back(src[i]);
    }

    return;
}



