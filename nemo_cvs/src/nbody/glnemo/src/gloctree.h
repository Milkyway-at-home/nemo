// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2007                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#ifndef GLOCTREE_H
#define GLOCTREE_H

#include <qgl.h>

#include "global_options.h"
#include "virtual_particles_select.h"
#include "gl_object.h"
#include "particles_data.h"
#include "frustumculling.h"
/**
@author Jean-Charles Lambert
*/
 
class Node {
  public:
    Node(int index_=-1, int obj_=-1);
    ~Node();
    int type;
    int index;
    int obj;
    Node * node[8];
};

class GLOctree : public GLObject {
  Q_OBJECT
  public:
    GLOctree(GlobalOptions * _options);

    ~GLOctree();
    void update(ParticlesData   *  ,
		ParticlesSelectVector * );
    void update();
    void displayPolygons(const double * mModel,const double * mProj, const GLuint texture,const float u_max,const float v_max);
  public slots:
    void buildDisplayList();
  private:
    Node * root;
    int level_max;     // level tree
    int level_requested;
    int size_max;      // max square size in power of 2
    bool new_data;     // true if new
    int nbody_keeped;  // #bodies keep after tree walk
    float rmid[3];
    ParticlesData * p_data;
    GlobalOptions * store_options;
    ParticlesSelectVector  * psv;
    bool first;
    void computeSizeMax();
    int build();
    void expandTree(float );
    void insertParticle(int, Node **, float *, int, int );
    int hackTreeDL(Node *, float *, int );
    int hackTreePolygons(Node *, float *, int, bool );
    int computePolygons(Node *);
    void init();

    const double * mModel;
    GLuint texture,textureX2;
    float  u_max;
    float  v_max;
    FrustumCulling frustum;
};

#endif
