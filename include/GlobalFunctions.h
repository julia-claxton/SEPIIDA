#ifndef __GLOBALFUNCTIONS_H__
#define __GLOBALFUNCTIONS_H__ 1

#include "G4RunManager.hh"

inline void printVector(G4ThreeVector v, G4double unit){
  G4cout << "(" << v.x()/unit << ", " << v.y()/unit << ", " << v.z()/unit << ")";
}

inline void printVector(G4ThreeVector v){
  G4cout << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
}

#endif