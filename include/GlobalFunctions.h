#ifndef __GLOBALFUNCTIONS_H__
#define __GLOBALFUNCTIONS_H__ 1

#include "G4RunManager.hh"

#define __DEBUG_PING__ G4cout << __FILE__ << ": " << __LINE__ << G4endl

inline void printVector(G4ThreeVector v, G4double unit){
  G4cout << "(" << v.x()/unit << ", " << v.y()/unit << ", " << v.z()/unit << ")";
}

inline void printVector(G4ThreeVector v){
  G4cout << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
}

inline void println(G4String line){
  G4cout << line << G4endl;
}





inline G4double HACKY_LATITUDE = -999.0;
inline G4String HACKY_FIELDMODEL = "throw an error";










#endif