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

inline G4ThreeVector rotateVector(G4ThreeVector startingVector, G4ThreeVector rotateAbout, G4double rotationAngleDeg){
  G4ThreeVector v = startingVector; // To match nomenclature on wikipedia to reduce chances of a typo
  G4ThreeVector k = rotateAbout / rotateAbout.mag(); // Unit vector for axis of rotation
  G4double rotationAngleRad = rotationAngleDeg * M_PI / 180.0;

  // Euler/Rodrigues's finite rotation formula
  // https://en.wikipedia.org/wiki/Rodrigues'_rotation_formula#Statement
  G4ThreeVector rotatedVector = 
    (v * std::cos(rotationAngleRad))
    + (k.cross(v) * std::sin(rotationAngleRad))
    + (k * k.dot(v) * (1 - std::cos(rotationAngleRad)))
  ;
  return rotatedVector;
}


// Passing around variables between .cc files is hard, so I use the hacky solution of globals
inline G4double HACKY_LATITUDE = -999.0;
inline G4String HACKY_FIELDMODEL = "throw an error";
inline G4bool CACHED_MAGNETIC_FIELD;

#endif