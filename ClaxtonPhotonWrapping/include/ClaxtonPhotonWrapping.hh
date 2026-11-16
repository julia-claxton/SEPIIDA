#ifndef ClaxtonPhotonWrapping_h
#define ClaxtonPhotonWrapping_h 1

#include "G4Gamma.hh"
#include "G4VEmProcess.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4ParticleDefinition;
class G4VEmModel;
class G4PropagatorInField;

class ClaxtonPhotonWrapping : public G4VDiscreteProcess
{
  public:  // with description
    explicit ClaxtonPhotonWrapping(const G4String& processName = "photonWrapping",
                              G4ProcessType type = fElectromagnetic);
    ~ClaxtonPhotonWrapping() override = default;

    G4double GetMeanFreePath(const G4Track& track, G4double previousStepSize,
                             G4ForceCondition* condition) override;

    G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& Step) override;

    G4bool IsApplicable(const G4ParticleDefinition&) override;
    void BuildPhysicsTable(const G4ParticleDefinition&) override;

    void ProcessDescription(std::ostream&) const override;
    void DumpInfo() const override { ProcessDescription(G4cout); }

    void SetSurfaceRoughness(const G4double value);

    G4double Reflectivity(const G4double GamEner, const G4double SinIncidentAngle,
                          const G4Material* theMat) const;

    G4int ReadHenkeXrayData(std::string ElName, std::vector<G4double>& Ephot,
                            std::vector<G4double>& f1, std::vector<G4double>& f2);
    void SaveHenkeDataAsMaterialProperty();  // Save Henke data as Material propoerty

  private:
    G4VPhysicalVolume* fLastVolume = nullptr;
    G4ThreeVector fSurfaceNormal;
    static G4double fSurfaceRoughness;  // same value in all threads
};

#endif
