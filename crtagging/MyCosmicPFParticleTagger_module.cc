////////////////////////////////////////////////////////////////////////
// Class:       MyCosmicPFParticleTagger
// Plugin Type: producer (art v2_05_00)
// File:        MyCosmicPFParticleTagger_module.cc
//
// Generated at Tue Jun  6 06:38:34 2017 by Andrew Smith using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/geo.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"

#include "larreco/RecoAlg/SpacePointAlg.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "TVector3.h"

#include <Eigen/Dense>

namespace crtagging {
  class MyCosmicPFParticleTagger;
}

class crtagging::MyCosmicPFParticleTagger : public art::EDProducer {
public:
  explicit MyCosmicPFParticleTagger(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MyCosmicPFParticleTagger(MyCosmicPFParticleTagger const &) = delete;
  MyCosmicPFParticleTagger(MyCosmicPFParticleTagger &&) = delete;
  MyCosmicPFParticleTagger & operator = (MyCosmicPFParticleTagger const &) = delete;
  MyCosmicPFParticleTagger & operator = (MyCosmicPFParticleTagger &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  enum FaceType {
    NONE       = 0,
    TOP        = 1,
    BOTTOM     = 2,  
    UPSTREAM   = 3,
    DOWNSTREAM = 4
  };

  // Declare member data here.

  // To be set from fhicl configuration
  std::string  m_PFParticleModuleLabel;  ///< Label for the PFParticle producer
  std::string  m_TrackModuleLabel;       ///< A label for a Track producer (is only used to match the functionality of the CRHitRemoval module)
  unsigned int m_MinHitsToConsider;      ///< Minimum number of 3D hits (space points) that a PFP must have for it to be considered by the tagger
  double       m_FaceAssociationThresh;  ///< Distance in cm that an endpoint must be to a detector face to be classed as "associated"
  double       m_OutOfTimeThresh;        ///< Distance in cm that an endpoint must be outside the detector to be classed as "out of time"

  // Information about the detector geometry
  double m_TPCExtentX;      ///< The width  of the detector along the X-axis ( anode - cathode )
  double m_TPCExtentY;      ///< The height of the detector along the Y-axis ( bottom - top )
  double m_TPCExtentZ;      ///< The length of the detector along the Z-axis ( upstream - downstream )

  double m_TPCCentralX;     ///< The central coordinate along the X-axis
  double m_TPCCentralY;     ///< The central coordinate along the Y-axis
  double m_TPCCentralZ;     ///< The central coordinate along the Z-axis

  double m_TPCFaceMinX;     ///< The X-coordinate of the face at minimum X ( anode )
  double m_TPCFaceMaxX;     ///< The X-coordinate of the face at maximum X ( cathode )
  double m_TPCFaceMinY;     ///< The Y-coordinate of the face at minimum Y ( bottom )
  double m_TPCFaceMaxY;     ///< The Y-coordinate of the face at maximum Y ( top )
  double m_TPCFaceMinZ;     ///< The Z-coordinate of the face at minimum Z ( upstream )
  double m_TPCFaceMaxZ;     ///< The Z-coordinate of the face at maximum Z ( downstream )

  // Private member functions


  /**
   *  @breif  GetPrimaryId  gets the ID of the primary PFP associated with the PFP index supplied
   *
   *  @param  pfParticleHandle  input handle to the list of PFParticles
   *  @param  thisId            input ID of the PFP in questions
   *
   *  @return the ID of the primary PFP associated with the PFP with thisId
   */
  size_t GetPrimaryId( art::Handle<std::vector<recob::PFParticle> > & pfParticleHandle, size_t thisId ) const;
  
  /**
   *  @brief  IsOutOfTime  checks if a PFParticle is out of time
   *
   *  @param  endpoint1  input position vector of one endpoint of a PFP
   *  @param  endpoint2  input position vector of the other endpoint of a PFP
   *  @param  thresh     input threshold ( in cm ) - to be out of time, and endpoint must be at least this distance out of the detector volume
   *
   *  @return if any of the end points are out of time
   */
  bool IsOutOfTime( TVector3 & endpoint1, TVector3 & endpoint2, double thresh ) const;

  /**
   *  @brief  GetAssociatedFaces  gives a list of faces associated to a given endpoint
   *
   *  @param  endpoint  input position vector of an end-point
   *  @param  enddir    input direction vector at the end-point
   *  @param  thresh    input distance threshold ( in cm ) to associate an endpoint to a face
   *  @param  faceList  output list of assocaited faces
   */ 
  void GetAssociatedFaces( TVector3 & endpoint, TVector3 & enddir, double thresh, std::vector< FaceType > & faceList ) const;

  /**
   *  @brief  CheckFace  check if the face given is in the face list given
   *
   *  @param  faceList  a list of detector faces associated with a PFP
   *  @param  type      the type of face to check
   *
   *  @return if a face of the supplied type is in the face list
   */
  bool CheckFace( std::vector< FaceType > & faceList, FaceType type ) const;

  /**
   *  @brief  CheckFaces  check if two endpoints given are associated with the two faces given
   *
   *  @param  faceList1  a list of detector faces associated with the first endpoint of a PFP
   *  @param  endpoint1  the first endpoint position of the PFP
   *  @param  faceList2  a list of detector faces associated with the second endpoint of a PFP
   *  @param  endpoint2  the seconed endpoint position of the PFP
   *  @param  typeUpper  the detector face to which the uppermost endpoint should be associated
   *  @param  typeLower  the detector face to which the lowermost endpoint should be associated
   *
   *  @return if both endpoints are associated to the faces supplied
   */
  bool CheckFaces( std::vector< FaceType > & faceList1, TVector3 & endpoint1, std::vector< FaceType > & faceList2, TVector3 & endpoint2, FaceType typeUpper, FaceType typeLower ) const;
  
  /**
   *  @brief  GetCosTheta  gets the opening cosine of the angle with the vertical by the line joining the two endpoints supplied
   *
   *  @param  endpoint1  position of the first endpoint of a PFP
   *  @param  endDir1    direction at the first endpoint of a PFP
   *  @param  endpoint2  position of the second endpoint of a PFP
   *  @param  endDir2    direction at the second endpoint of a PFP
   *
   *  @return the opening angle with the vertical
   */
  double GetCosTheta( TVector3 & endpoint1, TVector3 & endDir1, TVector3 & endpoint2, TVector3 & endDir2) const;

  /**
   *  @brief  GetEndpointVector  gets a std::vector of 3 floats corresponding to the position vector of the PFP with the ID supplied
   *
   *  @param  pfpToEndpointMap  input mapping between PFP id and endpoint position
   *  @param  thisId            input ID of the PFP in question
   *  @param  outputVect        output float vector - required for the hit removal algorithm
   */
  void GetEndpointVector( std::map< size_t, TVector3 > & pfpToEndpointMap, size_t thisId, std::vector<float> & outputVect ) const;

  /**
   *  @brief  SlidingFit class gives the result of a sliding fit to a supplied vector of space points
   *  
   *  An instance of this class performs a sliding fit to a supplied vector of space points.
   *
   *  First, a linear fit is performed on all of the space points to define a primary axis.
   *  A window size is chosen (as close as possible to the supplied size) which breaks the
   *  primary axis into an integer number of pieces. 
   *
   *  The window then slides along the primary axis, moving a distance half the window size
   *  with each iteration. For every position, a linear fit is performed on the space points 
   *  lying within the window. The central position and direction of the window is 
   *  stored along with the fit quality.
   */
  class SlidingFit {
  public:
    /**
     *  @brief  Constructor
     *
     *  @param  spacePointVect vector of space points to perform fit upon
     *  @param  fitWindow      sliding fit window ( in cm ) 
     */ 
    SlidingFit( std::vector< art::Ptr< recob::SpacePoint > > const & spacePointVect, double fitWindow = 3.0);
    
    /**
     *  @brief Destructor
     */ 
    ~SlidingFit() {}

    /**
     *  @brief  Get the first endpoint position
     */
    TVector3 GetEndpoint1() const { return m_endpoint1; }
    
    /**
     *  @brief  Get the first endpoint direction
     */
    TVector3 GetEndDir1() const { return m_enddir1; }
    
    /**
     *  @brief  Get the second endpoint position
     */
    TVector3 GetEndpoint2() const { return m_endpoint2; }
    
    /**
     *  @brief  Get the second endpoint direction
     */
    TVector3 GetEndDir2() const { return m_enddir2; }

    /**
     *  @brief  Get the curvature
     */
    double GetCurvature() const { return m_curvature; }

  private:
    // Private member variables
    std::vector< TVector3 > m_points;         ///< Full vector of space points, converted to TVector3
    std::vector< TVector3 > m_rotatedPoints;  ///< Transformed vector of points, centred at (0, 0, 0) with the primary direction along the z-axis
    double                  m_fitWindow;      ///< Fit window size ( in cm )
    TVector3                m_endpoint1;      ///< Position of the first end point
    TVector3                m_enddir1;        ///< Direction of the first end point    
    TVector3                m_endpoint2;      ///< Position of the second end point
    TVector3                m_enddir2;        ///< Direction of the second end point    
    std::vector< TVector3 > m_fitStart;       ///< Vector of fitted start points to each window of the sliding fit
    std::vector< TVector3 > m_fitCentre;      ///< Vector of fitted central points to each window of the sliding fit
    std::vector< TVector3 > m_fitEnd;         ///< Vector of fitted end points to each window of the sliding fit
    std::vector< TVector3 > m_fitDir;         ///< Vector of fitted directions to each window of the sliding fit
    std::vector< double >   m_fitSSR;         ///< Vector of sum of squared residuals for each window of the sliding fit
    double                  m_curvature;      ///< An average measure of amount of angular change over each fit window

    // Private member functions

    /**
     *  @brief  PerformLinearFit
     *
     *  @param  points              input vector of 3D positions on which to perform the fit
     *  @param  start               output 3D position at the start of the fit
     *  @param  centre              output 3D position at the centre of the fit
     *  @param  end                 output 3D position at the end of the fit
     *  @param  direction           output direction unit vector at the centre of the fit
     *  @param  sumSquareResiduals  output sum of squared residuals
     *  @param  curvature           output curvature of the fit
     */ 
    void PerformLinearFit( std::vector< TVector3 > & points, TVector3 & start, TVector3 & centre, TVector3 & end, TVector3 & direction, double & sumSquareResiduals ) const;
  };
};


crtagging::MyCosmicPFParticleTagger::MyCosmicPFParticleTagger(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  this->reconfigure( p );

  // Get the geometry information about the detector
  auto const * geo = lar::providerFrom<geo::Geometry>();

  // Get the detector extent
  m_TPCExtentX = 2 * geo->DetHalfWidth();
  m_TPCExtentY = 2 * geo->DetHalfHeight();
  m_TPCExtentZ =     geo->DetLength();

  // Get the central detector coordinate (using the MicroBooNE coordinate system)
  m_TPCCentralX = m_TPCExtentX / 2;
  m_TPCCentralY = 0;
  m_TPCCentralZ = m_TPCExtentZ / 2;
 
  // Calculate the positions of the detector faces
  m_TPCFaceMinX = m_TPCCentralX - m_TPCExtentX / 2;
  m_TPCFaceMaxX = m_TPCCentralX + m_TPCExtentX / 2;
  m_TPCFaceMinY = m_TPCCentralY - m_TPCExtentY / 2;
  m_TPCFaceMaxY = m_TPCCentralY + m_TPCExtentY / 2;
  m_TPCFaceMinZ = m_TPCCentralZ - m_TPCExtentZ / 2;
  m_TPCFaceMaxZ = m_TPCCentralZ + m_TPCExtentZ / 2;

  // Call appropriate produces<>() functions here.
  produces< std::vector<anab::CosmicTag>>();
  produces< art::Assns<anab::CosmicTag,   recob::Track>>();
  produces< art::Assns<recob::PFParticle, anab::CosmicTag>>();
}

void crtagging::MyCosmicPFParticleTagger::produce(art::Event & e)
{
  // Implementation of required member function here.

  // Instansiate the output
  // ----------------------------------------------------------------------------------------
  
  // A vector of the cosmic tags produced. Each tag consists of two end points, a cosmic score, and a tag id
  std::unique_ptr< std::vector< anab::CosmicTag > >                  cosmicTagTrackVector(       new std::vector<anab::CosmicTag>                  );

  // Associations between these tags and the reconstructed Tracks / PFParticles
  std::unique_ptr< art::Assns<anab::CosmicTag  , recob::Track    > > assnOutCosmicTagTrack(      new art::Assns<anab::CosmicTag,   recob::Track   >);
  std::unique_ptr< art::Assns<recob::PFParticle, anab::CosmicTag > > assnOutCosmicTagPFParticle( new art::Assns<recob::PFParticle, anab::CosmicTag>);

  // ----------------------------------------------------------------------------------------

  // Get the input collections
  // ----------------------------------------------------------------------------------------

  // Get the PFParticle handle
  art::Handle<std::vector<recob::PFParticle> > pfParticleHandle;
  e.getByLabel( m_PFParticleModuleLabel, pfParticleHandle);
  
  // Get the PFP --> SpacePoint associations
  art::FindManyP<recob::SpacePoint>  PFPToSpacePointListAssoc( pfParticleHandle, e, m_PFParticleModuleLabel );
 
  // Get the list of tracks - only required to match the functionality of the hit removal module 
  art::FindManyP<recob::Track> pfPartToTrackAssns(pfParticleHandle, e, m_TrackModuleLabel);
  
  // ----------------------------------------------------------------------------------------

  // Mapping between primary PFP ids and if out of time
  std::map< size_t, bool > pfpOutOfTimeMap;

  // Mapping between primary PFP ids and if they are top->bottom
  std::map< size_t, bool > pfpTopToBottomMap;
  
  // Mapping between primary PFP ids and their angle with the vertical
  std::map< size_t, double > pfpCosThetaMap;

  // Mapping between primary PFP ids and their curvature
  std::map< size_t, double > pfpCurvatureMap;

  // Mapping between PFP ids and their end points
  std::map< size_t, TVector3 > pfpToEndpoint1Map;
  std::map< size_t, TVector3 > pfpToEndpoint2Map;
   
  // Loop over all PFPs
  for( size_t pfPartIdx = 0; pfPartIdx != pfParticleHandle->size(); pfPartIdx++) {

    
    // Get the required objects from the handle
    // ----------------
    // Get a pointer to this PFP
    art::Ptr<recob::PFParticle> pPFP(pfParticleHandle, pfPartIdx);
 
    // Get a vector of pointers to the associated space points
    std::vector< art::Ptr< recob::SpacePoint > > spacePointList = PFPToSpacePointListAssoc.at( pfPartIdx );


    // In order to consider the PFP, we insist that it has a sufficient number of space points associated
    if ( spacePointList.size() < m_MinHitsToConsider ) continue;
    // ----------------

 
    // Now use the sliding fit class to get a geometrical description of the spacepoint List 
    // ----------------
    // Small fit window to get accurate endpoint
    SlidingFit fit ( spacePointList, 3.0 );
    TVector3 end1    = fit.GetEndpoint1();
    TVector3 end2    = fit.GetEndpoint2();

    pfpToEndpoint1Map.insert( std::make_pair( pPFP->Self(), end1 ) );
    pfpToEndpoint2Map.insert( std::make_pair( pPFP->Self(), end2 ) );

    // Larger fit window to get accurate end direction
    SlidingFit fitDir ( spacePointList, 10.0 );
    TVector3 endDir1 = -fitDir.GetEndDir1();
    TVector3 endDir2 =  fitDir.GetEndDir2();
    // ----------------
   
    // Now check if the PFP is out of time 
    // ----------------
    // Check if out of time
    bool outOfTime = this->IsOutOfTime( end1, end2, m_OutOfTimeThresh );
  
    // Get the id of the primary PFP from which this PFP originates
    size_t primaryId = this->GetPrimaryId( pfParticleHandle, pPFP->Self() );
   
    // If this PFP is out of time, mark the associated primary as out of time too, 
    // so it and all of its daughters can be tagged at a later stage
    if ( pfpOutOfTimeMap.count( primaryId ) == 0 ) {
      pfpOutOfTimeMap.insert( std::make_pair( primaryId, outOfTime ) ); 
    }
    else{
      pfpOutOfTimeMap[ primaryId ] = outOfTime || pfpOutOfTimeMap[ primaryId ];
    }
    // ----------------
 

    // From now on, consider only primary PFPs. 
    // Check their association to detector faces, and geometrical properties 
    // and those most likely to be cosmic ray muons
    // ----------------
    // Now we will only consider a PFPs geometry if they are primary
    if ( ! pPFP->IsPrimary() ) continue;

    // Associate PFP to detector faces
    std::vector< FaceType > faceList1;
    this->GetAssociatedFaces( end1, endDir1, m_FaceAssociationThresh, faceList1 );

    std::vector< FaceType > faceList2;
    this->GetAssociatedFaces( end2, endDir2, m_FaceAssociationThresh, faceList2 );
    
    // Check if the PFPs travels from the top to the bottom of the detector
    bool topToBottom = this->CheckFaces( faceList1, end1, faceList2, end2, TOP, BOTTOM);
    pfpTopToBottomMap.insert( std::make_pair( primaryId, topToBottom ) ); 

    // Check if the PFP is associated with the top
    bool topAssociated = ( this->CheckFace( faceList1, TOP ) || this->CheckFace( faceList2, TOP ) );

    // Check if the PFP is strongly vertically pointing
    double costheta = this->GetCosTheta( end1, endDir1, end2, endDir2 );
    pfpCosThetaMap.insert( std::make_pair( primaryId, costheta ) );

    // Get the curvature of the fit 
    double curvature = fit.GetCurvature();
    pfpCurvatureMap.insert( std::make_pair( primaryId, curvature) );
 
    std::cout << "PFP " << pPFP->Self() << "  -  " << spacePointList.size() << "  -  " << pPFP->IsPrimary() << std::endl;
    std::cout << "  Top?      : " << topAssociated << std::endl;
    std::cout << "  costheta  : " << costheta      << std::endl;
    std::cout << "  curvature : " << curvature     << std::endl;
    std::cout << std::endl;

    // ----------------
  }

  // Now we have mappings between primary PFPs and reasons for them to be tagged,
  // actually go through and do the tagging!
  // ----------------

  // Now loop over all PFPs again and make the ouput associations
  for(size_t pfPartIdx = 0; pfPartIdx != pfParticleHandle->size(); pfPartIdx++) {
    
    // Define the cosmic tag to use as a placeholder
    anab::CosmicTagID_t tag = anab::CosmicTagID_t::kNotTagged;
    float             score = 1.0;

    // Get the ID of the primary associated to this PFP
    art::Ptr<recob::PFParticle> pPFP(pfParticleHandle, pfPartIdx);
    size_t primaryId = this->GetPrimaryId( pfParticleHandle, pPFP->Self() );

    // Check if a map entry exists for that primary 
    if  ( pfpOutOfTimeMap.count( primaryId ) != 0 ) {
  
      // Check if out of time 
      if ( pfpOutOfTimeMap[ primaryId ] ) {
        tag = anab::CosmicTagID_t::kOutsideDrift_Partial;
      }
      // Check if top to bottom
      else if ( pfpTopToBottomMap[ primaryId ] ) {
        tag = anab::CosmicTagID_t::kGeometry_YY;
      }
      // Straight and downwards
      else if ( pfpCosThetaMap[ primaryId ] > 0.6 && pfpCurvatureMap[ primaryId ] < 0.05 ) {
        tag = anab::CosmicTagID_t::kGeometry_Y;
      }
    }

    // Add the cosmic tag to the vector
    std::vector<float> endPt1, endPt2;
    this->GetEndpointVector( pfpToEndpoint1Map, pPFP->Self(), endPt1 );
    this->GetEndpointVector( pfpToEndpoint2Map, pPFP->Self(), endPt2 );

    // Add this tag to the ouput
    cosmicTagTrackVector->emplace_back( endPt1, endPt2, score, tag);

    std::vector<art::Ptr<recob::Track> > trackVec = pfPartToTrackAssns.at(pfPartIdx);
    util::CreateAssn(*this, e, *cosmicTagTrackVector, trackVec, *assnOutCosmicTagTrack );

    util::CreateAssn(*this, e, *cosmicTagTrackVector, pPFP, *assnOutCosmicTagPFParticle);

  }
  
  e.put( std::move(cosmicTagTrackVector)      );
  e.put( std::move(assnOutCosmicTagTrack)     );
  e.put( std::move(assnOutCosmicTagPFParticle));

  // ----------------
  
}

// ------------------------------------------------------------------------------------------------------------------------------

size_t crtagging::MyCosmicPFParticleTagger::GetPrimaryId( art::Handle<std::vector<recob::PFParticle> > & pfParticleHandle, size_t thisId ) const
{
  // Get this PFP 
  art::Ptr<recob::PFParticle> pPFP(pfParticleHandle, thisId);
  
  // If primary, return ID 
  if ( pPFP->IsPrimary() ) return pPFP->Self();
  
  // Otherwise move to the parent
  return ( this->GetPrimaryId( pfParticleHandle, pPFP->Parent() ) );
}

// ------------------------------------------------------------------------------------------------------------------------------

void crtagging::MyCosmicPFParticleTagger::GetAssociatedFaces( TVector3 & endpoint, TVector3 & enddir, double thresh, std::vector< FaceType > & faceList ) const 
{
  // Check associations with each face
  if ( ( endpoint + enddir * thresh ).Y() > m_TPCFaceMaxY ) faceList.push_back( TOP        );
  if ( ( endpoint + enddir * thresh ).Y() < m_TPCFaceMinY ) faceList.push_back( BOTTOM     );
  if ( ( endpoint + enddir * thresh ).Z() > m_TPCFaceMaxZ ) faceList.push_back( DOWNSTREAM );
  if ( ( endpoint + enddir * thresh ).Z() < m_TPCFaceMinZ ) faceList.push_back( UPSTREAM   );
}

// ------------------------------------------------------------------------------------------------------------------------------

bool crtagging::MyCosmicPFParticleTagger::CheckFace( std::vector< FaceType > & faceList, FaceType type ) const
{
  return ( std::find( faceList.begin(), faceList.end(), type ) != faceList.end() );
}

// ------------------------------------------------------------------------------------------------------------------------------

bool crtagging::MyCosmicPFParticleTagger::CheckFaces( std::vector< FaceType > & faceList1, TVector3 & endpoint1, std::vector< FaceType > & faceList2, TVector3 & endpoint2, FaceType typeUpper, FaceType typeLower ) const
{
  // Find the uppper and lower endpoint
  std::vector< FaceType > faceListUpper = endpoint1.Y() > endpoint2.Y() ? faceList1 : faceList2;
  std::vector< FaceType > faceListLower = endpoint1.Y() > endpoint2.Y() ? faceList2 : faceList1;

  return ( this->CheckFace( faceListUpper, typeUpper ) && this->CheckFace( faceListLower, typeLower ) );
}

// ------------------------------------------------------------------------------------------------------------------------------

bool crtagging::MyCosmicPFParticleTagger::IsOutOfTime( TVector3 & endpoint1, TVector3 & endpoint2, double thresh ) const 
{
  return ( endpoint1.X() < m_TPCFaceMinX - thresh || 
           endpoint1.X() > m_TPCFaceMaxX + thresh ||
           endpoint2.X() < m_TPCFaceMinX - thresh ||
           endpoint2.X() > m_TPCFaceMaxX + thresh );
}

// ------------------------------------------------------------------------------------------------------------------------------

double crtagging::MyCosmicPFParticleTagger::GetCosTheta( TVector3 & endpoint1, TVector3 & endDir1, TVector3 & endpoint2, TVector3 & endDir2) const
{
  return ( endpoint1.Y() > endpoint2.Y() ? endDir1.Y() : endDir2.Y() );
}

// ------------------------------------------------------------------------------------------------------------------------------

void crtagging::MyCosmicPFParticleTagger::GetEndpointVector( std::map< size_t, TVector3 > & pfpToEndpointMap, size_t thisId, std::vector<float> & outputVect ) const
{
  if ( pfpToEndpointMap.count( thisId ) != 0 ) {
    TVector3 thisEndpoint = pfpToEndpointMap[ thisId ];
    outputVect.push_back( thisEndpoint.X() );
    outputVect.push_back( thisEndpoint.Y() );
    outputVect.push_back( thisEndpoint.Z() );
  }
  else{
    outputVect.push_back( std::numeric_limits<float>::max() );
    outputVect.push_back( std::numeric_limits<float>::max() );
    outputVect.push_back( std::numeric_limits<float>::max() );
  }
}

// ------------------------------------------------------------------------------------------------------------------------------

void crtagging::MyCosmicPFParticleTagger::beginJob()
{
  // Implementation of optional member function here.
}

void crtagging::MyCosmicPFParticleTagger::endJob()
{
  // Implementation of optional member function here.
}

void crtagging::MyCosmicPFParticleTagger::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  
  m_PFParticleModuleLabel = p.get< std::string >("PFParticleModuleLabel");
  m_TrackModuleLabel      = p.get< std::string >("TrackModuleLabel");
  m_MinHitsToConsider     = p.get< unsigned int >("MinHitsToConsider", 15);
  m_FaceAssociationThresh = p.get< double >("FaceAssociationThreshold", 20);
  m_OutOfTimeThresh       = p.get< double >("OutOfTimeThreshold", 5);
}

// -----------------------------------------------------------------------------------------
// Implementation of the SlidingFit class
// -----------------------------------------------------------------------------------------

crtagging::MyCosmicPFParticleTagger::SlidingFit::SlidingFit( std::vector< art::Ptr< recob::SpacePoint > > const & spacePointVect, double fitWindow ) {

  // Load the spacepoint into a vector of TVector3 objects
  for ( art::Ptr< recob::SpacePoint > const & point : spacePointVect ) {
    TVector3 p( point->XYZ()[0], point->XYZ()[1], point->XYZ()[2] );
    m_points.push_back( p );
  }

  // Perform a linear fit every space point to get the axis
  TVector3 axisStart;
  TVector3 axisCentre;
  TVector3 axisEnd;
  TVector3 axisDirection;
  double   axisSSR;
  this->PerformLinearFit( m_points, axisStart, axisCentre, axisEnd, axisDirection, axisSSR );

  // Find the angle of rotation such that the axis direction lies in the YZ plane
  TVector3 axisXY( axisDirection.X(), axisDirection.Y(), 0 );
  TVector3 yDir( 0, 1, 0 );
  double thetaZ = axisXY.Angle(yDir);

  TVector3 newAxisDir( axisDirection );
  newAxisDir.RotateZ( thetaZ );

  // Find the angle of rotation new axis lies along the Z-direction
  TVector3 zDir( 0, 0, 1 );
  double thetaX = newAxisDir.Angle( zDir );
  
  // Produce a new set of points, centered around 0,0,0 that point along the z-axis
  for ( TVector3 & point : m_points ) {
    TVector3 p( point );
    p -= axisCentre;
    p.RotateZ( thetaZ );
    p.RotateX( thetaX );
    m_rotatedPoints.push_back( p );   
  }

  // Sort by new z-position
  std::sort( m_rotatedPoints.begin( ), m_rotatedPoints.end( ), [ ]( const TVector3 & lhs, const TVector3 & rhs ) {
    return lhs.Z() < rhs.Z();
  });

  // Get the straight line length
  double minZ = m_rotatedPoints[ 0 ].Z();
  double maxZ = m_rotatedPoints[ m_rotatedPoints.size() - 1 ].Z(); 
  double L = maxZ - minZ;

  // Find the number of windows that will fit into this length, and the corresponding window length
  int nFitInLength = std::ceil( L / fitWindow );
  m_fitWindow = L / double( nFitInLength );
  
  // Find the number of points to sample
  int nSamples = 2*nFitInLength + 1;

  // Now break up the points into these windows
  unsigned int iLower = 0;
  unsigned int iUpper = 0;

  m_curvature = 0;
  int nSamplesUsed = 0;

  for ( int w=0; w < nSamples; w++ ) {
    double windowCentre = minZ + ( double (w) / double ( nSamples-1 ) ) * L;
    double windowLower  = windowCentre - m_fitWindow / 2;
    double windowUpper  = windowCentre + m_fitWindow / 2;

    while ( ( m_rotatedPoints[iLower].Z() < windowLower ) && ( iLower < m_rotatedPoints.size()-1 ) ) iLower++;
    iUpper = iLower;
    while ( ( m_rotatedPoints[iUpper].Z() < windowUpper ) && ( iUpper < m_rotatedPoints.size()-1 ) ) iUpper++;
    if ( iUpper != m_rotatedPoints.size() - 1 )  iUpper -= 1;
    
    std::vector< TVector3 > pointsInWindow; 
    for ( unsigned int i=iLower; i <= iUpper; i++ ) {
      // Rotate the point back to the world frame
      TVector3 p = m_rotatedPoints[i];
      p.RotateX( -thetaX );
      p.RotateZ( -thetaZ );
      p += axisCentre;
      pointsInWindow.push_back( p );
    }
    
    // Perform a linear fit on this window ( if it contains sufficient points )
    if ( pointsInWindow.size() > 1 ) {
      TVector3 winStart;
      TVector3 winCentre;
      TVector3 winEnd;
      TVector3 winDirection;
      double   winSSR;
      this->PerformLinearFit( pointsInWindow, winStart, winCentre, winEnd, winDirection, winSSR );
  
      // There is a degenerecy in the direction ( +- PI ), so we want to ensure
      // we always choose the same option
      double a = winDirection.Dot( axisDirection );
      winDirection *= ( a > 0 ) - ( a < 0 );
    
      m_fitStart.push_back( winStart ); 
      m_fitCentre.push_back( winCentre ); 
      m_fitEnd.push_back( winEnd ); 
      m_fitDir.push_back( winDirection ); 
      m_fitSSR.push_back( winSSR ); 

      m_curvature += (axisDirection - winDirection).Mag();
      nSamplesUsed++;
    }
  }

  m_curvature /= double ( nSamplesUsed );

  // Now set the end-points
  m_endpoint1 = m_fitStart[0];
  m_enddir1   = m_fitDir[0];
  m_endpoint2 = m_fitEnd[ m_fitEnd.size() - 1 ];
  m_enddir2   = m_fitDir[ m_fitDir.size() - 1 ];
}
    
void crtagging::MyCosmicPFParticleTagger::SlidingFit::PerformLinearFit( std::vector< TVector3 > & points, TVector3 & start, TVector3 & centre, TVector3 & end, TVector3 & direction, double & sumSquareResiduals ) const {

  // Perform PCA in 3D to find the direction
  Eigen::MatrixXd X( points.size(), 3 ); 

  // Fill the matrix with the variable
  int row = 0;
  for ( TVector3 & point : points ) {
    X( row, 0 ) = point.X();
    X( row, 1 ) = point.Y();
    X( row, 2 ) = point.Z();
    row++;
  }

  auto mean = X.colwise().mean();
  Eigen::MatrixXd centered = X.rowwise() - mean;
  Eigen::MatrixXd cov = centered.adjoint() * centered;

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);

  centre.SetX( mean(0) );
  centre.SetY( mean(1) );
  centre.SetZ( mean(2) );

  direction.SetX( eig.eigenvectors().col(2)(0) );
  direction.SetY( eig.eigenvectors().col(2)(1) );
  direction.SetZ( eig.eigenvectors().col(2)(2) );

  sumSquareResiduals = points.size() * eig.eigenvalues()(2);

  // Now find the start and end positions along the fit
  double maxD = - std::numeric_limits<double>::max();
  double minD = + std::numeric_limits<double>::max();
  for ( TVector3 & point : points ) {
    double D = (point - centre).Dot(direction);
    if ( D > maxD ) {
      maxD  = D;
      end   = centre + direction * D;  
    }
    if ( D < minD ) {
      minD  = D;
      start = centre + direction * D;  
    }
  }
  
}


DEFINE_ART_MODULE(crtagging::MyCosmicPFParticleTagger)
