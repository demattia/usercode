#include "DataFormats/Common/interface/Wrapper.h"
#include "AnalysisExamples/AnalysisObjects/interface/AnalyzedTrack.h"
#include "AnalysisExamples/AnalysisObjects/interface/AnalyzedCluster.h"
#include <vector>
#include <map>

namespace {
  namespace {
    anaobj::AnalyzedTrackCollection tkc1;
    edm::Wrapper<anaobj::AnalyzedTrackCollection> wtkc1;
    anaobj::AnalyzedTrackCollection::const_iterator cittkc1;
    edm::Wrapper<anaobj::AnalyzedTrackCollection::const_iterator> wcittkc1;
    anaobj::AnalyzedTrackCollection::iterator ittkc1;
    edm::Wrapper<anaobj::AnalyzedTrackCollection::iterator> wittkc1;
    anaobj::AnalyzedTrackRef rtk1;
    edm::Wrapper<anaobj::AnalyzedTrackRef> wrtk1;
    anaobj::AnalyzedTrackRefProd rptk1;
    edm::Wrapper<anaobj::AnalyzedTrackRefProd> wrptk1;
    anaobj::AnalyzedTrackRefVector rvtk1;
    edm::Wrapper<anaobj::AnalyzedTrackRefVector> wvtk1;
    anaobj::AnalyzedTrackRefVector::const_iterator cittk1;
    edm::Wrapper<anaobj::AnalyzedTrackRefVector::const_iterator> wcittk1;
    anaobj::AnalyzedTrackRefVector::iterator ittk1;
    edm::Wrapper<anaobj::AnalyzedTrackRefVector::iterator> wittk1;
    std::auto_ptr<anaobj::AnalyzedTrackRef> aptk1;

    anaobj::AnalyzedClusterCollection cluc1;
    edm::Wrapper<anaobj::AnalyzedClusterCollection> wcluc1;
    anaobj::AnalyzedClusterCollection::const_iterator citcluc1;
    edm::Wrapper<anaobj::AnalyzedClusterCollection::const_iterator> wcitcluc1;
    anaobj::AnalyzedClusterCollection::iterator itcluc1;
    edm::Wrapper<anaobj::AnalyzedClusterCollection::iterator> witcluc1;
    anaobj::AnalyzedClusterRef rc1;
    edm::Wrapper<anaobj::AnalyzedClusterRef> wrc1;
    anaobj::AnalyzedClusterRefProd rpc1;
    edm::Wrapper<anaobj::AnalyzedClusterRefProd> wrpc1;
    anaobj::AnalyzedClusterRefVector rvc1;
    edm::Wrapper<anaobj::AnalyzedClusterRefVector> wvc1;
    anaobj::AnalyzedClusterRefVector::const_iterator citc1;
    edm::Wrapper<anaobj::AnalyzedClusterRefVector::const_iterator> wcitc1;
    anaobj::AnalyzedClusterRefVector::iterator itc1;
    edm::Wrapper<anaobj::AnalyzedClusterRefVector::const_iterator> witc1;
    std::auto_ptr<anaobj::AnalyzedClusterRef> apc1;

    std::pair<int, int> pair_int;
    edm::Wrapper<std::pair<int, int> > wpair_int;
    std::map<int, int> map_int;
    edm::Wrapper<std::map<int, int> > wmap_int;
    std::map<int, int>::iterator map_int_iter;
    std::map<int, int>::const_iterator map_int_const_iter;

    std::pair<const int, int> cpair_int;
    std::map<const int, int> cmap_int;
    edm::Wrapper<std::map<const int, int> > cwmap_int;
    std::map<const int, int>::iterator cmap_int_iter;
    std::map<const int, int>::const_iterator cmap_int_const_iter;

    std::pair<int, float> pair_float;
    edm::Wrapper<std::pair<int, float> > wpair_float;
    std::map<int, float> map_float;
    edm::Wrapper<std::map<int, float> > wmap_float;
    std::map<int, float>::iterator map_float_iter;
    std::map<int, float>::const_iterator map_float_const_iter;

    std::pair<const int, float> cpair_float;
    std::map<const int, float> cmap_float;
    edm::Wrapper<std::map<const int, float> > cwmap_float;
    std::map<const int, float>::iterator cmap_float_iter;
    std::map<const int, float>::const_iterator cmap_float_const_iter;
  }
}
