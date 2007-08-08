#include "DataFormats/Common/interface/Wrapper.h"
#include "AnalysisExamples/PixelJet/interface/PixelJet.h"
#include <vector>
#include <map>

namespace {
  namespace {
    PixelJetCollection pjc1;
    edm::Wrapper<PixelJetCollection> wpjc1;
    PixelJetCollection::const_iterator citpjc1;
    edm::Wrapper<PixelJetCollection::const_iterator> wcitpjc1;
    PixelJetCollection::iterator itpjc1;
    edm::Wrapper<PixelJetCollection::iterator> witpjc1;
    PixelJetRef rpj1;
    edm::Wrapper<PixelJetRef> wrpj1;
    PixelJetRefProd rppj1;
    edm::Wrapper<PixelJetRefProd> wrppj1;
    PixelJetRefVector rvpj1;
    edm::Wrapper<PixelJetRefVector> wvrpj1;
    PixelJetRefVector::const_iterator citpj1;
    edm::Wrapper<PixelJetRefVector::const_iterator> wcitpj1;
    PixelJetRefVector::iterator itpj1;
    edm::Wrapper<PixelJetRefVector::iterator> witpj1;
    std::auto_ptr<PixelJetRef> appj1;

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
