#ifndef SUMMARY_H
#define SUMMARY_H

namespace anaobj {

  class Summary {
  public:
    /// Empty constructor, needed only for classes.h
    Summary() {}
    Summary( const unsigned int EVENTTYPE, const int JETNUM, const bool HTOBB, const bool SEMILEPTONIC ) {
      eventType_ = EVENTTYPE;
      jetNum_ = JETNUM;
      hTobb_ = HTOBB;
      semileptonic_ = SEMILEPTONIC;
    }
    unsigned int eventType() const { return eventType_; }
    int jetNum() const { return jetNum_; }
    bool hTobb() const { return hTobb_; }
    bool semileptonic() const { return semileptonic_; }
    void setEventType( const unsigned int EVENTTYPE ) { eventType_ = EVENTTYPE; }
    void setJetNum( const int JETNUM ) { jetNum_ = JETNUM; }
    void setHTobb( const bool HTOBB ) { hTobb_ = HTOBB; }
    void setSemileptonic( const bool SEMILEPTONIC) { semileptonic_ = SEMILEPTONIC; }
  protected:
    unsigned int eventType_;
    int jetNum_;
    bool hTobb_;
    bool semileptonic_;
  };

}

#endif //SUMMARY
