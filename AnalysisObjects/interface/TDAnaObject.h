#ifndef TDANAOBJECT_H
#define TDANAOBJECT_H

namespace anaobj {

  /**
   * Stores all the variables produced in TDAna.
   * Saved in the event by TDAnaProducer.
   * ATTENTION: the genreflex seems to fail if there are too many parameters passed
   * to the constructor. Thus pass only the first 24+1 variables to it and set the
   * remaining with the set methods.
   *
   * Author M. De Mattia - 29/4/2008
   *
   */

  class TDAnaObject {
  public:
    /// Default constructor, only needed for classes.h
    TDAnaObject() {
      // 25 variables
      c8_ = 0.;
      m8_ = 0.;
      c6_ = 0.;
      m6_ = 0.;
      met_ = 0.;
      metsig_ = 0.;
      corrSumEt_ = 0.;
      goodHt_ = 0.;
      hbestcomb_ = 0.;
      chi2_ = 0.;
      mbbnohmax_ = 0.;
      dpbbnohmax_ = 0.;
      sumhed4_ = 0.;
      sumhed6_ = 0.;
      dpmin_ = 0.;
      dp1st_ = 0.;
      dp2nd_ = 0.;
      m_others_ = 0.;
      et6_ = 0.;
      scprodthbest_ = 0.;
      thdetabest_ = 0.;
      m5_ = 0.;
      tbestcomb_ = 0.;
      wbestcomb_ = 0.;
      ttms1_ = 0.;
      // This is another version
      metsignew_ = 0.;
      // Additional variables
      goodIc5Jets_ = 0;
      uncorrHt_ = 0.;
      corrHt_ = 0.;
      goodHt2_ = 0.;
      uncorrSumEt_ = 0.;
      goodSumEt_ = 0.;
      uncorrMEtSig_ = 0.;
      corrMEtSig_ = 0.;
      m45best_ = 0.;
      chi2ext_ = 0.;
      sumet_ = 0.;
      metphi_ = 0.;
      dp12_ = 0.;
      dpbb_ = 0.;
      m45bestall_ = 0.;
      chi2extall_ = 0.;
      dpbball_ = 0.;
      mwmin_ = 0.;
      drpairbestall_ = 0.;
      m3a_ = 0.;
      mwa_ = 0.;
      ttms2_ = 0.;
      ttms3_ = 0.;
      sumhpd4_ = 0.;
      sumhpd6_ = 0.;
      dp3rd_ = 0.;
      ptot_ = 0.;
      ptotE2_ = 0.;
      l1Trigger_ = false;
    }

    TDAnaObject( const double & C8, const double & M8, const double & C6, const double & M6,
                 const double & MET, const double & METSIG, const double & CORRSUMET, const double & GOODHT,
                 const double & HBESTCOMB, const double & CHI2, const double & MBBNOHMAX, const double & DPBBNOHMAX,
                 const double & SUMHED4, const double & SUMHED6, const double & DPMIN, const double & DP1ST,
                 const double & DP2ND, const double & M_OTHERS, const double & ET6, const double & SCPRODTHBEST,
                 const double & THDETABEST, const double & M5, const double & TBESTCOMB, const double & WBESTCOMB,
                 const double & TTMS1, const double & METSIGNEW )
//                  const int GOODIC5JETS, const double & UNCORRHT,
//                  const double & CORRHT, const double & GOODHT2, const double & UNCORRSUMET, const double & GOODSUMET,
//                  const double & UNCORRMETSIG, const double & CORRMETSIG, const double & M45BEST, const double & CHI2EXT )
//                  const double & SUMET, const double & METPHI, const double & DP12, const double & DPBB,
//                  const double & M45BESTALL, const double & CHI2EXTALL, const double & DPBBALL, const double & MWMIN,
//                  const double & DRPAIRBESTALL, const double & M3A, const double & MWA, const double & TTMS2,
//                  const double & TTMS3, const double & SUMHPD4, const double & SUMHPD6, const double & DP3RD,
//                  const double & PTOT, const double & PTOTE2 )

//     TDAnaObject( const double & C8, const double & M8, const double & C6, const double & M6,
//                  const double & MET, const double & METSIG, const double & CORRSUMET, const double & GOODHT,
//                  const double & HBESTCOMB, const double & CHI2, const double & MBBNOHMAX, const double & DPBBNOHMAX,
//                  const double & SUMHED4, const double & SUMHED6, const double & DPMIN, const double & DP1ST,
//                  const double & DP2ND, const double & M_OTHERS, const double & ET6, const double & SCPRODTHBEST,
//                  const double & THDETABEST, const double & M5, const double & TBESTCOMB, const double & WBESTCOMB,
//                  const double & TTMS1, const double & METSIGNEW, const int GOODIC5JETS, const double & UNCORRHT,
//                  const double & CORRHT, const double & GOODHT2, const double & UNCORRSUMET, const double & GOODSUMET,
//                  const double & UNCORRMETSIG, const double & CORRMETSIG, const double & M45BEST, const double & CHI2EXT, 
//                  const double & SUMET, const double & METPHI, const double & DP12, const double & DPBB,
//                  const double & M45BESTALL, const double & CHI2EXTALL, const double & DPBBALL, const double & MWMIN,
//                  const double & DRPAIRBESTALL, const double & M3A, const double & MWA, const double & TTMS2,
//                  const double & TTMS3, const double & SUMHPD4, const double & SUMHPD6, const double & DP3RD,
//                  const double & PTOT, const double & PTOTE2 )
      : c8_(C8),
        m8_(M8),
        c6_(C6),
        m6_(M6),
        met_(MET),
        metsig_(METSIG),
        corrSumEt_(CORRSUMET),
        goodHt_(GOODHT),
        hbestcomb_(HBESTCOMB),
        chi2_(CHI2),
        mbbnohmax_(MBBNOHMAX),
        dpbbnohmax_(DPBBNOHMAX),
        sumhed4_(SUMHED4),
        sumhed6_(SUMHED6),
        dpmin_(DPMIN),
        dp1st_(DP1ST),
        dp2nd_(DP2ND),
        m_others_(M_OTHERS),
        et6_(ET6),
        scprodthbest_(SCPRODTHBEST),
        thdetabest_(THDETABEST),
        m5_(M5),
        tbestcomb_(TBESTCOMB),
        wbestcomb_(WBESTCOMB),
        ttms1_(TTMS1),
        metsignew_(METSIGNEW),
//         goodIc5Jets_(GOODIC5JETS),
//         uncorrHt_(UNCORRHT),
//         corrHt_(CORRHT),
//         goodHt2_(GOODHT2),
//         uncorrSumEt_(UNCORRSUMET),
//         goodSumEt_(GOODSUMET),
//         uncorrMEtSig_(UNCORRMETSIG),
//         corrMEtSig_(CORRMETSIG),
//         m45best_(M45BEST),
//         chi2ext_(CHI2EXT),
//         sumet_(SUMET),
//         metphi_(METPHI),
//         dp12_(DP12),
//         dpbb_(DPBB),
//         m45bestall_(M45BESTALL),
//         chi2extall_(CHI2EXTALL),
//         dpbball_(DPBBALL),
//         mwmin_(MWMIN),
//         drpairbestall_(DRPAIRBESTALL),
//         m3a_(M3A),
//         mwa_(MWA),
//         ttms2_(TTMS2),
//         ttms3_(TTMS3),
//         sumhpd4_(SUMHPD4),
//         sumhpd6_(SUMHPD6),
//         dp3rd_(DP3RD),
//         ptot_(PTOT),
//         ptotE2_(PTOTE2)



        goodIc5Jets_(0),
        uncorrHt_(0.),
        corrHt_(0.),
        goodHt2_(0.),
        uncorrSumEt_(0.),
        goodSumEt_(0.),
        uncorrMEtSig_(0.),
        corrMEtSig_(0.),
        m45best_(0.),
        chi2ext_(0.),
        sumet_(0.),
        metphi_(0.),
        dp12_(0.),
        dpbb_(0.),
        m45bestall_(0.),
        chi2extall_(0.),
        dpbball_(0.),
        mwmin_(0.),
        drpairbestall_(0.),
        m3a_(0.),
        mwa_(0.),
        ttms2_(0.),
        ttms3_(0.),
        sumhpd4_(0.),
        sumhpd6_(0.),
        dp3rd_(0.),
        ptot_(0.),
        ptotE2_(0.),
        l1Trigger_(false)
    {
    }
    // Getter methods
    double c8() const { return c8_; }
    double m8() const { return m8_; }
    double c6() const { return c6_; }
    double m6() const { return m6_; }
    double met() const { return met_; }
    double metsig() const { return metsig_; }
    double corrSumEt() const { return corrSumEt_; }
    double goodHt() const { return goodHt_; }
    double hbestcomb() const { return hbestcomb_; }
    double chi2() const { return chi2_; }
    double mbbnohmax() const { return mbbnohmax_; }
    double dpbbnohmax() const { return dpbbnohmax_; }
    double sumhed4() const { return sumhed4_; }
    double sumhed6() const { return sumhed6_; }
    double dpmin() const { return dpmin_; }
    double dp1st() const { return dp1st_; }
    double dp2nd() const { return dp2nd_; }
    double m_others() const { return m_others_; }
    double et6() const { return et6_; }
    double scprodthbest() const { return scprodthbest_; }
    double thdetabest() const { return thdetabest_; }
    double m5() const { return m5_; }
    double tbestcomb() const { return tbestcomb_; }
    double wbestcomb() const { return wbestcomb_; }
    double ttms1() const { return ttms1_; }
    // This is another version
    double metsignew() const { return metsignew_; }
    // Additional variables
    int goodIc5Jets() const { return goodIc5Jets_; }
    double uncorrHt() const { return uncorrHt_; }
    double corrHt() const { return corrHt_; }
    double goodHt2() const { return goodHt2_; }
    double uncorrSumEt() const { return uncorrSumEt_; }
    double goodSumEt() const { return goodSumEt_; }
    double uncorrMEtSig() const { return uncorrMEtSig_; }
    double corrMEtSig() const { return corrMEtSig_; }
    double m45best() const { return m45best_; }
    double chi2ext() const { return chi2ext_; }
    double sumet() const { return sumet_; }
    double metphi() const { return metphi_; }
    double dp12() const { return dp12_; }
    double dpbb() const { return dpbb_; }
    double m45bestall() const { return m45bestall_; }
    double chi2extall() const { return chi2extall_; }
    double dpbball() const { return dpbball_; }
    double mwmin() const { return mwmin_; }
    double drpairbestall() const { return drpairbestall_; }
    double m3a() const { return m3a_; }
    double mwa() const { return mwa_; }
    double ttms2() const { return ttms2_; }
    double ttms3() const { return ttms3_; }
    double sumhpd4() const { return sumhpd4_; }
    double sumhpd6() const { return sumhpd6_; }
    double dp3rd() const { return dp3rd_; }
    double ptot() const { return ptot_; }
    double ptotE2() const { return ptotE2_; }
    bool l1Trigger() const { return l1Trigger_; }

    // Setter methods
    void setC8(const double & C8) { c8_ = C8; }
    void setM8(const double & M8) { m8_ = M8; }
    void setC6(const double & C6) { c6_ = C6; }
    void setM6(const double & M6) { m6_ = M6; }
    void setMet(const double & MET) { met_ = MET; }
    void setMetsig(const double & METSIG) { metsig_ = METSIG; }
    void setCorrSumEt(const double & CORRSUMET) { corrSumEt_ = CORRSUMET; }
    void setGoodHt(const double & GOODHT) { goodHt_ = GOODHT; }
    void setHbestcomb(const double & HBESTCOMB) { hbestcomb_ = HBESTCOMB; }
    void setChi2(const double & CHI2) { chi2_ = CHI2; }
    void setMbbnohmax(const double & MBBNOHMAX) { mbbnohmax_ = MBBNOHMAX; }
    void setDpbbnohmax(const double & DPBBNOHMAX) { dpbbnohmax_ = DPBBNOHMAX; }
    void setSumhed4(const double & SUMHED4) { sumhed4_ = SUMHED4; }
    void setSumhed6(const double & SUMHED6) { sumhed6_ = SUMHED6; }
    void setDpmin(const double & DPMIN) { dpmin_ = DPMIN; }
    void setDp1st(const double & DP1ST) { dp1st_ = DP1ST; }
    void setDp2nd(const double & DP2ND) { dp2nd_ = DP2ND; }
    void setM_others(const double & M_OTHERS) { m_others_ = M_OTHERS; }
    void setEt6(const double & ET6) { et6_ = ET6; }
    void setScprodthbest(const double & SCPRODTHBEST) { scprodthbest_ = SCPRODTHBEST; }
    void setThdetabest(const double & THDETABEST) { thdetabest_ = THDETABEST; }
    void setM5(const double & M5) { m5_ = M5; }
    void setTbestcomb(const double & TBESTCOMB) { tbestcomb_ = TBESTCOMB; }
    void setWbestcomb(const double & WBESTCOMB) { wbestcomb_ = WBESTCOMB; }
    void setTtms1(const double & TTMS1) { ttms1_ = TTMS1; }
    // This is another version
    void setMetsignew(const double & METSIGNEW) { metsignew_ = METSIGNEW; }
    // Additional variables
    void setGoodIc5Jets(const int GOODIC5JETS) { goodIc5Jets_ = GOODIC5JETS; }
    void setUncorrHt(const double & UNCORRHT) { uncorrHt_ = UNCORRHT; }
    void setCorrHt(const double & CORRHT) { corrHt_ = CORRHT; }
    void setGoodHt2(const double & GOODHT2) { goodHt2_ = GOODHT2; }
    void setUncorrSumEt(const double & UNCORRSUMET) { uncorrSumEt_ = UNCORRSUMET; }
    void setGoodSumEt(const double & GOODSUMET) { goodSumEt_ = GOODSUMET; }
    void setUncorrMEtSig(const double & UNCORRMETSIG) { uncorrMEtSig_ = UNCORRMETSIG; }
    void setCorrMEtSig(const double & CORRMETSIG) { corrMEtSig_ = CORRMETSIG; }
    void setM45best(const double & M45BEST) { m45best_ = M45BEST; }
    void setChi2ext(const double & CHI2EXT) { chi2ext_ = CHI2EXT; }
    void setSumet(const double & SUMET) { sumet_ = SUMET; }
    void setMetphi(const double & METPHI) { metphi_ = METPHI; }
    void setDp12(const double & DP12) { dp12_ = DP12; }
    void setDpbb(const double & DPBB) { dpbb_ = DPBB; }
    void setM45bestall(const double & M45BESTALL) { m45bestall_ = M45BESTALL; }
    void setChi2extall(const double & CHI2EXTALL) { chi2extall_ = CHI2EXTALL; }
    void setDpbball(const double & DPBBALL) { dpbball_ = DPBBALL; }
    void setMwmin(const double & MWMIN) { mwmin_ = MWMIN; }
    void setDrpairbestall(const double & DRPAIRBESTALL) { drpairbestall_ = DRPAIRBESTALL; }
    void setM3a(const double & M3A) { m3a_ = M3A; }
    void setMwa(const double & MWA) { mwa_ = MWA; }
    void setTtms2(const double & TTMS2) { ttms2_ = TTMS2; }
    void setTtms3(const double & TTMS3) { ttms3_ = TTMS3; }
    void setSumhpd4(const double & SUMHPD4) { sumhpd4_ = SUMHPD4; }
    void setSumhpd6(const double & SUMHPD6) { sumhpd6_ = SUMHPD6; }
    void setDp3rd(const double & DP3RD) { dp3rd_ = DP3RD; }
    void setPtot(const double & PTOT) { ptot_ = PTOT; }
    void setPtotE2(const double & PTOTE2) { ptotE2_ = PTOTE2; }
    void setL1Trigger(const bool L1TRIGGER) { l1Trigger_ = L1TRIGGER; }
  protected:
    // 25 variables
    double c8_;
    double m8_;
    double c6_;
    double m6_;
    double met_;
    double metsig_;
    double corrSumEt_;
    double goodHt_;
    double hbestcomb_;
    double chi2_;
    double mbbnohmax_;
    double dpbbnohmax_;
    double sumhed4_;
    double sumhed6_;
    double dpmin_;
    double dp1st_;
    double dp2nd_;
    double m_others_;
    double et6_;
    double scprodthbest_;
    double thdetabest_;
    double m5_;
    double tbestcomb_;
    double wbestcomb_;
    double ttms1_;
    // This is another version
    double metsignew_;
    // Additional variables
    int goodIc5Jets_;
    double uncorrHt_;
    double corrHt_;
    double goodHt2_;
    double uncorrSumEt_;
    double goodSumEt_;
    double uncorrMEtSig_;
    double corrMEtSig_;
    double m45best_;
    double chi2ext_;
    double sumet_;
    double metphi_;
    double dp12_;
    double dpbb_;
    double m45bestall_;
    double chi2extall_;
    double dpbball_;
    //           for ( int i=0; i<iJ && i<NHSJ; i++ ) {
    //             HEDSSSW_->Fill(JHET[i],PTOTE2);
    //             HPDSSSW_->Fill(JHPT[i],PTOTE2);
    //           } 
    double mwmin_;
    double drpairbestall_;
    double m3a_;
    double mwa_;
    double ttms2_;
    double ttms3_;
    double sumhpd4_;
    double sumhpd6_;
    double dp3rd_;
    double ptot_;
    double ptotE2_;
    bool l1Trigger_;
  };

}

#endif //MCPARTICLE_H
