#ifndef _JET_H
#define _JET_H

#include "allvars.hpp"

/// \name Jets namespace 
//@{
namespace jets {

struct BriefJet {
    double     eta, phi, kt2, NN_dist;
    BriefJet * NN;
    int        _jets_index;
};

const double MaxRap = 1e5;
const double pseudojet_invalid_phi = -100.0;
const double pseudojet_invalid_rap = -1e200;

// PseudoJet class 
class PseudoJet {

private:
  double _px,_py,_pz,_E;
  mutable double _phi, _rap;
  double _kt2;
  void _set_rap_phi() const;
public:
  PseudoJet() : _px(0), _py(0), _pz(0), _E(0) {}
  PseudoJet(const double px, const double py, const double pz, const double E);
  ~PseudoJet() = default;
  double E()   const {return _E;}
  double px()  const {return _px;}
  double py()  const {return _py;}
  double pz()  const {return _pz;}
  double phi() const {_phi = _phi > pi ? _phi-twopi : _phi; return _phi;}
  double rap() const {return _rap;}
  double perp() const {return sqrt(_kt2);}
  double pseudorapidity() const;
  double eta() const {return pseudorapidity();}
  double m2() const {return (_E + _pz) * (_E - _pz) - _kt2;}

  double kt_distance(const PseudoJet & other) const;
  double plain_distance(const PseudoJet & other) const;
  inline double delta_R(const PseudoJet & other) const {return sqrt(plain_distance(other));}
  double delta_phi_to(const PseudoJet & other) const;
  double kt2() const {return _kt2;}
  double beam_distance() const {return _kt2;}
  std::valarray<double> four_mom() const;
  enum { X=0, Y=1, Z=2, T=3, NUM_COORDINATES=4, SIZE=NUM_COORDINATES };

  double operator () (int i) const ;
  inline double operator [] (int i) const { return (*this)(i); };

  void operator*=(double);
  void operator/=(double);
  void operator+=(const PseudoJet &);
  void operator-=(const PseudoJet &);

  PseudoJet operator+(const PseudoJet &);
  PseudoJet operator-(const PseudoJet &);
  PseudoJet operator*(const PseudoJet &);
  PseudoJet operator*(double);
  PseudoJet operator/(double);

  bool operator==(const PseudoJet &);
  bool operator!=(const PseudoJet & b) {return !(*this==b);}
  //bool operator==(const double val);
  //bool operator!=(const double & val) {return !(*this==val);}
};

inline double dot_product(const PseudoJet & a, const PseudoJet & b) {
  return a.E()*b.E() - a.px()*b.px() - a.py()*b.py() - a.pz()*b.pz();
}

// ClusterSequence class  
class ClusterSequence {
public:
  //temmplate <class J> std::vector<J> _jets;
  std::vector<PseudoJet> _jets;
  ClusterSequence () 
  { 
    // _jets;
    _R2=pi;
    _invR2 = 1.0/_R2;
  };
  ~ClusterSequence () = default;

  struct history_element{
    int parent1; /// index in _history where first parent of this jet
    int parent2; /// index in _history where second parent of this jet
    int child;   /// index in _history where the current jet is
		 /// recombined with another jet to form its child. It
		 /// is Invalid if this jet does not further
		 /// recombine.
    int jetp_index; /// index in the _jets vector where we will find the
    double dij;  /// the distance corresponding to the recombination
		 /// at this stage of the clustering.
    double max_dij_so_far; /// the largest recombination distance seen
			   /// so far in the clustering history.
  };

  template<class J> void simple_N2_cluster();
  double jet_scale_for_algorithm(const PseudoJet & jet) const ;
  template <class J> void bj_set_jetinfo(J * jeta, const int i) const;
  // template <class J> void bj_set_jetinfo(J &jeta, const int i) const;
  void do_ij_recombination_step(
      const int & jet_i, const int & jet_j, 
      const double & dij, int & newjet_k
      );
  void do_iB_recombination_step(const int & jet_i, const double & diB);

  template <class J> double bj_dist(const J * const jeta, const J * const jetb) const;
  template <class J> double bj_diJ(const J * const jeta) const;

  // template <class J> double bj_dist(const J &jeta, const J &jetb) const;
  // template <class J> double bj_diJ(const J &jeta) const;

  template <class J> void bj_set_NN_nocross(J * const jeta,
    J * const head, const J * const tail) const;
  template <class J> void bj_set_NN_crosscheck(J * const jeta,
    J * const head, const J * const tail) const;

  template <class J> void bj_set_NN_nocross(J &jeta,
    int ihead, int itail) const;
  template <class J> void bj_set_NN_crosscheck(J &jeta,
    int ihead, int itail) const;
  private: 
    double _R2, _invR2;
    // inline double _wrap_rapidity(double x) {
    //   return x;
      
    // };
    //inline double _wrap_phi(double x) {return x > twopi ? x-twopi : x; };
};

}
//@}

#endif