
#include "jet.hpp"

namespace jets {

PseudoJet::PseudoJet(const double px_in, const double py_in, const double pz_in, const double E_in) {
  _E  = E_in ;
  _px = px_in;
  _py = py_in;
  _pz = pz_in;
  _kt2 = _px * _px + _py * _py;
  _set_rap_phi();
}
void PseudoJet::_set_rap_phi() const {
    // set _phi
    if (_kt2 == 0.0) { _phi = 0.0;}
    else { _phi = atan2(_py,_px); }
    if (_phi < 0.0) {_phi += twopi;}
    else if (_phi >= twopi) {_phi -= twopi;} 
  
    // set rapidity
    auto abspz = _pz > 0 ? _pz : -_pz; 

    if ( _E == abspz && _kt2 == 0) {
        double MaxRapHere = MaxRap + abspz;
        if (_pz >= 0.0) _rap = MaxRapHere;
        else _rap = -MaxRapHere;
    }
    else 
    {
        double effective_m2 = std::max(0.0,m2()); // force non tachyonic mass
        double E_plus_pz    = _E + abspz; // the safer of p+, p-
        _rap = 0.5*log((_kt2 + effective_m2)/(E_plus_pz*E_plus_pz));
        if (_pz > 0) _rap = - _rap;
    }
}

std::valarray<double> PseudoJet::four_mom() const {
  std::valarray<double> mom(4);
  mom[0] = _px;
  mom[1] = _py;
  mom[2] = _pz;
  mom[3] = _E ;
  return mom;
}
double PseudoJet::operator () (int i) const {
  switch(i) {
  case X:
    return _px;
  case Y:
    return _py;
  case Z:
    return _pz;
  case T:
    return _E;
  default:
    std::ostringstream err;
    err << "PseudoJet subscripting: bad index (" << i << ")";
    //throw Error(err.str());
  }
  return 0.;
}
double PseudoJet::pseudorapidity() const {
  if (_px == 0.0 && _py == 0.0) return MaxRap;
  if (_pz == 0.0) return 0.0;
  double theta = atan(perp()/_pz);
  if (theta < 0) theta += pi;
  return -log(tan(theta/2.0));
}
PseudoJet PseudoJet::operator+ (const PseudoJet & jet2) {
  return PseudoJet(_px+jet2._px,
		   _py +jet2._py,
		   _pz +jet2._pz,
		   _E +jet2._E);
}
PseudoJet PseudoJet::operator-(const PseudoJet & jet2) {
  return PseudoJet(_px-jet2._px,
		   _py-jet2._py,
		   _pz-jet2._pz,
		   _E -jet2._E);
}
PseudoJet PseudoJet::operator*(double coeff) {
  return PseudoJet(_px * coeff, _py * coeff, _pz * coeff, _E  * coeff);
}
PseudoJet PseudoJet::operator/ (double coeff) {
  return operator*(1.0/coeff);
}
void PseudoJet::operator*=(double coeff) {
  _px *= coeff;
  _py *= coeff;
  _pz *= coeff;
  _E  *= coeff;
  _kt2*= coeff*coeff;
}
void PseudoJet::operator/=(double coeff) {
  (*this) *= 1.0/coeff;
}
void PseudoJet::operator+=(const PseudoJet & other_jet) {
  _px += other_jet._px;
  _py += other_jet._py;
  _pz += other_jet._pz;
  _E  += other_jet._E ;
  _kt2 = _px * _px + _py * _py;
  _set_rap_phi();
}
void PseudoJet::operator-=(const PseudoJet & other_jet) {
  _px -= other_jet.px();
  _py -= other_jet.py();
  _pz -= other_jet.pz();
  _E  -= other_jet.E() ;
  _kt2 = _px * _px + _py * _py;
  _set_rap_phi();
}
bool PseudoJet::operator==(const PseudoJet & b) {
  if (_px != b.px()) return false;
  if (_py != b.py()) return false;
  if (_pz != b.pz()) return false;
  if (_E  != b.E() ) return false;
  return true;
}



double ClusterSequence::jet_scale_for_algorithm(const PseudoJet & jet) const 
{
    return jet.kt2();
//   if (_jet_algorithm == kt_algorithm)             {return jet.kt2();}
//   else if (_jet_algorithm == cambridge_algorithm) {return 1.0;}
//   else if (_jet_algorithm == antikt_algorithm) {
//     double kt2=jet.kt2();
//     return kt2 > 1e-300 ? 1.0/kt2 : 1e300;
//   } else if (_jet_algorithm == genkt_algorithm) {
//     double kt2 = jet.kt2();
//     double p   = jet_def().extra_param();
//     if (p <= 0 && kt2 < 1e-300) kt2 = 1e-300; // dodgy safety check
//     return pow(kt2, p);
//   } else if (_jet_algorithm == cambridge_for_passive_algorithm) {
//     double kt2 = jet.kt2();
//     double lim = _jet_def.extra_param();
//     if (kt2 < lim*lim && kt2 != 0.0) {
//       return 1.0/kt2;
//     } else {return 1.0;}
//   } else {throw Error("Unrecognised jet algorithm");}
}

template <class J> void ClusterSequence::bj_set_jetinfo(
  J * jetA, const int _jets_index) const 
{
    jetA->eta  = _jets[_jets_index].eta();
    jetA->phi  = _jets[_jets_index].phi();
    jetA->kt2  = jet_scale_for_algorithm(_jets[_jets_index]);
    jetA->_jets_index = _jets_index;
    jetA->NN_dist = _R2;
    //jetA->NN      = NULL;
    jetA->NN      = nullptr;
}

// template <class J> void ClusterSequence::bj_set_jetinfo(
//   J &jetA, const int _jets_index) const 
// {
//     jetA.eta  = _jets[_jets_index].eta();
//     jetA.phi  = _jets[_jets_index].phi();
//     jetA.kt2  = jet_scale_for_algorithm(_jets[_jets_index]);
//     jetA._jets_index = _jets_index;
//     jetA.NN_dist = _R2;
//     jetA.NN      = nullptr;
// }

// distance between two jets from their phi and eta 
//template <class J> inline double ClusterSequence::_bj_dist(
template <class J> double ClusterSequence::bj_dist(
  const J * const jetA, const J * const jetB) const {
  //double dphi = std::abs(jetA->phi - jetB->phi);
  // if (dphi > pi) {dphi = twopi - dphi;}
  // get absolute value and wrap 
  double dphi = jetA->phi - jetB->phi;
  dphi = (dphi < 0) ? -dphi : dphi;
  dphi = (dphi > pi) ? twopi - dphi : dphi;
  double deta = (jetA->eta - jetB->eta);
  return dphi*dphi + deta*deta;
}

// template <class J> double ClusterSequence::bj_dist(
//   const J &jetA, const J &jetB) const {
//   double dphi = jetA.phi - jetB.phi;
//   double deta = (jetA.eta - jetB.eta);
//   dphi = (dphi < 0) ? -dphi : dphi;
//   dphi = (dphi > pi) ? twopi - dphi : dphi;
//   return dphi*dphi + deta*deta;
// }

// distance is either jet's kt2 value * either the NN's kt2 value 
// times the NN_distance (which if object has no NN, is _R2)
//template <class J> inline double ClusterSequence::_bj_diJ(
template <class J> double ClusterSequence::bj_diJ(
  const J * const jet) const {
  double kt2 = jet->kt2;
  //if (jet->NN != NULL) {if (jet->NN->kt2 < kt2) {kt2 = jet->NN->kt2;}}
  if (jet->NN != NULL) kt2 = (jet->NN->kt2 < kt2) ? jet->NN->kt2 : kt2;
  return jet->NN_dist * kt2;
}

// template <class J> double ClusterSequence::bj_diJ(
//   const J &jet) const {
//   double kt2 = jet.kt2;
//   if (jet.NN != NULL) kt2 = (jet.NN->kt2 < kt2) ? jet.NN->kt2 : kt2;
//   return jet.NN_dist * kt2;
// }

// for a given jet input and head and tails of array
// iterate over all jets between head and tail to calcaute 
// the jet distance and also set NN, setting only the jet's NN 
// values (hence nocross). 
// Can be parallelised
//template <class J> inline void ClusterSequence::_bj_set_NN_nocross(
template <class J> void ClusterSequence::bj_set_NN_nocross(
  J * const jet, J * const head, const J * const tail) const 
{
  double NN_dist = _R2;
  J * NN  = NULL;
  if (head < jet) {
    for (J * jetB = head; jetB != jet; jetB++) {
      double dist = bj_dist(jet,jetB);
      if (dist < NN_dist) {
	      NN_dist = dist;
	      NN = jetB;
      }
    }
  }
  if (tail > jet) {
    for (J * jetB = jet+1; jetB != tail; jetB++) {
      double dist = bj_dist(jet,jetB);
      if (dist < NN_dist) {
	      NN_dist = dist;
	      NN = jetB;
      }
    }
  }
  jet->NN = NN;
  jet->NN_dist = NN_dist;
}

// determine the NN of a jet, given an head and an tail
// loops over all jets from head to tail and gets distance between 
// jet and all the other jets. if jet is minimum set 
//template <class J> inline void ClusterSequence::_bj_set_NN_crosscheck(
template <class J> void ClusterSequence::bj_set_NN_crosscheck(
  J * const jet, J * const head, const J * const tail) const 
{
  double NN_dist = _R2;
  J * NN  = NULL;
  for (J * jetB = head; jetB != tail; jetB++) {
    double dist = bj_dist(jet,jetB);
    if (dist < NN_dist) {
      NN_dist = dist;
      NN = jetB;
    }
    if (dist < jetB->NN_dist) {
      jetB->NN_dist = dist;
      jetB->NN = jet;
    }
  }
  jet->NN = NN;
  jet->NN_dist = NN_dist;
}

/// still working on adding in recombination 
void ClusterSequence::do_ij_recombination_step(
    const int & jet_i, const int & jet_j,
    const double & dij,
    int & newjet_k) 
{
    PseudoJet newjet(false);
    _jet_def.recombiner()->recombine(_jets[jet_i], _jets[jet_j], newjet);
    _jets.push_back(newjet);
    newjet_k = _jets.size()-1;
    int newstep_k = _history.size();
    _jets[newjet_k].set_cluster_hist_index(newstep_k);
    int hist_i = _jets[jet_i].cluster_hist_index();
    int hist_j = _jets[jet_j].cluster_hist_index();
    _add_step_to_history(newstep_k, std::min(hist_i, hist_j), std::max(hist_i,hist_j), newjet_k, dij);
}
void ClusterSequence::do_iB_recombination_step(
    const int & jet_i, const double & diB) 
{
  int newstep_k = _history.size();
  _add_step_to_history(newstep_k,_jets[jet_i].cluster_hist_index(),BeamJet, Invalid, diB);
}


// unit tests for ktjet clustering algorithm that scales as N^2
template<class BJ> void ClusterSequence::simple_N2_cluster()
{
    int nthreads = omp_get_max_threads();
    //  allocate BJ jet class array
    // use jetA and jetB pointers to point to it 
    // and set properties
    int n = _jets.size();
    std::vector<BJ> briefjets(n), briefjets_NN(n);
    BJ * head = &briefjets[0], * tail = &briefjets[n-1];
    BJ *jetA, *jetB;
    std::vector<double> diJ(n);
    n=0;
    // loop over the briefjets and initialise
    #pragma omp parallel for \
    schedule(static) \
    default(none) shared(briefjets, n) private(jetA) \
    num_threads(nthreads) \
    if (nthreads > 1) 
    for (auto i=0;i<n;i++) {
        jetA = &briefjets[i];
        bj_set_jetinfo(jetA, i);
    }

    // generate loop to produce a crosscheck and NN and NN distance 
    // for briefjet and other jets between head and itself. 
    // relistically since this will update not just briefjet[i] but 
    // all jets from head to briefjet[i] cannot be trivially parallelised 
    // as written 
    // HOWEVER, since need to calculate distances between all possible pairs anyway
    // makes sense to actually just separate the distance calculations out from 
    // determining the NN and NN distances 
    for (auto i=1;i<n;i++) bj_set_NN_crosscheck(&briefjets[i], head, &briefjets[i]);

    // calculate distances of jet based on kt2 and NN_dist. 
    #pragma omp parallel for \
    schedule(static) \
    default(none) shared(briefjets, diJ, n) private(jetA) \
    num_threads(nthreads) \
    if (nthreads > 1) 
    for (auto i=0;i<n;i++) {
        jetA = &briefjets[i];
        diJ[i] = bj_diJ(jetA);
    } 

    // PJE still need to figure out what this is for history starts at the end (or tail)
    int history_location = n-1;
    // one starts at the active tail and current head 
    // and the history location of the tail position in the briefjets array 
    // 1) find the minimum NN distance for all active jets 
    // 2) set the active jetA to the jet with the minimum NN distance and jetB to its NN 
    // 3) check jetB (NN of jetA)
    // 3.1) if the NN of jetA is NOT NULL then
    // 3.1.1) if the pointer jetA < jetB, swap the memory between the two so that jetB is 
    //        at the lower index and jetA is at the higher index
    // 3.1.2) recombine the jets and set jetB (lower index) info 
    // 3.2) if the NN is NULL then "recombine" jet with iB, making it a real jet
    // 4) now the minimum distance jet has either been recombined with another briefjet
    //    or has been fully processed and become a real jet with _do_iB_recombination_step
    //    decrease the active tail, and nominal size of the briefjet array since
    //    have removed a jet
    //    set diJ[jetA - head] = diJ[tail-head];
    //    PJE: I am very confused by above since *jetA = *tail. jetA should be 
    //    pointing to the lower index address of the jet with the minimum NN distance
    //    yet the current active tail jet's properties are now copied to the position of jetA
    //    so although *jetA == *tail, jetA != tail. Tail has also been decreased since a jet 
    //    has been processed. 
    //    so above is diJ[diJ_min_jet - 0] = diJ[iActiveTail - 0];
    // 5) loop over active jet array from 0 to iActiveTail and 
    // 5.1) if jetI pointer does not point to jetA pointer and jetB pointer 
    //      that is ijet != diJ_min_jet && != briefjet[diJ_min_jet].NN index
    //      calculate _bj_set_NN_nocross(jetI, head, tail) and update distance of 
    //      briefjet[ijetI] with the distance based on kt2 and distance to NN (or _R2 if no NN)
    // 5.2) if the current jetB is not NULL and jetI != jetB 
    //      then also compute the jetI, jetB distance and update NNs if necessary
    // 5.3) update the NN of jetI (=briefjet[ijetI]) if the NN points to the current 
    //      active tail to the jetA (=briefjet[diJ_min_jet] (don't know why)
    // 6) update kt2 *NN_dist of jetB stored in diJ

    // So issue is the updates to the NN list via _do_ij_recombination_step and possibly _do_iB_recombination_step that causes distances to be recalculate for new pairs 


    // go through into tail and head 
    while (tail != head) {
        /// 1)
        //@{
        double diJ_min = diJ[0];
        int diJ_min_jet = 0;
        for (int i = 1; i < n; i++) {
            if (diJ[i] < diJ_min) {diJ_min_jet = i; diJ_min  = diJ[i];}
        }
        //@}

        /// 2)
        //@{
        // increment the history location (which at the start, would point to 
        // a bad address (has value of n)
        history_location++;
        // select minimum jet
        jetA = & briefjets[diJ_min_jet];
        // PJE: what does NN do? 
        jetB = static_cast<BJ *>(jetA->NN);
        // PJE: why multiple temporary variable by _invR2? what is _invR2?
        diJ_min *= _invR2;
        //@}

        /// 3)
        //@{
        /// 3.1) 
        //@{
        // PJE: if jetB which is min jet's NN is NULL 
        // call _do_iB_recombination_setp
        // otherwise some some more preliminary calculations 
        if (jetB != NULL) {
            // 3.1.1 
            //@{}
            // PJE: I am confused as why would you swap at all? 
            // jetA is a jet within the BJ array
            // but jetB is is just jetA -> NN. Is this guaranteed to be 
            // in the array? 
            if (jetA < jetB) {std::swap(jetA,jetB);}
            //@}
            // 3.1.2
            //@{
            int nn; // new jet index
            // PJE: what does this do?
            do_ij_recombination_step(jetA->_jets_index, jetB->_jets_index, diJ_min, nn);
            // PJE: what does this do?
            bj_set_jetinfo(jetB, nn);
            //@}
        } 
        //@}
        /// 3.2)
        //@{
        else {
            do_iB_recombination_step(jetA->_jets_index, diJ_min);
        }
        //@}
        //@}

        /// 4)
        //@{
        // PJE: if jetB is not null then jetB has been updated
        // as has diJ_min updated (possibly), though diJ_min_jet has not been 
        // updated so that unlikely

        // PJE: why decrease n? this is the size of the array
        // PJE: move tail towards head
        tail--; n--;
        // PJE: why copy the data residing in tail to jetA? 
        // jetA currently points to the jet with min diJ
        *jetA = *tail;
        // PJE: now with a possibly updated diJ_min 

        // PJE: I am confused as to why pointers are used to 
        // to access the double array of diJ
        // would it not be better to use integers taht 
        // store the index of tail and head? Certainly clearer
        diJ[jetA - head] = diJ[tail-head];
        //@}

        /// 5)
        //@{
        // PJE: use jetI to iterate over a now 
        // reduced head/tail set
        for (BJ * jetI = head; jetI != tail; jetI++) 
        {
            /// 5.1 
            //@{
            // PJE: if jetI is jetA or jetB then do 
            // computation of NN_no_cross and update the diJ array
            if (jetI->NN == jetA || jetI->NN == jetB) 
            {
                bj_set_NN_nocross(jetI, head, tail);
                diJ[jetI-head] = bj_diJ(jetI); // update diJ
            }
            //@}

            /// 5.2 
            //@{
            // PJE: now if jetB is not null (jetB = briefjets[diJ_min_jet]->NN)
            // alter jetI and the diJ array 
            if (jetB != NULL && jetI != jetB) {
                double dist = bj_dist(jetI,jetB);
                if (dist < jetI->NN_dist) {
                    jetI->NN_dist = dist;
                    jetI->NN = jetB;
                    diJ[jetI-head] = bj_diJ(jetI); // update diJ...
                }
                // PJE: if dist is smaller update jetB NN list
                if (dist < jetB->NN_dist) {
                    jetB->NN_dist = dist;
                    jetB->NN      = jetI; 
                    // PJE: why is _bj_diJ of jetB NOT update?! NN_dist has changed 
                    // so surely it should update the diJ?
                    // ANSWER: this because jetB distance updated after for loop
                }
            }
            //@}

            // PJE: if jetI's NN list contains tail, move NN to jetA 
            // which is weird because it is tail, though it is a deep copy
            if (jetI->NN == tail) {jetI->NN = jetA;}
        }
        //@}
        // PJE: end of for loop 

        /// 6
        //@{
        // PJE: if jet B is not null, then update again the diJ array 
        if (jetB != NULL) {diJ[jetB-head] = bj_diJ(jetB);}
        //@}
    }
    // end of while loop 
    //@}

    // PJE: NOTES on while loop
    // so key is that loop does calculate lots of _bj_diJ(jetA,jetB) distance pairs
    // for almost all combinations of jetA and jetB. 
    // thus it might be useful to just calculate all pairs on GPU 

    // // free memory 
    // delete[] diJ;
    // delete[] briefjets;
}

template void ClusterSequence::simple_N2_cluster<BriefJet>();

}