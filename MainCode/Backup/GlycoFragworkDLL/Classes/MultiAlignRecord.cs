using System;
using System.Collections.Generic;

namespace GlycoFragworkDLL.Classes
{
    public enum _ID_label {Verified=0 , Tentative , Unverified} ; 
    public class MultiAlignRecord
    {
        int _ID;
        int _Size; 
        double _Mass;
        double _NET;

        // A collection of peptide sequences and glycan compositions
        public string[] _CandidateProteinName;
        public string[] _CandidatePeptideSeq;
        public string[] _CandidateGlycanComposition;
        public string[] _CandidateNGlycoSite; 
        public double [] _CandidateGlycanMass;
        public double [] _CandidatePeptideMass;

        public GlycopeptideRecord[] _CandidateGlycopeptideRecords; 

        // Just the best matching one 
        string _BestMatchProteinName;
        string _BestMatchPeptideSeq;
        string _BestMatchGlycanComposition;
        string _BestMatchNGlycoSite; 
        double _BestMatchGlycanMass;
        double _BestMatchPeptideMass;
        bool _BestMatchFalseHit;
        double _BestMatchParentMz;
        double _BestMatchParentScanTime;
        string _BestMatchGlycanSequence; 
                
        
        public GlypID.Peaks.clsPeak[] _RepresentativeCIDPeaks;
        public GlypID.Peaks.clsPeak[] _RepresentativeHCDPeaks;
        public GlypID.Peaks.clsPeak[] _RepresentativeETDPeaks;

        public int _RepresentativeDatasetID_CID = -1;
        public int _RepresentativeDatasetID_HCD = -1;
        public int _RepresentativeDatasetID_ETD = -1;

        int _RepresentativeCIDScore = 0;
        double _RepresentativeHCDScore = 1.0;
        double _RepresentativeETDScore = 0.0; 
        double _RepresentativeCIDSequencingScore = 0 ; 

        bool _IsResolved = false ; 


        public GlypID.HornTransform.clsHornTransformResults _RepresentativeTransformResult; 

        // A collection of _umc_records is associated with a single multi align record
        // Each  umc record comes from one dataset        
        public List<UMCRecord> _AssociatedUMCRecords;

        public List<FragEvents> _ClusterRepFragEvents; 
        public List<string> _ClusterNames; 

        
      
        _ID_label _id_label = _ID_label.Unverified ; 

        public MultiAlignRecord()
        {
            _ID = -1;
            _Size = 0;
            _Mass = 0;
            _NET = 0;
           
            _AssociatedUMCRecords = new List<UMCRecord>();

            _BestMatchProteinName = "";
            _BestMatchPeptideSeq = "";
            _BestMatchGlycanComposition = "";
            _BestMatchGlycanSequence = "";
            _BestMatchNGlycoSite = "" ;
            _BestMatchGlycanSequence = "";
            _BestMatchPeptideMass = 0;           
            _BestMatchGlycanMass = 0;
            _BestMatchFalseHit = false;
            _BestMatchParentMz = 0;
            _BestMatchParentScanTime = 0; 

            _ClusterNames = new List<string>();
            _ClusterRepFragEvents = new List<FragEvents>(); 

        }

        public MultiAlignRecord(MultiAlignRecord m)
        {
            _ID = m._ID;
            _Size = m._Size;
            _Mass = m._Mass;
            _NET = m._NET;
            _AssociatedUMCRecords = new List<UMCRecord>();

            if ((m._CandidatePeptideSeq != null) && (m._CandidatePeptideSeq.Length > 0))
            {
                _CandidatePeptideSeq = m._CandidatePeptideSeq;
                _CandidateProteinName = m._CandidateProteinName;
                _CandidatePeptideMass = m._CandidatePeptideMass;
                _CandidateGlycanMass = m._CandidateGlycanMass;
                _CandidateGlycanComposition = m._CandidateGlycanComposition;
                _CandidateNGlycoSite = m._CandidateNGlycoSite; 
               
            }
            AllocateNumberDatasets(m._AssociatedUMCRecords.Count);
            for (int i = 0; i < m._AssociatedUMCRecords.Count; i++)
            {
                _AssociatedUMCRecords.Add(m._AssociatedUMCRecords[i]);
            }
        }

        public bool HasUniformChargeStates()
        {
            bool has_unif_cs = true;
            if (_AssociatedUMCRecords.Count > 1)
            {               
                short cs = _AssociatedUMCRecords[0].UMCRepCS; 
                for (int i = 1; i < _AssociatedUMCRecords.Count; i++)
                {
                    if (_AssociatedUMCRecords[i].UMCRepCS != cs)
                    {
                        has_unif_cs = false;
                        break;
                    }
                }
            }

            return has_unif_cs; 
                
        }
       

        public bool HasUniformIDsAcrossClusters()
        {
            bool has_unif_ids = true;

            if (_ClusterRepFragEvents.Count > 1)
            {
                int begin_search = 0;
                while (begin_search < _ClusterRepFragEvents.Count)
                {
                    if ((_ClusterRepFragEvents[begin_search].ETDScore > 0) && (_ClusterRepFragEvents[begin_search].FDR < 0.05) && (_ClusterRepFragEvents[begin_search].GP_Record.IsDecoy == false))
                        break;
                    begin_search++;
                }

                if (begin_search < _ClusterRepFragEvents.Count)
                {
                    string protein = _ClusterRepFragEvents[begin_search].GP_Record.Sequence.proteinName;
                    string site = _ClusterRepFragEvents[begin_search].GP_Record.Sequence.nGlycoSite;
                    string glycan = _ClusterRepFragEvents[begin_search].GP_Record.Glycan.composition;
                    for (int i = begin_search + 1; i < _ClusterRepFragEvents.Count; i++)
                    {
                        // This way only ids which have been scored will get allocated
                        if ((_ClusterRepFragEvents[i].ETDScore > 0) && (_ClusterRepFragEvents[i].FDR < 0.05) && (_ClusterRepFragEvents[begin_search].GP_Record.IsDecoy == false))
                        {

                            if ((_ClusterRepFragEvents[i].GP_Record.Sequence.proteinName == protein) && (_ClusterRepFragEvents[i].GP_Record.Glycan.composition == glycan))
                            {
                                int i1 = site.IndexOf(_ClusterRepFragEvents[i].GP_Record.Sequence.nGlycoSite);
                                int i2 = _ClusterRepFragEvents[i].GP_Record.Sequence.nGlycoSite.IndexOf(site);
                                if ((i1 == -1) || (i2 == -1)) // This is very unlikely
                                {
                                    has_unif_ids = false;
                                    break;
                                }
                            }
                            else
                            {
                                has_unif_ids = false;
                                break;
                            }
                        }
                    }
                }
                //else
                  //  has_unif_ids = false; // All were decoy hits, should ignore
            }
           

            return has_unif_ids; 
        }

        public void SetIdsBasedOnMassAndGlycanType(double ppm_tolerance, bool look_in_clustered)
        {            
            int best_index = -1;
            double min_ppm = ppm_tolerance ;
            double min_hcd = 1; 
            int min_hcd_index = -1;
            double max_cid = 0;
            int max_cid_index = -1; 
            GlycoFragworkDLL.Utils.Utilities _utils = new GlycoFragworkDLL.Utils.Utilities();
            FragEvents e = new FragEvents() ; 
            if (look_in_clustered)
            {
                for (int i = 0; i < _ClusterRepFragEvents.Count; i++)
                {
                    if (_ClusterRepFragEvents[i].HCDScore < min_hcd)
                    {
                        min_hcd = _ClusterRepFragEvents[i].HCDScore;
                        min_hcd_index = i;
                    }
                    if (_ClusterRepFragEvents[i].CIDScore > max_cid)
                    {
                        max_cid = _ClusterRepFragEvents[i].CIDScore;
                        max_cid_index = i;
                    }
                }
                if (min_hcd < 1)
                     _RepresentativeHCDScore = min_hcd;
                if (max_cid > 0)
                    _RepresentativeCIDScore = (int) max_cid;
                _RepresentativeETDScore = 0;
                _RepresentativeCIDSequencingScore = 0; 
                e = new FragEvents() ;           
                e = _ClusterRepFragEvents[min_hcd_index];
                
            }
            else
            {
                int min_hcd_umc_index = -1; 
                int max_cid_umc_index = -1;
                for (int i = 0; i < _AssociatedUMCRecords.Count; i++)
                {
                    for (int j = 0; j < _AssociatedUMCRecords[i]._AssociatedFragEvents.Count; j++)
                    {
                        if (_AssociatedUMCRecords[i]._AssociatedFragEvents[j].HCDScore < min_hcd)
                        {
                            min_hcd = _AssociatedUMCRecords[i]._AssociatedFragEvents[j].HCDScore;
                            min_hcd_index = j;
                            min_hcd_umc_index = i;
                        }
                        if (_AssociatedUMCRecords[i]._AssociatedFragEvents[j].CIDScore > max_cid)
                        {
                            max_cid = _AssociatedUMCRecords[i]._AssociatedFragEvents[j].CIDScore;
                            max_cid_umc_index = i;
                            max_cid_index = j;
                        }
                    }
                }
                if (min_hcd < 1)
                    _RepresentativeHCDScore = min_hcd;
                if (max_cid > 0)
                    _RepresentativeCIDScore = (int)max_cid;
                _RepresentativeETDScore = 0;
                _RepresentativeCIDSequencingScore = 0;
                e = new FragEvents();
                e = _AssociatedUMCRecords[min_hcd_umc_index]._AssociatedFragEvents[min_hcd_index];        
            }


           
            bool check_sialylated = false ;
            if ((look_in_clustered) && (e.GP_Record.Glycan.numNeuAc > 0))
            {
                check_sialylated = true;
            }
            else
            {
                if ((e.GlycanType == GlypID.enmGlycanType.CS) || e.GlycanType == GlypID.enmGlycanType.HY)
                    check_sialylated = true;
            }

            for (int k = 0; k < _CandidateGlycopeptideRecords.Length; k++)
            {
                bool is_sialylated = false ;
                if (_CandidateGlycopeptideRecords[k].Glycan.numNeuAc>0)
                    is_sialylated = true ; 
                if( check_sialylated == is_sialylated)
                {
                    double ppm_diff = _utils.CalculateDelMassPPM(_CandidateGlycopeptideRecords[k].GP_Mono_Mass, _Mass); 
                    if (ppm_diff < min_ppm)
                    {
                        best_index = k;
                        min_ppm = ppm_diff;
                    }
                }
            }
            if (best_index > -1)
            {
                _BestMatchProteinName = _CandidateGlycopeptideRecords[best_index].Sequence.proteinName;
                _BestMatchPeptideSeq = _CandidateGlycopeptideRecords[best_index].Sequence.sequence;
                _BestMatchNGlycoSite = _CandidateGlycopeptideRecords[best_index].Sequence.nGlycoSite;
                _BestMatchGlycanComposition = _CandidateGlycopeptideRecords[best_index].Glycan.composition;
                _BestMatchGlycanMass = _CandidateGlycopeptideRecords[best_index].GlycanMonoMass;
                _BestMatchPeptideMass = _CandidateGlycopeptideRecords[best_index].SequenceMonoMass;
                _BestMatchFalseHit = _CandidateGlycopeptideRecords[best_index].IsDecoy;
                _BestMatchParentScanTime = e.ParentScanTime;
                _BestMatchParentMz = e.ParentMz;
                

                if (e.CIDPeaks.Length > 1)
                {
                    _RepresentativeCIDPeaks = e.CIDPeaks; 
                }
                if (e.HCDPeaks.Length > 1)
                {
                    _RepresentativeHCDPeaks = e.HCDPeaks; 
                }                
                
            }
        }

        public void SetScoresBasedOnLowestHCD()
        {
            double min_hcd = 1;
            int min_hcd_index = -1;

            for (int i = 0; i < _ClusterRepFragEvents.Count; i++)
            {
                if (_ClusterRepFragEvents[i].HCDScore < min_hcd)
                {
                    min_hcd = _ClusterRepFragEvents[i].HCDScore;
                    min_hcd_index = i;
                }
            }

            if (min_hcd < 1)
            {
                IDLabel = _ID_label.Unverified;
                _RepresentativeCIDScore = (int)_ClusterRepFragEvents[min_hcd_index].CIDScore;
                _RepresentativeETDScore = _ClusterRepFragEvents[min_hcd_index].ETDScore;
                _RepresentativeHCDScore = _ClusterRepFragEvents[min_hcd_index].HCDScore;
                _RepresentativeCIDSequencingScore = _ClusterRepFragEvents[min_hcd_index].CIDSequencingScore;
                _BestMatchParentScanTime = _ClusterRepFragEvents[min_hcd_index].ParentScanTime;
                _BestMatchParentMz = _ClusterRepFragEvents[min_hcd_index].ParentMz;

                if (_ClusterRepFragEvents[min_hcd_index].CIDPeaks.Length > 1)
                {
                    _RepresentativeCIDPeaks = _ClusterRepFragEvents[min_hcd_index].CIDPeaks;
                }
                if (_ClusterRepFragEvents[min_hcd_index].HCDPeaks.Length > 1)
                {
                    _RepresentativeHCDPeaks = _ClusterRepFragEvents[min_hcd_index].HCDPeaks;
                }
                if (_ClusterRepFragEvents[min_hcd_index].ETDPeaks.Length > 1)
                {
                    _RepresentativeETDPeaks = _ClusterRepFragEvents[min_hcd_index].ETDPeaks;
                } 
            }
        }

        public void SetIdsBasedOnHighestETD()
        {
            double max_etd = 0;
            int max_etd_index = -1;

            for (int i = 0; i < _ClusterRepFragEvents.Count; i++)
            {
                if (_ClusterRepFragEvents[i].ETDScore > max_etd)
                {
                    max_etd = _ClusterRepFragEvents[i].ETDScore;
                    max_etd_index = i;                          
                }
            }

            if (max_etd > 0)
            {
                _BestMatchProteinName = _ClusterRepFragEvents[max_etd_index].GP_Record.Sequence.proteinName;
                _BestMatchPeptideSeq = _ClusterRepFragEvents[max_etd_index].GP_Record.Sequence.sequence;
                _BestMatchNGlycoSite = _ClusterRepFragEvents[max_etd_index].GP_Record.Sequence.nGlycoSite;
                _BestMatchGlycanComposition = _ClusterRepFragEvents[max_etd_index].GP_Record.Glycan.composition;
                _BestMatchGlycanMass = _ClusterRepFragEvents[max_etd_index].GP_Record.GlycanMonoMass;
                _BestMatchPeptideMass = _ClusterRepFragEvents[max_etd_index].GP_Record.SequenceMonoMass;
                _BestMatchFalseHit = _ClusterRepFragEvents[max_etd_index].GP_Record.IsDecoy;
                _BestMatchGlycanSequence = _ClusterRepFragEvents[max_etd_index].GP_Record.GlycanSequence; 
                if (_ClusterRepFragEvents[max_etd_index].FDR < 0.05)
                    IDLabel = _ID_label.Verified;
                else
                    IDLabel = _ID_label.Tentative;
                _RepresentativeCIDScore = (int) _ClusterRepFragEvents[max_etd_index].CIDScore;
                _RepresentativeETDScore = _ClusterRepFragEvents[max_etd_index].ETDScore;
                _RepresentativeHCDScore = _ClusterRepFragEvents[max_etd_index].HCDScore;
                _RepresentativeCIDSequencingScore = _ClusterRepFragEvents[max_etd_index].CIDSequencingScore;
                _BestMatchParentScanTime = _ClusterRepFragEvents[max_etd_index].ParentScanTime;
                _BestMatchParentMz = _ClusterRepFragEvents[max_etd_index].ParentMz;

                if (_ClusterRepFragEvents[max_etd_index].CIDPeaks.Length > 1)
                {
                    _RepresentativeCIDPeaks = _ClusterRepFragEvents[max_etd_index].CIDPeaks;
                }
                if (_ClusterRepFragEvents[max_etd_index].HCDPeaks.Length >1)
                {
                    _RepresentativeHCDPeaks = _ClusterRepFragEvents[max_etd_index].HCDPeaks; 
                }
                if (_ClusterRepFragEvents[max_etd_index].ETDPeaks.Length >1)
                {
                    _RepresentativeETDPeaks = _ClusterRepFragEvents[max_etd_index].ETDPeaks; 
                }
                  
            }
        }

        //void SetIdsBasedOnLowestHCD()
        //{


        public short MaxChargeStateObserved()
        {
            short max_cs = 0;

            for (int i = 0; i < _AssociatedUMCRecords.Count; i++)
            {
                for (int j = 0; j < _AssociatedUMCRecords[i]._AssociatedFragEvents.Count; j++)
                {
                    if (_AssociatedUMCRecords[i]._AssociatedFragEvents[j].TransformResult.mshort_cs > max_cs)
                        max_cs = _AssociatedUMCRecords[i]._AssociatedFragEvents[j].TransformResult.mshort_cs;
                }
            }
            return max_cs; 
        }

        public short MinChargeStateObserved()
        {
            short min_cs = short.MaxValue;
            for (int i = 0; i < _AssociatedUMCRecords.Count; i++)
            {
                for (int j = 0; j < _AssociatedUMCRecords[i]._AssociatedFragEvents.Count; j++)
                {
                    if (_AssociatedUMCRecords[i]._AssociatedFragEvents[j].TransformResult.mshort_cs < min_cs)
                        min_cs = _AssociatedUMCRecords[i]._AssociatedFragEvents[j].TransformResult.mshort_cs;
                }
            }
            return min_cs; 
        }


        public void AllocateNumberDatasets(int i)
        {
            if (_AssociatedUMCRecords.Count > 0)
            {
                _AssociatedUMCRecords.Clear();
            }
            _AssociatedUMCRecords = new List<UMCRecord>(i) ; 
        }

       /* public double MinHCDScore(bool among_clustered_only)
        {
            double min_score = 1; 
            if (among_clustered_only)

            /*for (int i = 0; i < _AssociatedUMCRecords.Count; i++)
            {
                if (_AssociatedUMCRecords[i].HCDScore < min_score)
                    min_score = _AssociatedUMCRecords[i].HCDScore; 
            }
            return min_score; 
        }

       


        public double MaxCIDScore()
        {
            double max_score = 0;
            for (int i = 0; i < _AssociatedUMCRecords.Count; i++)
            {
                if (_AssociatedUMCRecords[i].CIDScore > max_score)
                    max_score = _AssociatedUMCRecords[i].CIDScore;
            }
            return max_score; 
        }*/

       /* public bool HasProteinDs()
        {
            bool has_id = false;
            for (int i = 0; i < _AssociatedUMCRecords.Count; i++)
            {
                 if (_AssociatedUMCRecords[i].ProteinName != "")

            return has_id; 
        }*/
        public int ID
        {
            get { return _ID; }
            set { _ID = value; }
        }
        public int Size
        {
            get { return _Size; }
            set { _Size = value; }
        }
        public double Mass
        {
            get { return _Mass; }
            set { _Mass = value; }
        }
        public double NET
        {
            get { return _NET; }
            set { _NET= value; }
        }

        public double BestMatchPeptideMass
        {
            get { return _BestMatchPeptideMass; }
            set { _BestMatchPeptideMass = value; }
        }
        public double BestMatchGlycanMass
        {
            get { return _BestMatchGlycanMass; }
            set { _BestMatchGlycanMass = value; }
        }

        public string BestMatchProteinName
        {
            get { return _BestMatchProteinName; }
            set { _BestMatchProteinName = value; }
        }
        public string BestMatchPeptideSeq
        {
            get { return _BestMatchPeptideSeq; }
            set { _BestMatchPeptideSeq = value; }
        }
        public string BestMatchGlycanComposition
        {
            get { return _BestMatchGlycanComposition; }
            set { _BestMatchGlycanComposition = value; }
        }
        public string BestMatchGlycanSequence
        {
            get { return _BestMatchGlycanSequence; }
            set { _BestMatchGlycanSequence = value; }
        }
        public string BestMatchNGlycoSite
        {
            get { return _BestMatchNGlycoSite; }
            set { _BestMatchNGlycoSite = value; }
        }

        public bool BestMatchFalseHit
        {
            get { return _BestMatchFalseHit; }
            set { _BestMatchFalseHit = value; }
        }
        public bool IsResolved
        {
            get { return _IsResolved ; }
            set { _IsResolved = value ; }
        }

        public int RepCIDScore
        {
            get { return _RepresentativeCIDScore ; }
            set { _RepresentativeCIDScore = value; }
        }
        public double RepCIDSequencingScore
        {
            get { return _RepresentativeCIDSequencingScore; }
            set { _RepresentativeCIDSequencingScore = value; }
        }
        public double RepHCDScore
        {
            get { return _RepresentativeHCDScore; }
            set { _RepresentativeHCDScore = value; }
        }

        public double RepETDScore
        {
            get { return _RepresentativeETDScore; }
            set { _RepresentativeETDScore = value; }
        }

        public double BestMatchParentMz
        {
            get { return _BestMatchParentMz; }
            set { _BestMatchParentMz = value; }
        }

        public double BestMatchParentScanTime
        {
            get { return _BestMatchParentScanTime; }
            set { _BestMatchParentScanTime = value; }
        }
        public _ID_label IDLabel
        {
            get { return _id_label ;}
            set { _id_label = value ; }

        }
    }
}
