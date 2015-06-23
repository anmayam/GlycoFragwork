using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GlycoFragworkDLL.Tasks
{
  /*  public class MapScorer
    {
     
        GlypID.CIDScoring.clsCIDScoring _CIDScoring;
        GlypID.HCDScoring.clsHCDScoring _HCDScoring;
        GlypID.ETDScoring.clsETDScoring _ETDScoring;
       
        Utils.Utilities _utils;

        public string output_folder; 
        List <string> peptides;


       

        public MapScorer()
        {
            
            _CIDScoring = new GlypID.CIDScoring.clsCIDScoring();
            _HCDScoring = new GlypID.HCDScoring.clsHCDScoring();
            _ETDScoring = new GlypID.ETDScoring.clsETDScoring();
            _utils = new GlycoFragworkDLL.Utils.Utilities();
            peptides = new List<string>() ; 
            
        }

        public void SetOptions(GlycoFragworkDLL.Classes.Params scoreParams)
        {
            _HCDScoring.ScoringParameters  = scoreParams.ScoringParams  ;           
            _ETDScoring.ScoringParameters = scoreParams.ScoringParams;
            _CIDScoring.ScoringParameters = scoreParams.ScoringParams;
           /* _HCDScoring.TransformParameters = scoreParams.TransformParams;            
            _ETDScoring.TransformParameters = scoreParams.TransformParams;            
            _CIDScoring.TransformParameters = scoreParams.TransformParams; 
        }

        public void SetNglycopeptides(List<string> sequences)
        {
            peptides.Clear();
            foreach (string s in sequences)
            {
                peptides.Add(s);
            }

        }

        public void ScoreMap(ref Classes.MapRecord _glycoMap)
        {
           
            /*for (int i = 0; i < _glycoMap._AllMLNRecords.Count; i++)
            {
                Classes.MultiAlignRecord m = _glycoMap._AllMLNRecords[i];
                int num_records = m._AssociatedUMCRecords.Count;

               
                for (int j = 0; j < num_records; j++)
                {
                    Classes.UMCRecord u = m._AssociatedUMCRecords[j];

                    // Get out parent
                    if (u.DatasetID == m._RepresentativeDatasetID_CID)
                    {
                        GlypID.HornTransform.clsHornTransformResults[] TransformResults = new GlypID.HornTransform.clsHornTransformResults[1];
                        TransformResults[0] = u.TransformResult;
                        GlypID.Peaks.clsPeak[] CIDPeaks = m._RepresentativeCIDPeaks;

                        float[] mzs = new float[CIDPeaks.Length];
                        float[] intensities = new float[CIDPeaks.Length]; 
                        int index = 0 ;
                        foreach (GlypID.Peaks.clsPeak peak in CIDPeaks)
                        {
                            mzs[index] = Convert.ToSingle(CIDPeaks[index].mdbl_mz);
                            intensities[index] = Convert.ToSingle(CIDPeaks[index].mdbl_intensity);
                        }

                        GlypID.CIDScoring.clsCIDScoringScanResults[] CIDScoreResults = new GlypID.CIDScoring.clsCIDScoringScanResults[1];
                        bool found_score = _CIDScoring.ScoreCIDSpectra(ref CIDPeaks, ref mzs, ref intensities, ref TransformResults, ref CIDScoreResults);

                        // Do sequencing
                      /*  foreach (string s in peptides)
                        {
                            string PeptideSeq = s;
                            GlycoSeq.AminoAcidWeight AAMW = new GlycoSeq.AminoAcidWeight();
                            float PeptideMass = AAMW.GetMonoMW(PeptideSeq, true);
                            int y1cs = 2;
                            float y1mz = (PeptideMass + GlycoSeq.GlycanMass.GetMZ(GlycoSeq.Glycan.Types.HexNAc, 1) + (float)_utils._CC_MASS * y1cs) / y1cs;

                            double y1mtest = _utils.CalculateMz((PeptideMass + 203.007937), 2);



                            GlycoSeq.MSScan _msScan = new GlycoSeq.MSScan(CIDPeaks, mzs, intensities, Convert.ToSingle(TransformResults[0].mdbl_mz), Convert.ToInt32(TransformResults[0].mshort_cs), 1);
                            GlycoSeq.GlycanSequencing _Gs = new GlycoSeq.GlycanSequencing(_msScan, y1mz, y1cs, true, true, false, true, false, output_folder, true, 50, 500);

                            int FullStructureCount = _Gs.FullSequencedStructures.Count;
                            if (FullStructureCount > 0)
                            {
                                bool debug = true;
                            }
                            int PartialStructureCount = _Gs.SequencedStructures.Count ;
                            if (PartialStructureCount > 0)
                            {
                                bool debug = true; 
                            }
                        }

                    }
                }



                // Precursor init
               /* GlypID.HornTransform.clsHornTransformResults[] TransformResults = new GlypID.HornTransform.clsHornTransformResults[1]; 
                TransformResults[0] = m._RepresentativeTransformResult;
                
                //CID scoring
                GlypID.Peaks.clsPeak[] CIDPeaks = m._RepresentativeCIDPeaks;
                
                // HCD Scoring
                GlypID.Peaks.clsPeak[] HCDPeaks = m._RepresentativeHCDPeaks;
                GlypID.HCDScoring.clsHCDScoringScanResults[] HCDScoreResults = new GlypID.HCDScoring.clsHCDScoringScanResults[1];
                double hcdScore = _HCDScoring.ScoreHCDSpectra(ref HCDPeaks, ref mz, ref intensities, ref TransformResults, ref HCDScoreResults); 

                // ETD Scoring
                GlypID.Peaks.clsPeak[] ETDPeaks = m._RepresentativeETDPeaks ; 
                GlypID.ETDScoring.clsETDScoringScanResults ETDScoreResults = new GlypID.ETDScoring.clsETDScoringScanResults() ; 
                double score = _ETDScoring.ScoreETDSpectra(ref ETDPeaks, ETDScoreResults) ;
                
            }
        }

        public void ScoreMapBestMatch(ref Classes.MapRecord _glycoMap)
        {
           /* for (int i = 0; i < _glycoMap._AllMLNRecords.Count; i++)
            {
                Classes.MultiAlignRecord m = _glycoMap._AllMLNRecords[i];
                int num_records = m._AssociatedUMCRecords.Count;
                double min_ppm = _ETDScoring.ScoringParameters.PPMTolerance;
                int min_ppm_index = -1;
                for (int j = 0; j < num_records; j++)
                {
                    double ppm = _utils.CalculateDelMassPPM(m.Mass, m._AssociatedUMCRecords[j].TransformResult.mdbl_mono_mw);
                    if (ppm < min_ppm)
                    {
                        min_ppm_index = j;
                        min_ppm = ppm;
                    }
                }
                Classes.UMCRecord u  = m._AssociatedUMCRecords[min_ppm_index];
            }
        }

        public int DetermineBestMatchBasedOnMassAndGlycanType(GlypID.ETDScoring.clsETDScoringScanResults[] etdCandidates, double gp_mass, GlypID.enmGlycanType gType)
        {
           /* int best_index = -1;
            double min_ppm = _ETDScoring.ScoringParameters.PPMTolerance; 
            for (int k = 0; k < etdCandidates.Length; k++)
            {
                GlypID.Glycan.clsGlycan gc = new GlypID.Glycan.clsGlycan();
                //bool process_result = false;
                bool process_result = true; 
                if (gType == GlypID.enmGlycanType.CS)
                {
                    bool is_sialylated = gc.GlycanCompositionHasNeuAC(etdCandidates[k].mstr_glycan_composition);
                    if (is_sialylated)
                        process_result = true;
                    else
                        process_result = false;
                }
                else
                    process_result = true;
                if (process_result)
                {
                    double glycomass = _utils.CalculateGPMass(etdCandidates[k].mdbl_seq_mass, etdCandidates[k].mdbl_glycan_mass);
                    double ppm_diff = _utils.CalculateDelMassPPM(glycomass, gp_mass);
                    if (ppm_diff < min_ppm)
                    {
                        best_index = k;
                        min_ppm = ppm_diff;
                    }
                }
            }
            return best_index; 
        }
        
       

        private void ScoreHCDPeaks(ref Classes.UMCRecord u)
        {
            // Get out parent
           /* GlypID.HornTransform.clsHornTransformResults[] TransformResults = new GlypID.HornTransform.clsHornTransformResults[1];
            TransformResults[0] = u.TransformResult;
            if (u.HCDPeaks != null)
            {
                GlypID.Peaks.clsPeak[] HCDPeaks = u.HCDPeaks;
                float[] hcdmzs = u.HCDMzs;
                float[] hcdintensities = u.HCDIntensities;
                GlypID.HCDScoring.clsHCDScoringScanResults[] HCDScoreResults = new GlypID.HCDScoring.clsHCDScoringScanResults[1];
                double score = _HCDScoring.ScoreHCDSpectra(ref HCDPeaks, ref hcdmzs, ref hcdintensities, ref TransformResults, ref HCDScoreResults);
                u.HCDScore = HCDScoreResults[0].mdbl_hcd_score;
                u.GlycanType = (GlypID.enmGlycanType) HCDScoreResults[0].menm_glycan_type; 
            }
        }
        
                 
        private void ScoreETDPeaks(ref Classes.UMCRecord u, GlypID.ETDScoring.clsETDScoringScanResults[] thisETDScoreResults )
        {
           /* if (u.ETDPeaks[0] != null)
            {               
                // Score each candidate
                int best_scoring_index = -1;
                double max_score = 0;
                string max_score_peptide = "";
                string max_score_protein = "";
                string max_score_glycan = "";
                for (int k = 0; k < thisETDScoreResults.Length; k++)
                {
                    GlypID.ETDScoring.clsETDScoringScanResults thisResult = thisETDScoreResults[k];
                 
                    bool process_result = false;
                    if (u.GlycanType == GlypID.enmGlycanType.CS)
                    {
                        GlypID.Glycan.clsGlycan gc = new GlypID.Glycan.clsGlycan();
                        bool is_sialylated = gc.GlycanCompositionHasNeuAC(thisResult.mstr_glycan_composition);
                        if (is_sialylated)
                            process_result = true;
                        else
                            process_result = false;
                    }
                    else
                        process_result = true;


                    if (thisResult.mstr_pep_seq != "" && process_result)
                    {
                        GlypID.Peaks.clsPeak[] ETDPeaks = u.ETDPeaks;
                        double etd_score = _ETDScoring.ScoreETDSpectra(ref ETDPeaks, thisResult);
                        if (etd_score > max_score)
                        {
                            thisETDScoreResults[k].mdbl_etd_score = etd_score;                            
                            best_scoring_index = k;
                            max_score = etd_score;
                            max_score_peptide = thisResult.mstr_pep_seq;
                            max_score_protein = thisResult.mstr_pro_seq_name;
                            max_score_glycan = thisResult.mstr_glycan_composition;
                        }
                    }
                }
                u.ETDScore = max_score;
                u.GlycanComposition = max_score_glycan;
                u.PeptideSeq = max_score_peptide;
                u.ProteinName = max_score_protein;
            }
        }

        
        private void ScoreCIDPeaks(ref Classes.UMCRecord u)
        {
             // Get out parent
           /* GlypID.HornTransform.clsHornTransformResults[] TransformResults = new GlypID.HornTransform.clsHornTransformResults[1];
            TransformResults[0] = u.TransformResult;

            // Now score the CID
            if (u.CIDPeaks[0] != null)
            {
                GlypID.Peaks.clsPeak[] CIDPeaks = u.CIDPeaks;
                float[] mzs = new float[CIDPeaks.Length];
                float[] intensities = new float[CIDPeaks.Length];
                int index = 0;
                foreach (GlypID.Peaks.clsPeak peak in CIDPeaks)
                {
                    mzs[index] = Convert.ToSingle(CIDPeaks[index].mdbl_mz);
                    intensities[index] = Convert.ToSingle(CIDPeaks[index].mdbl_intensity);
                }

                GlypID.CIDScoring.clsCIDScoringScanResults[] CIDScoreResults = new GlypID.CIDScoring.clsCIDScoringScanResults[1];
                bool found_score = _CIDScoring.ScoreCIDSpectra(ref CIDPeaks, ref mzs, ref intensities, ref TransformResults, ref CIDScoreResults);
                u.CIDScore = CIDScoreResults[0].mdbl_cid_score;
            }
        }

        


        public void SequenceCIDPeaks(ref Classes.MultiAlignRecord m, List<string> peptides)
        {

           /* GlypID.HornTransform.clsHornTransformResults[] TransformResults = new GlypID.HornTransform.clsHornTransformResults[1];

            

            for (int i = 0; i < m._AssociatedUMCRecords.Count; i++)
            {
                Classes.UMCRecord u = m._AssociatedUMCRecords[i];
                if (u.DatasetID == m._RepresentativeDatasetID_CID)
                {
                    TransformResults[0] = u.TransformResult;
                    GlypID.Peaks.clsPeak[] CIDPeaks = u.CIDPeaks;

                    float[] mzs = new float[CIDPeaks.Length];
                    float[] intensities = new float[CIDPeaks.Length];
                    int index = 0;
                    foreach (GlypID.Peaks.clsPeak peak in CIDPeaks)
                    {
                        mzs[index] = Convert.ToSingle(CIDPeaks[index].mdbl_mz);
                        intensities[index] = Convert.ToSingle(CIDPeaks[index].mdbl_intensity);
                    }

                    GlypID.CIDScoring.clsCIDScoringScanResults[] CIDScoreResults = new GlypID.CIDScoring.clsCIDScoringScanResults[1];
                    bool found_score = _CIDScoring.ScoreCIDSpectra(ref CIDPeaks, ref mzs, ref intensities, ref TransformResults, ref CIDScoreResults);

                    // get yx ions
                    //GlycoSeq.GetYxPeaks(CIDPeaks); 

                    //Calculate peptide masses.
                    



                    // Do sequencing
                    foreach (string s in peptides)
                    {
                        string PeptideSeq = s;
                        GlycoSeq.AminoAcidWeight AAMW = new GlycoSeq.AminoAcidWeight();
                        float PeptideMass = AAMW.GetMonoMW(PeptideSeq, true);
                        int y1cs = 2;
                        float y1mz = (PeptideMass + GlycoSeq.GlycanMass.GetMZ(GlycoSeq.Glycan.Types.HexNAc, 1) + (float)_utils._CC_MASS * y1cs) / y1cs;

                        double y1mtest = _utils.CalculateMz((PeptideMass + 203.007937), 2);

                        GlycoSeq.MSScan _msScan = new GlycoSeq.MSScan(CIDPeaks, mzs, intensities, Convert.ToSingle(TransformResults[0].mdbl_mz), Convert.ToInt32(TransformResults[0].mshort_cs), 1);

                        bool hasNeuAc = true;
                        if (u.GlycanType == GlypID.enmGlycanType.CS)
                            hasNeuAc = true;
                        else
                        {
                            if (!(u.GlycanType == GlypID.enmGlycanType.NA))
                                hasNeuAc = false;
                        }

                        GlycoSeq.GlycanSequencing _Gs = new GlycoSeq.GlycanSequencing(_msScan, y1mz, y1cs);
                        List<GlycoSeq.GlycanStructure> fullstructures = _Gs.FullSequencedStructures;
                        fullstructures[0].Score; 
                        

                       /* int FullStructureCount = _Gs.FullSequencedStructures.Count;
                        if (FullStructureCount > 0)
                        {
                            bool debug = true;
                        }
                        int PartialStructureCount = _Gs.SequencedStructures.Count;
                        if (PartialStructureCount > 0)
                        {
                            bool debug = true;
                        }*/
                   /* }
                }
            }    *    
            

        }



        public void ScoreMapIndividual(ref Classes.MapRecord _glycoMap, ref Classes.Params _params)
        {
           /* float[] mz = { 100, 200, 300 };
            float[] intensities = { 100, 50, 50 };
          

            for (int i = 0; i < _glycoMap._AllMLNRecords.Count; i++)
            {
                Classes.MultiAlignRecord m = _glycoMap._AllMLNRecords[i];

                if (m.ID == 1)
                {
                    bool debug = true;
                }

               
                // Get out all candidates
                GlypID.ETDScoring.clsETDScoringScanResults[] ETDScoreResults = new GlypID.ETDScoring.clsETDScoringScanResults[m._PeptideSeq.Length];

                List <string> peptides_as_string = new List<string>() ; 
                
                for (int k=0 ; k < m._PeptideSeq.Length; k++)
                {
                    ETDScoreResults[k] = new GlypID.ETDScoring.clsETDScoringScanResults() ; 
                    ETDScoreResults[k].mstr_pep_seq = m._PeptideSeq[k] ;
                    ETDScoreResults[k].mstr_pro_seq_name = m._ProteinName[k]; 
                    ETDScoreResults[k].mdbl_seq_mass = m._PeptideMass[k] ; 
                    ETDScoreResults[k].mdbl_glycan_mass = m._GlycanMass[k] ;
                    ETDScoreResults[k].mstr_glycan_composition = m._GlycanComposition[k];
                    ETDScoreResults[k].mstr_glyco_site = m._NGlycoSite[k]; 
                    peptides_as_string.Add(ETDScoreResults[k].mstr_pep_seq) ; 
                }

                int num_records = m._AssociatedUMCRecords.Count;
                for (int j = 0; j < num_records; j++)
                {
                    Classes.UMCRecord u = m._AssociatedUMCRecords[j];                  

                    // First score HCD
                    if (_glycoMap._IsHCD[u.DatasetID])
                    {
                        if (u.HCDPeaks[0] != null)
                        {
                                ScoreHCDPeaks(ref u) ;                               
                        }
                    }
                    
                    // score up all th ETD 
                    if (u.HCDScore < 1)
                    {
                        if (_glycoMap._IsETD[u.DatasetID])
                        {
                            ScoreETDPeaks(ref u, ETDScoreResults);
                        }
                        else
                        {
                            // no ETD peaks in the dataset, so use closest match based on mass error                                              
                            //u.GlycanType = GlypID.enmGlycanType.NA; Taken out July 9, 2012, need to see if this changes PooledTargeted

                            int best_matching_index = DetermineBestMatchBasedOnMassAndGlycanType(ETDScoreResults, m.Mass, u.GlycanType);
                            u.ETDScore = 0;
                            if (best_matching_index >= 0)
                            {
                                u.GlycanComposition = ETDScoreResults[best_matching_index].mstr_glycan_composition;
                                u.PeptideSeq = ETDScoreResults[best_matching_index].mstr_pep_seq;
                                u.ProteinName = ETDScoreResults[best_matching_index].mstr_pro_seq_name;
                                u.PeptideMass = ETDScoreResults[best_matching_index].mdbl_seq_mass;
                                u.GlycanMass = ETDScoreResults[best_matching_index].mdbl_glycan_mass;
                                u.NGlycoSite = ETDScoreResults[best_matching_index].mstr_glyco_site;
                            }
                        }
                    }  
                              
                    if (_glycoMap._IsCID[u.DatasetID])
                    {
                        ScoreCIDPeaks(ref u);
                        if (u.DatasetID == m._RepresentativeDatasetID_CID)
                        {
                            m.RepCIDScore = (int) u.CIDScore;
                        }
                    }


                }                  
            }
        }
    }*/
}
