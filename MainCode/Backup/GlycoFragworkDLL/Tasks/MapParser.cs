using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
//using System.Drawing;
using GlycoFragworkDLL.Classes;
using GlycoFragworkDLL.Utils;

namespace GlycoFragworkDLL.Tasks
{
    public class MapParser
    {
        Utilities _utils;
        Utils.SpectralUtilities _CIDSpectralUtilities;
        Utils.SpectralUtilities _HCDSpectralUtilities;
        Utils.SpectralUtilities _ETDSpectralUtilities;
        GlycoFragworkDLL.Classes.Params _params;


        GlypID.Scoring.clsScoringParameters _ScoringParameters;

        public MapParser()
        {
            _utils = new GlycoFragworkDLL.Utils.Utilities();
            _ScoringParameters = new GlypID.Scoring.clsScoringParameters();
            _CIDSpectralUtilities = new SpectralUtilities();
            _HCDSpectralUtilities = new SpectralUtilities();
            _ETDSpectralUtilities = new SpectralUtilities();
            _params = new Params();
        }

        /// <summary>
        ///  Function to calculate a representative HCD spectrum.
        ///  For every mln record, the set od hcd records with the minimum score is chosen.
        ///  The spectrum with the greates SNR is reported.
        /// </summary>
        /// <param name="_glycoMap"></param>
        /// <param name="_params"></param>
        public void CalculateHCDRepresentativeFragmentationSpectra(ref Classes.MapRecord _glycoMap, ref Classes.Params _params)
        {
           /* int Id = 0;

            for (int i = 0; i < _glycoMap._AllMLNRecords.Count; i++)
            {
                MultiAlignRecord m = new MultiAlignRecord();
                m = _glycoMap._AllMLNRecords[i];
                int num_records = m._AssociatedUMCRecords.Count;
                double min_score = m.MinHCDScore();
                if (min_score < 1)
                {
                    GlycoFragworkDLL.Utils.SpectralUtilities _HCDSpectralUtilities = new GlycoFragworkDLL.Utils.SpectralUtilities();

                    List<string> cluster_names = new List<string>();
                    string cum_spectra_names = null;
                    for (int j = 0; j < m._AssociatedUMCRecords.Count; j++)
                    {
                        Classes.UMCRecord u = m._AssociatedUMCRecords[j];
                        int spectra_num = 0;
                        if (u.HCDScore == min_score)
                        {
                            string spectra_name = "Spectra_" + u.DatasetID + "_" + spectra_num + "_" + u.HCDScanRep;
                            spectra_num++;
                            GlypID.Peaks.clsPeak[] thisHCDPeaks = new GlypID.Peaks.clsPeak[u.HCDPeaks.Length];
                            Array.Copy(u.HCDPeaks, thisHCDPeaks, u.HCDPeaks.Length);
                            _HCDSpectralUtilities.AddPeaksToList(ref thisHCDPeaks, spectra_name);
                            if (cum_spectra_names != null)
                                cum_spectra_names = cum_spectra_names + "-" + spectra_name;
                            else
                                cum_spectra_names = spectra_name;
                        }
                    }


                    if (cum_spectra_names.Length > 0)
                    {
                        cluster_names.Add(cum_spectra_names);
                        _HCDSpectralUtilities.AssignClusters(ref cluster_names);
                        GlypID.Peaks.clsPeak[] repHCDPeaks = new GlypID.Peaks.clsPeak[0];
                        int repOrigIndex = _HCDSpectralUtilities.GetRepresentativePeaksFromCluster(0, ref repHCDPeaks, _params.ScoringParams.MinHCDMz, _params.ScoringParams.MaxHCDMz, true);
                        if (repHCDPeaks.Length > 1)
                        {
                            m._RepresentativeHCDPeaks = repHCDPeaks;
                            m._RepresentativeDatasetID_HCD = repOrigIndex;
                            m.RepHCDScore = min_score;
                        }
                    }
                }
            }*/
        }

       /* public void CalculateCIDRepresentativeFragmentationSpectra(ref Classes.UMCRecord _u, ref Classes.Params _params)
        {
            Classes.UMCRecord _tempU = new UMCRecord();
            int spectra_num = 0;
            string cum_spectra_names = null;
            List<string> cluster_names = new List<string>();
            List<int> orphan_ids = new List<int>();
            for (int i = 0; i < _u._AssociatedFragEvents.Count; i++)
            {
                if (_u._AssociatedFragEvents[i].CIDPeaks[0] != null)
                {
                    string spectra_name = "Spectra_" + _u.DatasetID + "_" + spectra_num + "_" + _u._AssociatedFragEvents[i].CIDScan;
                    GlypID.Peaks.clsPeak[] thisCIDPeaks = new GlypID.Peaks.clsPeak[_u._AssociatedFragEvents[i].CIDPeaks.Length];
                    Array.Copy(_u._AssociatedFragEvents[i].CIDPeaks, thisCIDPeaks, _u._AssociatedFragEvents[i].CIDPeaks.Length);
                    _CIDSpectralUtilities.AddPeaksToList(ref thisCIDPeaks, spectra_name);
                    spectra_num++;
                    if (cum_spectra_names != null)
                        cum_spectra_names = cum_spectra_names + "-" + spectra_name;
                    else
                        cum_spectra_names = spectra_name;
                }
                else
                {
                    // TO check
                    // This happens which means there was no fragmentation event then (or) HCD score was bad
                    if (_u._AssociatedFragEvents[i].HCDScore < 1)
                        orphan_ids.Add(_u.DatasetID);
                }



            }

        }*/

        public void AssignFDR(ref Classes.MapRecord _glycoMap, ref Classes.Params _params)
        {
            bool use_combined_score = false; // To do change this based on _params
            Classes.FragEvents e = new FragEvents();
            
            for (int i = 0 ; i < _glycoMap._AllMLNRecords.Count ; i++)
            {                
                MultiAlignRecord m = new MultiAlignRecord();
                m = _glycoMap._AllMLNRecords[i];
             
                for (int j = 0; j < m._ClusterRepFragEvents.Count; j++)
                {
                    if (m._ClusterRepFragEvents[j].ETDScore > 0)
                    {
                        if (m._ClusterRepFragEvents[j].GP_Record.IsDecoy)
                            _glycoMap._AllFalseHitsFDRScore.Add(m._ClusterRepFragEvents[j].ETDScore);
                        else
                            _glycoMap._AllTrueHitsFDRScore.Add(m._ClusterRepFragEvents[j].ETDScore); 

                    }
                }
            }
            
            
            // Assign ETD score based on ETD type
            for (int i = 0; i < _glycoMap._AllMLNRecords.Count; i++)
            {
                MultiAlignRecord m = new MultiAlignRecord();
                m = _glycoMap._AllMLNRecords[i];

               
                for (int j = 0; j < m._ClusterRepFragEvents.Count; j++)
                {
                    e = new FragEvents();
                    e = m._ClusterRepFragEvents[j];
                    if (e.ETDScore > 0)
                    {
                        int num_false_hits = 0;
                        int num_true_hits = 0;
                        foreach (double score in _glycoMap._AllFalseHitsFDRScore)
                        {
                            if (score > e.ETDScore)
                                num_false_hits++;
                        }
                        foreach (double score in _glycoMap._AllTrueHitsFDRScore)
                        {
                            if (score > e.ETDScore)
                                num_true_hits++;
                        }
                        if (num_true_hits == 0)
                            num_true_hits = 1; // to avoid divide by 0 for the highest.
                        e.FDR = (float)num_false_hits / num_true_hits;
                    }

                }
            }
                   /* int num_records = m._AssociatedUMCRecords.Count;
                    for (int j = 0; j < num_records; j++)
                    {
                        for (int k = 0; k < m._AssociatedUMCRecords[j]._AssociatedFragEvents.Count; k++)
                        {
                            e = new FragEvents();
                            e = m._AssociatedUMCRecords[j]._AssociatedFragEvents[k];

                            if (e.ETDScore > 0)
                            {
                                int num_false_hits = 0;
                                int num_true_hits = 0;
                                foreach (double score in _glycoMap._AllFalseHitsFDRScore)
                                {
                                    if (score > e.ETDScore)
                                        num_false_hits++;
                                }
                                foreach (double score in _glycoMap._AllTrueHitsFDRScore)
                                {
                                    if (score > e.ETDScore)
                                        num_true_hits++;
                                }
                                e.FDR = (float)num_false_hits / num_true_hits;
                            }
                        }
                    }*/
             
        }


        public void GetRepresentatives(ref Classes.MapRecord _glycoMap, ref Classes.Params _params)
        {
            Classes.MapRecord _tempMap = new MapRecord();
            _tempMap._AssociatedDatasetNames = _glycoMap._AssociatedDatasetNames;
            _tempMap._IsCID = _glycoMap._IsCID;
            _tempMap._IsETD = _glycoMap._IsETD;
            _tempMap._IsHCD = _glycoMap._IsHCD;
            _tempMap._AllFalseHitsFDRScore = _glycoMap._AllFalseHitsFDRScore;
            _tempMap._AllTrueHitsFDRScore = _glycoMap._AllTrueHitsFDRScore;

            bool use_etd = false ; 
            if (_glycoMap._IsETD.Contains(true))
                use_etd = true ; 
            

            FragEvents tempE = new FragEvents(); 
            for (int i = 0; i < _glycoMap._AllMLNRecords.Count; i++)
            {
                MultiAlignRecord m = new MultiAlignRecord();
                m = _glycoMap._AllMLNRecords[i];

                int num_clusters = m._ClusterNames.Count;

                if (m.ID == 2290)
                {

                    bool test = true;
                }
                for (int c = 0; c < num_clusters; c++)
                {
                    string clustername = m._ClusterNames[c];
                    string[] spectra = clustername.Split('-');

                    Dictionary <int, List<int>> umcs_frag_ids  = new Dictionary<int, List<int>>() ; 
                    List<int>frag_to_look  = new List<int> () ; 


                    foreach (string s in spectra)
                    {
                        string[] parts = s.Split('_');
                        int umc_id = Convert.ToInt32(parts[1]);
                        int frag_id = Convert.ToInt32(parts[2]) ; 
                        if (umcs_frag_ids.ContainsKey(umc_id))                        
                        {
                            umcs_frag_ids[umc_id].Add(frag_id); 
                        }
                        else
                        {
                            frag_to_look.Clear() ; 
                            frag_to_look.Add(frag_id) ; 
                            umcs_frag_ids.Add(umc_id, frag_to_look); 
                        }                        
                    }


                    double max_etd_score = 0;
                    double min_hcd_score = 0 ; 
                    FragEvents maxFragEvent = new FragEvents(); 
                    for (int j = 0; j < m._AssociatedUMCRecords.Count; j++)
                    {
                         if (umcs_frag_ids.ContainsKey(m._AssociatedUMCRecords[j].DatasetID))
                         {
                             tempE = new FragEvents();
                             if (use_etd)
                             {
                                 tempE = m._AssociatedUMCRecords[j].FragEventWithHighestETDScore(true, umcs_frag_ids[m._AssociatedUMCRecords[j].DatasetID]);
                                 if (tempE.ETDScore > max_etd_score)
                                 {
                                     maxFragEvent = tempE;
                                 }
                             }
                             else
                             {
                                 if (!_params.ProcessOnlyNGlycopeptides)
                                     maxFragEvent = m._AssociatedUMCRecords[j].FragEventWithLowestHCDScore(umcs_frag_ids[m._AssociatedUMCRecords[j].DatasetID]);
                             }
                             

                         }
                    }
                    if ((maxFragEvent.ETDScore > 0) || (!_params.ProcessOnlyNGlycopeptides))
                        m._ClusterRepFragEvents.Add(maxFragEvent);
                    /*else
                    {
                        FragEvents minHCDFragEvent = new FragEvents();
                        double min_hcd_score = 1; 
                        for (int j= 0 ; j < m._AssociatedUMCRecords.Count ; j++)
                        {
                            if (umcs_frag_ids.ContainsKey(m._AssociatedUMCRecords[j].DatasetID))
                            {
                                tempE = new FragEvents();
                                tempE = m._AssociatedUMCRecords[j].FragEventWithLowestHCDScore(umcs_frag_ids[m._AssociatedUMCRecords[j].DatasetID]);
                                if (tempE.HCDScore < min_hcd_score)
                                    minHCDFragEvent = tempE;
                            }
                        }
                        m._ClusterRepFragEvents.Add(minHCDFragEvent); 
                    }      */              
                }
            }

           

        }


        public void ResolveIdentificationsAcrossClusters(ref Classes.MapRecord _glycoMap, ref Classes.Params _params)
        {
            Classes.MapRecord _tempMap = new MapRecord();
            _tempMap._AssociatedDatasetNames = _glycoMap._AssociatedDatasetNames;
            _tempMap._IsCID = _glycoMap._IsCID;
            _tempMap._IsETD = _glycoMap._IsETD;
            _tempMap._IsHCD = _glycoMap._IsHCD;
            _tempMap._AllFalseHitsFDRScore = _glycoMap._AllFalseHitsFDRScore;
            _tempMap._AllTrueHitsFDRScore = _glycoMap._AllTrueHitsFDRScore;
            int num_clashes = 0; 

            for (int i = 0; i < _glycoMap._AllMLNRecords.Count; i++)
            {
                MultiAlignRecord m = new MultiAlignRecord();
                m = _glycoMap._AllMLNRecords[i];


                if (_glycoMap._IsETD.Contains(true))
                {
                    if (m._ClusterRepFragEvents.Count > 0)
                    {
                        if (m.HasUniformIDsAcrossClusters())
                        {
                            m.SetIdsBasedOnHighestETD();
                            if (m.IDLabel == _ID_label.Unverified)
                                m.SetIdsBasedOnMassAndGlycanType(_params.ScoringParams.PPMTolerance, true);

                        }
                        else
                        {
                            // Need to think about this
                            bool debug = true;
                            num_clashes++; 
                        }
                    }
                    else
                    {
                        m.SetIdsBasedOnMassAndGlycanType(_params.ScoringParams.PPMTolerance, false);
                    }
                }
                else
                {
                    //m.
                }

                 
            }


            
        }

        public void SequenceCIDPeaks(ref Classes.MapRecord _glycoMap, ref Classes.Params _params, string sequencingFolder)
        {
            Classes.MapRecord _tempMap = new MapRecord();
            _tempMap._AssociatedDatasetNames = _glycoMap._AssociatedDatasetNames;
            _tempMap._IsCID = _glycoMap._IsCID;
            _tempMap._IsETD = _glycoMap._IsETD;
            _tempMap._IsHCD = _glycoMap._IsHCD;

            Classes.FragEvents e = new FragEvents();
            if (sequencingFolder == "")
             sequencingFolder = @"c:\sequencing\"; 

            for (int i = 0; i < _glycoMap._AllMLNRecords.Count; i++)
            {
                MultiAlignRecord m = new MultiAlignRecord();
                m = _glycoMap._AllMLNRecords[i];
                int num_records = m._AssociatedUMCRecords.Count;

                 if (m.ID == 2290)
                {
                    
                    bool test = true;
                }

                //for (int j = 0; j < num_records; j++)
                //{
                   // for (int k = 0; k < m._AssociatedUMCRecords[j]._AssociatedFragEvents.Count; k++)
                for (int k = 0; k < m._ClusterRepFragEvents.Count; k++)
                {
                    e = new FragEvents();
                    e = m._ClusterRepFragEvents[k]; // m._AssociatedUMCRecords[j]._AssociatedFragEvents[k];

                   
                    if (e.ETDScore > 0)
                    {
                        float PeptideMass = 0;
                        float GlcNAcMass = 0;
                        COL.MassLib.MSScan _msScan = new COL.MassLib.MSScan(e.CIDMzs, e.CIDIntensities, Convert.ToSingle(e.TransformResult.mdbl_mz),
                            Convert.ToSingle(e.TransformResult.mdbl_mono_mw), Convert.ToSingle(e.TransformResult.mdbl_average_mw), Convert.ToInt32(e.TransformResult.mshort_cs));
                        if (e.TransformResult.mdbl_average_mw > 0)
                        {
                            PeptideMass = Convert.ToSingle(e.GP_Record.SequenceAverageMass); // this means proper deisotoping has occured so use average compostion
                            GlcNAcMass = COL.GlycoLib.GlycanMass.GetGlycanAVGMass(COL.GlycoLib.Glycan.Type.HexNAc);
                        }
                        else
                        {
                            PeptideMass = Convert.ToSingle(e.GP_Record.SequenceMonoMass); // CS has been assigned, in which case both y1 and precursor should be just mono
                            GlcNAcMass = COL.GlycoLib.GlycanMass.GetGlycanMass(COL.GlycoLib.Glycan.Type.HexNAc);
                        }

                        short y1cs = e.TransformResult.mshort_cs;
                        y1cs--; 
                        while (y1cs > 0)
                        {
                            float y1Mz = ((PeptideMass + GlcNAcMass) + (float)_utils._CC_MASS * y1cs) / y1cs;
                            COL.GlycoSequence.GlycanSequencing _Gs = new COL.GlycoSequence.GlycanSequencing(_msScan, y1Mz, y1cs, e.GP_Record.Glycan.numHex,
                                e.GP_Record.Glycan.numHexNAc, e.GP_Record.Glycan.numDeHex, e.GP_Record.Glycan.numNeuAc, 0, sequencingFolder, true, 0.8f, 60);

                           
                            _Gs.NumbersOfPeaksForSequencing = 140;
                            _Gs.CreatePrecursotMZ = true;
                            _Gs.RewardForCompleteStructure = 3; 

                            if (e.TransformResult.mdbl_average_mw > 0)
                                _Gs.UseAVGMass = true;
                            else
                                _Gs.UseAVGMass = false;
                            int structure_count = _Gs.StartSequencing();
                            if (structure_count > 0)
                            {
                                List<COL.GlycoLib.GlycanStructure> topstructures = _Gs.GetTopRankScoreStructre(1); 
                                e.CIDSequencingScore = topstructures[0].Score;
                                e.GP_Record.GlycanSequence = topstructures[0].IUPACString;
                                
                                // Printing out the sequences
                                string opfile = sequencingFolder + m.ID.ToString() + "_" + e.CIDScan.ToString() + ".txt";
                                Utils.CSVFileHandler cidString = new CSVFileHandler(opfile, CSVFileHandler.WRITE_ONLY);
                                cidString.openFile(); 
                                for (int zz = 0; zz < topstructures.Count; zz++)
                                {
                                    cidString.writeLine(topstructures[zz].IUPACString);
                                    /*COL.GlycoLib.GlycansDrawer _gdraw = new COL.GlycoLib.GlycansDrawer(topstructures[zz].IUPACString, false);
                                    System.Drawing.Image img1 = _gdraw.GetImage();*/
                                }
                                cidString.closeFile(); 
                                break;
                            }
                            y1cs--;
                        }
                    }

                }
                

            }


        }


        /// <summary>
        /// Function to cluster CID spectra in each record.  
        /// </summary>
        /// <param name="_glycoMap"></param>
        /// <param name="_params"></param>
        public void ClusterRecordsOnCID(ref Classes.MapRecord _glycoMap, ref Classes.Params _params)
        {
            Classes.MapRecord _tempMap = new MapRecord();
            _tempMap._AssociatedDatasetNames = _glycoMap._AssociatedDatasetNames;
            _tempMap._IsCID = _glycoMap._IsCID;
            _tempMap._IsETD = _glycoMap._IsETD;
            _tempMap._IsHCD = _glycoMap._IsHCD;

            int Id = 0;

            FragEvents tempE = new FragEvents(); 

            for (int i = 0; i < _glycoMap._AllMLNRecords.Count; i++)
            {
                MultiAlignRecord m = new MultiAlignRecord();
                m = _glycoMap._AllMLNRecords[i];
                int num_records = m._AssociatedUMCRecords.Count;


                short mincs = m.MinChargeStateObserved();
                short maxcs = m.MaxChargeStateObserved();

                for (short thiscs = mincs; thiscs <= maxcs; thiscs++)
                {
                    int spectra_num = 0;
                    string cum_spectra_names = null;
                    List<string> cluster_names = new List<string>();
                    List<int> orphan_ids = new List<int>();
                    for (int j = 0; j < num_records; j++)
                    {
                        if (m._AssociatedUMCRecords[j]._AssociatedFragEvents.Count == 0)
                        {
                            // Not fragmented in this UMC but keep track of this  but

                            // // This happens which means there was no fragmentation event then (or) HCD score was bad
                            if (m._AssociatedUMCRecords[j].Abundance > 0)
                                orphan_ids.Add(m._AssociatedUMCRecords[j].DatasetID);
                        }
                        for (int k = 0; k < m._AssociatedUMCRecords[j]._AssociatedFragEvents.Count; k++)
                        {
                            if ((m._AssociatedUMCRecords[j]._AssociatedFragEvents[k].CIDPeaks[0] != null) &&
                                (m._AssociatedUMCRecords[j]._AssociatedFragEvents[k].TransformResult.mshort_cs == thiscs))
                            {
                                string spectra_name = "Spectra_" + m._AssociatedUMCRecords[j].DatasetID + "_" + m._AssociatedUMCRecords[j]._AssociatedFragEvents[k].ID + "_" +
                                    spectra_num + "_" + m._AssociatedUMCRecords[j]._AssociatedFragEvents[k].CIDScan;

                                GlypID.Peaks.clsPeak[] thisCIDPeaks = new GlypID.Peaks.clsPeak[m._AssociatedUMCRecords[j]._AssociatedFragEvents[k].CIDPeaks.Length];
                                Array.Copy(m._AssociatedUMCRecords[j]._AssociatedFragEvents[k].CIDPeaks, thisCIDPeaks, m._AssociatedUMCRecords[j]._AssociatedFragEvents[k].CIDPeaks.Length);
                                _CIDSpectralUtilities.AddPeaksToList(ref thisCIDPeaks, spectra_name);
                                spectra_num++;
                                if (cum_spectra_names != null)
                                    cum_spectra_names = cum_spectra_names + "-" + spectra_name;
                                else
                                    cum_spectra_names = spectra_name;
                            }
                        }
                    }

                    int num_clusters = _CIDSpectralUtilities.ClusterSpectraInList();
                    _CIDSpectralUtilities.GetClusterNames(ref cluster_names);

                    if (num_clusters > 1)
                    {
                        // Indicates glycoforms
                        bool debug = true;
                        debug = true;
                    }

                    for (int c = 0; c < num_clusters; c++)
                    {  
                        m._ClusterNames.Add(cluster_names[c]);
                    }
                       

                    _CIDSpectralUtilities.Clear();
                }
                _tempMap.AddRecord(m);
            }
                 // Restore them. 
            _glycoMap.ClearRecords();
            for (int i = 0; i < _tempMap._AllMLNRecords.Count; i++)
            {
                _glycoMap.AddRecord(_tempMap._AllMLNRecords[i]);
            }
                    


                 

                    /*else
                    {
                        /* MultiAlignRecord _tempM = new MultiAlignRecord(m);
                         _tempM.ID = Id;
                         _tempMap.AddRecord(_tempM);
                         Id++; */

                    /* if (cum_spectra_names != null)
                     {
                         // Choose the spectra with the greatest SNR
                         num_clusters = 1;
                         cluster_names.Add(cum_spectra_names);
                         _CIDSpectralUtilities.AssignClusters(ref cluster_names);
                     }
                 }*/



                   /* for (int k = 0; k < num_clusters; k++)
                    {
                        // ----- Get a representative spectrum index for each cluster ---//
                        MultiAlignRecord _tempM = new MultiAlignRecord(m);
                        _tempM.ID = Id;
                        _tempM._AssociatedUMCRecords.Clear();

                        GlypID.Peaks.clsPeak[] repCIDPeaks = new GlypID.Peaks.clsPeak[0];
                        int repOrigIndex = _CIDSpectralUtilities.GetRepresentativePeaksFromCluster(k, ref repCIDPeaks, _params.ScoringParams.MinCIDMz, _params.ScoringParams.MaxCIDMz, true);
                        if (repCIDPeaks.Length > 1)
                        {
                            _tempM._RepresentativeCIDPeaks = repCIDPeaks;
                            _tempM._RepresentativeDatasetID_CID = repOrigIndex;
                        }
                        // Attach UMCs corresponding to that cluster.
                        List<int> allOrigIDs = new List<int>();
                        _CIDSpectralUtilities.GetOriginalIDFromCluster(k, ref allOrigIDs);
                        for (int j = 0; j < num_records; j++)
                        {
                            int id = m._AssociatedUMCRecords[j].DatasetID;
                            if (allOrigIDs.Exists(element => element == id))
                            {
                                UMCRecord _tempUMC = new UMCRecord();
                                _tempUMC = m._AssociatedUMCRecords[j];
                                _tempM._AssociatedUMCRecords.Add(_tempUMC);
                            }
                            else if (orphan_ids.Exists(element => element == id)) // This takes care of non fragmentation but still has stuff present
                            {
                                UMCRecord _tempUMC = new UMCRecord();
                                _tempUMC = m._AssociatedUMCRecords[j];
                                _tempM._AssociatedUMCRecords.Add(_tempUMC);
                            }

                        }

                        _tempMap.AddRecord(_tempM);
                        Id++;
                    }*/
                


           

        }


        /// <summary>
        /// Function to resolve all identifications in a mln record.  
        /// If all identifications are consistent, the BestRepresentative information is updated
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        private bool ResolveIdentifications(ref Classes.MultiAlignRecord m)
        {
            bool is_consistent = false;
           /* string set_protein = "";
            string set_site = "";
            string set_glycan = "";
            Classes.UMCRecord finalU = new GlycoFragworkDLL.Classes.UMCRecord();

            for (int i = 0; i < m._AssociatedUMCRecords.Count; i++)
            {
                Classes.UMCRecord u = m._AssociatedUMCRecords[i];
                if (u.PeptideMass > 0)
                {
                    if (set_protein == "")
                    {
                        set_protein = u.ProteinName;
                        set_site = u.NGlycoSite;
                        set_glycan = u.GlycanComposition;
                        is_consistent = true;
                        finalU = u;
                    }
                    else
                    {
                        if ((set_protein != u.ProteinName) || (set_site != u.NGlycoSite) || (set_glycan != u.GlycanComposition))
                        {
                            is_consistent = false;
                            break;
                        }
                    }
                }
            }
            if (is_consistent && set_protein != "")
            {
                m.BestMatchGlycanComposition = finalU.GlycanComposition;
                m.BestMatchGlycanMass = finalU.GlycanMass;
                m.BestMatchPeptideSeq = finalU.PeptideSeq;
                m.BestMatchPeptideMass = finalU.PeptideMass;
                m.BestMatchNGlycoSite = finalU.NGlycoSite;
                m.BestMatchProteinName = finalU.ProteinName;
            }*/
            return is_consistent;
        }


        public void GetRepresentativeIdentificationInformation(ref Classes.MapRecord _glycoMap, ref Classes.Params _params)
        {
            for (int i = 0; i < _glycoMap._AllMLNRecords.Count; i++)
            {
                Classes.MultiAlignRecord m = _glycoMap._AllMLNRecords[i];
                


            }


            // TO DO More
            // IF ETD then resolve things a way that is still to be determined

        }

        public void GetRepresentativesAcrossClusters(ref Classes.MapRecord _glycoMap, ref Classes.Params _params)
        {
            for (int i = 0; i < _glycoMap._AllMLNRecords.Count; i++)
            {
                Classes.MultiAlignRecord m = _glycoMap._AllMLNRecords[i];
                m.SetScoresBasedOnLowestHCD(); 
            }
        }




    }
}
