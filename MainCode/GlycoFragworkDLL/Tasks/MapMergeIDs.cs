//using System;
//using System.Collections.Generic;
//using System.Collections;
//using System.Linq;
//using System.Text;
//using GlycoFragworkDLL.Classes;
//using GlycoFragworkDLL.Utils;

//namespace GlycoFragworkDLL.Tasks
//{
//    public class MapMergeIDs
//    {
//        Classes.MapRecord _GlycoMap;
//        Utilities _utils;
//        List<NonGlycopeptideRecord> _nonglycopeptides;
//        Dictionary<int, List<NonGlycopeptideRecord>> _dNonglycopeptides;



//        GlypID.HornTransform.clsHornTransformParameters _TransformParameters;
//        GlypID.HornTransform.clsHornTransform _HornTransform;
//        GlypID.HornTransform.clsHornTransformResults[] _ParentTransformResults;

//        GlypID.Peaks.clsPeakProcessor _ParentPeakProcessor;
//        GlypID.Peaks.clsPeakProcessorParameters _ParentPeakProcessorParams;
//        GlypID.Peaks.clsPeak[] _ParentPeaks;


//        GlypID.Peaks.clsPeakProcessor _MSMSPeakProcessor;
//        GlypID.Peaks.clsPeakProcessorParameters _MSMSPeakProcessorParams;
//        GlypID.Peaks.clsPeak[] _MSMSPeaks;

//        GlypID.Scoring.clsScoringParameters _ScoringParameters;


//        double _PPM_Diff = 20;
//        double _DA_Diff = 5; 

//        GlypID.Readers.clsRawData _RawData;

//        public MapMergeIDs()
//        {
//            _utils = new Utilities(); 
//            _nonglycopeptides = new List<NonGlycopeptideRecord>();

//            _TransformParameters = new GlypID.HornTransform.clsHornTransformParameters();
//            _ParentPeakProcessorParams = new GlypID.Peaks.clsPeakProcessorParameters();
//            _MSMSPeakProcessorParams = new GlypID.Peaks.clsPeakProcessorParameters();
//            _ScoringParameters = new GlypID.Scoring.clsScoringParameters();

//            _HornTransform = new GlypID.HornTransform.clsHornTransform();
//            _ParentPeakProcessor = new GlypID.Peaks.clsPeakProcessor();
//            _MSMSPeakProcessor = new GlypID.Peaks.clsPeakProcessor();

//            _ParentPeaks = new GlypID.Peaks.clsPeak[1];
//            _MSMSPeaks = new GlypID.Peaks.clsPeak[1];
//            _ParentTransformResults = new GlypID.HornTransform.clsHornTransformResults[1];
//            _dNonglycopeptides = new Dictionary<int, List<NonGlycopeptideRecord>>();
//        }


//        public void LoadParameters(GlycoFragworkDLL.Classes.Params thisParams)
//        {
//            // To Do
//            _TransformParameters = thisParams.TransformParams;
//            _ParentPeakProcessorParams = thisParams.PeakParams;
//            _MSMSPeakProcessorParams = thisParams.PeakParams;
//            InitParameters();
//            _PPM_Diff = thisParams.ScoringParams.PPMTolerance;

//        }

//        public void InitParameters()
//        {
//            // To Do
//            _HornTransform.TransformParameters = _TransformParameters;
//            _ParentPeakProcessor.SetOptions(_ParentPeakProcessorParams);
//            _MSMSPeakProcessor.SetOptions(_MSMSPeakProcessorParams);
//            _ParentPeakProcessor.ProfileType = GlypID.enmProfileType.PROFILE;
//            _MSMSPeakProcessor.ProfileType = GlypID.enmProfileType.CENTROIDED;
//            _ParentPeakProcessor.ProfileType = GlypID.enmProfileType.PROFILE;
//            _MSMSPeakProcessor.ProfileType = GlypID.enmProfileType.CENTROIDED;
//        }



//        /// <summary>
//        /// Function to create a hash table/Dictionary format of all glycopeptides that makes searching faster
//        /// </summary>
//        public void CreateHashTableNonGlycopeptides(double min_mass, double max_mass, bool filter_mass)
//        {

//            _dNonglycopeptides = new Dictionary<int, List<NonGlycopeptideRecord>>();

//            foreach (NonGlycopeptideRecord g in _nonglycopeptides)
//            {
//                List<NonGlycopeptideRecord> value_nongp_list = new List<NonGlycopeptideRecord>();
//                int key_mass = (int)Math.Floor(g.Mass);

//                int min_mass_int = (int)Math.Floor(min_mass);
//                int max_mass_int = (int)Math.Floor(max_mass);

//                bool proceed = false;
//                if (filter_mass)
//                {
//                    if (key_mass >= min_mass_int && key_mass <= max_mass_int)
//                        proceed = true;
//                }
//                else
//                    proceed = true;

//                if (proceed)
//                {
//                    if (!_dNonglycopeptides.ContainsKey(key_mass))
//                    {
//                        if (value_nongp_list.Count > 0)
//                            value_nongp_list.Clear();
//                        value_nongp_list.Add(g);
//                        _dNonglycopeptides.Add(key_mass, value_nongp_list);
//                    }
//                    else
//                    {
//                        value_nongp_list.Clear();
//                        value_nongp_list = _dNonglycopeptides[key_mass];
//                        value_nongp_list.Add(g);
//                        _dNonglycopeptides[key_mass] = value_nongp_list;

//                    }
//                }

//            }
//        }

//        public bool SearchNonGlycopeptides(double mass, ref List<NonGlycopeptideRecord> ngps)
//        {
//            bool found_match = false;
//            ngps = new List<NonGlycopeptideRecord>();

//            double min_ppm = _PPM_Diff;

//            foreach (NonGlycopeptideRecord g in _nonglycopeptides)
//            {
//                if (_utils.CalculateDelMassPPM(g.Mass, mass) < min_ppm)
//                {
//                    found_match = true;
//                    ngps.Add(g); 
//                    g.unassigned = false ;                          
//                }
//                if (g.Mass > mass + 0.5)
//                    break; 
//            }

            

//            return found_match; 

//        }


//         public void LoadIdentificationsFromFile(string infile, List <string> datasetnames)
//         {
//             _nonglycopeptides = new List<NonGlycopeptideRecord>();
//             List<string> datanames_stripped = new List<string>();

//             double min_mass = Double.MaxValue;
//             double max_mass = 0; 

//             for(int i=0 ; i < datasetnames.Count ; i++)
//             {
//                 int folderpos = datasetnames[i].LastIndexOf("\\");
//                 int dotpos = datasetnames[i].LastIndexOf(".");
//                 string temp = datasetnames[i].Substring(folderpos + 1, dotpos - folderpos);
//                 datanames_stripped.Add(temp); 
//             }

//             // Read in a formatted Build summary output
//            CSVFileHandler CsvFileHandler2 = new CSVFileHandler(infile, CSVFileHandler.READ_ONLY);
//            CsvFileHandler2.openFile();
//            String[] Attributes2;
//            CsvFileHandler2.readLine();
//            while ((Attributes2 = CsvFileHandler2.readLine()) != null)
//            {
//                double mass = double.Parse(Attributes2[2]);

//                GlypID.Sequence.clsSequence sequence = new GlypID.Sequence.clsSequence();

//                if (datanames_stripped.Exists(element => element == Attributes2[0]))
//                {
//                    int index = datanames_stripped.FindIndex(element => element == Attributes2[0]);
//                    NonGlycopeptideRecord ngp = new NonGlycopeptideRecord();
                    
//                    sequence.sequence = Attributes2[3];
//                    sequence.proteinName = Attributes2[6];

                    
//                    ngp.Sequence = sequence;
//                    ngp.DatasetScan = int.Parse(Attributes2[1]);
//                    ngp.DatasetID = index;
//                    ngp.Mass = double.Parse(Attributes2[2]);
//                    ngp.CS = short.Parse(Attributes2[5]);
//                    _nonglycopeptides.Add(ngp);

//                    if (mass < min_mass)
//                        min_mass = mass;
//                    if (mass > max_mass)
//                        max_mass = mass; 
//                }
//            }

           
            

//            _nonglycopeptides.Sort(delegate(NonGlycopeptideRecord g1, NonGlycopeptideRecord g2) { return g1.Mass.CompareTo(g2.Mass); });
//            CreateHashTableNonGlycopeptides(min_mass, max_mass, true); 
//         }

//         public void MergeIDsIntoMap(ref Classes.MapRecord _glycoMap, ref Classes.Params _parameters)
//         {
//             Classes.MapRecord _tMap = new MapRecord(); 


//             _glycoMap._AllMLNRecords.ForEach(delegate(MultiAlignRecord m)
//             {
//                 bool process_mrecord = false;
//                 Console.WriteLine("UMC = " + m.ID);

//                 List<NonGlycopeptideRecord> matched_nongps = new List<NonGlycopeptideRecord>(); 
//                 process_mrecord = SearchNonGlycopeptides(m.Mass, ref matched_nongps);

//                 if (process_mrecord)
//                 {
//                     for (int i = 0; i < m._AssociatedUMCRecords.Count; i++)
//                     {
//                         // See if UMC was present in this dataset
//                         //+ "\tDataset = " + i); 
//                         UMCRecord u = m._AssociatedUMCRecords[i];
//                         int datasetid = u.DatasetID;
//                         bool foundId = false;
//                         double min_ppm = _PPM_Diff;

//                         if (u.ScanRep > 0)
//                         {
//                             for (int k = 0; k < matched_nongps.Count; k++)
//                             {
//                                 if ((matched_nongps[k].DatasetID == datasetid) &
//                                     (matched_nongps[k].DatasetScan >= u.ScanStart) &
//                                     (matched_nongps[k].DatasetScan <= u.ScanEnd))
//                                 {
//                                     if (foundId)
//                                     {
//                                         if (u.PeptideSeq != matched_nongps[k].Sequence.sequence)
//                                         {
//                                             // To do
//                                             bool error = false;
//                                             throw new Exception("More work needed");

//                                         }
//                                     }
//                                     else
//                                     {
//                                         foundId = true;
//                                         u.ProteinName = matched_nongps[k].Sequence.proteinName;
//                                         u.PeptideSeq = matched_nongps[k].Sequence.sequence;                                         
//                                     }

//                                 }
//                             }
//                         }
//                     }

//                     if (ResolveIds(ref m))
//                     {
//                        _tMap.AddRecord(m);
//                     }


//                 }


//             }) ;

//             if (_tMap._AllMLNRecords.Count != _glycoMap._AllMLNRecords.Count)
//             {
//                 _glycoMap.ClearRecords();
//                 for (int k = 0; k < _tMap._AllMLNRecords.Count; k++)
//                 {
//                     _glycoMap.AddRecord(_tMap._AllMLNRecords[k]);
//                 }
//             }

//         }


//         private bool ResolveIds(ref Classes.MultiAlignRecord m)
//         {
//            bool is_consistent = false;
//            string set_protein = "";
//            string set_peptide = "";            
//            Classes.UMCRecord finalU = new GlycoFragworkDLL.Classes.UMCRecord() ; 

//            for (int i = 0; i < m._AssociatedUMCRecords.Count; i++)
//            {
//                Classes.UMCRecord u = m._AssociatedUMCRecords[i] ;
//                if (u.PeptideSeq.Length > 0)
//                {
//                    if (set_protein == "")
//                    {
//                        set_protein = u.ProteinName;
//                        set_peptide = u.PeptideSeq;                        
//                        is_consistent = true; 
//                        finalU = u ; 
//                    }
//                    else
//                    {
//                        if ((set_protein != u.ProteinName) || (set_peptide != u.PeptideSeq))
//                        {
//                            is_consistent = false;                            
//                            break;
//                        }
//                    }
//                }               
//            }
//            if (is_consistent && set_protein!="")
//            {
//                m.BestMatchPeptideSeq = finalU.PeptideSeq ; 
//                m.BestMatchPeptideMass = finalU.PeptideMass ;                 
//                m.BestMatchProteinName = finalU.ProteinName; 
//            }
//            return is_consistent;                 
//        }
         



//    }
//}