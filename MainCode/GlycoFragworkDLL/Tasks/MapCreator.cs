using System;
using System.Collections.Generic;
using System.Collections; 
using System.Linq;
using System.Text;
using GlycoFragworkDLL.Classes;
using GlycoFragworkDLL.Utils; 


namespace GlycoFragworkDLL.Tasks
{
    public class MapCreator
    {
        Classes.MapRecord _GlycoMap;
        Utilities _utils; 

        //GlypID.Readers.clsRawData _RawData; AM: changing raw file reader dependency
        Classes.RawFileReader _RawData;

        List <GlycopeptideRecord> _glycopeptides;
       
        Dictionary<int, List<GlycopeptideRecord>> _dGlycopeptides ;
       

        GlypID.HornTransform.clsHornTransformParameters _TransformParameters;
        GlypID.HornTransform.clsHornTransform _HornTransform; 
        GlypID.HornTransform.clsHornTransformResults [] _ParentTransformResults ; 

        GlypID.Peaks.clsPeakProcessor _ParentPeakProcessor ; 
        GlypID.Peaks.clsPeakProcessorParameters _ParentPeakProcessorParams ; 
        GlypID.Peaks.clsPeak [] _ParentPeaks ;
        

        GlypID.Peaks.clsPeakProcessor _MSMSPeakProcessor ; 
        GlypID.Peaks.clsPeakProcessorParameters _MSMSPeakProcessorParams ;
        GlypID.Peaks.clsPeak[] _MSMSPeaks;

        GlypID.Scoring.clsScoringParameters _ScoringParameters;

        GlypID.CIDScoring.clsCIDScoring _CIDScoring;
        GlypID.HCDScoring.clsHCDScoring _HCDScoring;
        GlypID.ETDScoring.clsETDScoring _ETDScoring; 

        

        double _PPM_Diff = 20;
        double _DA_Diff = 6; 

      

        public MapCreator()
        {
            _GlycoMap = new Classes.MapRecord() ;
            _utils = new Utilities(); 
             _glycopeptides = new List<GlycopeptideRecord> ();
             

            _TransformParameters = new GlypID.HornTransform.clsHornTransformParameters() ; 
            _ParentPeakProcessorParams = new GlypID.Peaks.clsPeakProcessorParameters() ; 
            _MSMSPeakProcessorParams = new GlypID.Peaks.clsPeakProcessorParameters() ;
            _ScoringParameters = new GlypID.Scoring.clsScoringParameters(); 

            _HornTransform = new GlypID.HornTransform.clsHornTransform() ; 
            _ParentPeakProcessor = new GlypID.Peaks.clsPeakProcessor() ; 
            _MSMSPeakProcessor = new GlypID.Peaks.clsPeakProcessor() ; 

            _ParentPeaks = new GlypID.Peaks.clsPeak[1] ;
            _MSMSPeaks = new GlypID.Peaks.clsPeak[1];
            _ParentTransformResults = new GlypID.HornTransform.clsHornTransformResults[1];
            _dGlycopeptides = new Dictionary<int, List<GlycopeptideRecord>>();

            _CIDScoring = new GlypID.CIDScoring.clsCIDScoring();
            _HCDScoring = new GlypID.HCDScoring.clsHCDScoring();
            _ETDScoring = new GlypID.ETDScoring.clsETDScoring(); 
            

        }

        public void LoadParameters(GlycoFragworkDLL.Classes.Params thisParams)
        {            
            _TransformParameters = thisParams.TransformParams; 
            _ParentPeakProcessorParams = thisParams.PeakParams ;
            _MSMSPeakProcessorParams = thisParams.PeakParams;

            _HCDScoring.ScoringParameters = _CIDScoring.ScoringParameters = _ETDScoring.ScoringParameters = thisParams.ScoringParams;            

            InitParameters();
            _PPM_Diff = thisParams.ScoringParams.PPMTolerance; 
            
            
            
        }

        public void CleanUp()
        {
            _glycopeptides.Clear();
            _dGlycopeptides.Clear();
            Reset(); 
        }

        public void InitParameters()
        {
            // To Do
            _HornTransform.TransformParameters = _TransformParameters; 
            _ParentPeakProcessor.SetOptions(_ParentPeakProcessorParams);
            _MSMSPeakProcessor.SetOptions(_MSMSPeakProcessorParams);
            _ParentPeakProcessor.ProfileType = GlypID.enmProfileType.PROFILE;
            _MSMSPeakProcessor.ProfileType = GlypID.enmProfileType.CENTROIDED;
            _ParentPeakProcessor.ProfileType = GlypID.enmProfileType.PROFILE;
            _MSMSPeakProcessor.ProfileType = GlypID.enmProfileType.CENTROIDED; 
        }

        

         public void LoadMap(ref Classes.MapRecord _glycoMap, string MLNMapFilename, List<string> DatasetNames,
            List<bool> isCID, List <bool> isHCD, List<bool> isETD, List <string> dataType)
        {
            _glycoMap = new GlycoFragworkDLL.Classes.MapRecord(); 

            int _numdatasets = DatasetNames.Count;

            // read in MLN file
            //_glycoMap.ReadInMLNFile(MLNMapFilename, _numdatasets);            
            _glycoMap.ReanInMLNV2File(MLNMapFilename, _numdatasets); 

            // Assign datasetnames and dataset maps
            _glycoMap._AssociatedDatasetNames = DatasetNames ; 
            _glycoMap._IsCID = isCID; 
            _glycoMap._IsETD = isETD ;
            _glycoMap._IsHCD = isHCD;
            _glycoMap._AssociatedDatasetTypes = dataType; 
        }


        

         /// <summary>
         /// Function to read in glycopeptides dictionary from file
         /// </summary>
         /// <param name="infile"></param>
         public void LoadGlycopeptidesFromFile(string infile, double min_mass, double max_mass, bool create_hash)
         {
                       
            CSVFileHandler CsvFileHandler2 = new CSVFileHandler(infile, CSVFileHandler.READ_ONLY);
            CsvFileHandler2.openFile();
            String[] Attributes2;
            CsvFileHandler2.readLine();

            if (create_hash)
                _dGlycopeptides = new Dictionary<int, List<GlycopeptideRecord>>();
            else
                _glycopeptides = new List<GlycopeptideRecord>();   

                 

            while ((Attributes2 = CsvFileHandler2.readLine()) != null)
            {                
                double gpmass = double.Parse(Attributes2[0]); 
                
                if (gpmass >= min_mass && gpmass <= max_mass)
                {
                    GlycopeptideRecord gp = new GlycopeptideRecord();
                    gp.GP_Mono_Mass = gpmass ; 
                    gp.GP_Average_Mass = double.Parse(Attributes2[1]) ;
                    gp.Sequence.proteinName = Attributes2[2] ; 
                    gp.Sequence.sequence = Attributes2[3]; 
                    gp.Sequence.is_decoy = bool.Parse(Attributes2[4]);
                    gp.Sequence.nGlycoSite = Attributes2[5] ; 
                    gp.SequenceMonoMass = double.Parse(Attributes2[6]) ; 
                    gp.SequenceAverageMass = double.Parse(Attributes2[7]) ; 
                    gp.Glycan.SetMonosaccharideCompostion(Attributes2[8]) ;
                    gp.Glycan.composition = Attributes2[8]; 
                    gp.Glycan.is_decoy = bool.Parse(Attributes2[9]) ; 
                    gp.GlycanMonoMass = double.Parse(Attributes2[10]) ; 
                    gp.GlycanAverageMass = double.Parse(Attributes2[11]) ;

                    gp.IsDecoy = gp.Sequence.is_decoy | gp.Glycan.is_decoy;

                    if (create_hash)
                    {
                        List<GlycopeptideRecord> value_gp_list = new List<GlycopeptideRecord>();
                        int key_mass = (int)Math.Floor(gp.GP_Mono_Mass);

                        int min_mass_int = (int)Math.Floor(min_mass);
                        int max_mass_int = (int)Math.Floor(max_mass);

                        if (!_dGlycopeptides.ContainsKey(key_mass))
                        {
                            if (value_gp_list.Count > 0)
                                value_gp_list.Clear();
                            value_gp_list.Add(gp);
                            _dGlycopeptides.Add(key_mass, value_gp_list);
                        }
                        else
                        {
                            value_gp_list.Clear();
                            value_gp_list = _dGlycopeptides[key_mass];
                            value_gp_list.Add(gp);
                            _dGlycopeptides[key_mass] = value_gp_list;

                        }
                    }
                    else
                        _glycopeptides.Add(gp);                    
                }
            }
            CsvFileHandler2.closeFile(); 
         }


         public void LoadFasta(string fastaFile, List<String> peptides, bool set_decoy)
         {
             GlypID.Sequence.clsSequence[] sequences = new GlypID.Sequence.clsSequence[1]; 
             GlypID.Glycopeptide.clsGlycopeptide GlycoPeptide = new GlypID.Glycopeptide.clsGlycopeptide();
             GlycoPeptide.LoadNGlycopeptidesFromFasta(fastaFile, ref sequences, set_decoy);
             for (int i = 0; i < sequences.Length; i++)
             {
                 peptides.Add(sequences[i].sequence); 
             }

         }

         public void Reset()
        {
            if (_ParentTransformResults.Length > 0)           
                Array.Clear(_ParentTransformResults, 0, _ParentTransformResults.Length);
            if (_ParentPeaks.Length > 0)
                Array.Clear(_ParentPeaks, 0, _ParentPeaks.Length);
            if (_MSMSPeaks.Length > 0)
                Array.Clear(_MSMSPeaks, 0, _MSMSPeaks.Length); 
             
        }

        /// <summary>
        ///  Function to store glycopeptides
        /// </summary>
        /// <param name="sequences"> Peptide sequences</param>
        /// <param name="glycans"> Glycans</param>
        /// <param name="min_mass">Min Mass of the map</param>
        /// <param name="max_mass">Max mass of the map</param>
        public void SetGlycoPeptides(ref GlypID.Sequence.clsSequence[] sequences, ref GlypID.Glycan.clsGlycan[] glycans)
        {
            int num_gps = sequences.Length * glycans.Length; 
            _glycopeptides = new List<GlycopeptideRecord>();    
       
             
            for (int i = 0; i < sequences.Length; i++)
            {
                for (int j=0 ; j < glycans.Length; j++)
                {  
                    GlycopeptideRecord gp = new GlycopeptideRecord();
                    
                    gp.Sequence = sequences[i];
                    gp.SequenceAverageMass = sequences[i].CalculateSequenceMass(false) ; 
                    gp.SequenceMonoMass = sequences[i].CalculateSequenceMass(true) ; 
                                        
                    gp.Glycan = glycans[j];
                    gp.GlycanAverageMass = glycans[j].CalculateGlycanMass(false) ; 
                    gp.GlycanMonoMass = glycans[j].CalculateGlycanMass(true) ;
                    
                    

                    gp.GP_Average_Mass = _utils.CalculateGPMass(gp.SequenceAverageMass, gp.GlycanAverageMass) ; 
                    gp.GP_Mono_Mass = _utils.CalculateGPMass(gp.SequenceMonoMass, gp.GlycanMonoMass) ;

                    gp.IsDecoy = gp.Glycan.is_decoy | gp.Sequence.is_decoy; 
                    
                   _glycopeptides.Add(gp); 
                }
            }
            //-- Sort by mono mass and create a hash table -- //
            _glycopeptides.Sort(delegate(GlycopeptideRecord g1, GlycopeptideRecord g2) { return g1.GP_Mono_Mass.CompareTo(g2.GP_Mono_Mass); });
           
        }
        /// <summary>
        /// Function to create a hash table/Dictionary format of all glycopeptides that makes searching faster
        /// </summary>
        public void CreateHashTableGlycopeptides(double min_mass, double max_mass, bool filter_mass)
        {

            _dGlycopeptides = new Dictionary<int, List<GlycopeptideRecord>>();

            foreach (GlycopeptideRecord g in _glycopeptides)
            {
                List<GlycopeptideRecord> value_gp_list = new List<GlycopeptideRecord>();
                int key_mass = (int)Math.Floor(g.GP_Mono_Mass);

                int min_mass_int = (int)Math.Floor(min_mass);
                int max_mass_int = (int)Math.Floor(max_mass);

                bool proceed = false; 
                if (filter_mass)
                {
                    if (key_mass >= min_mass_int && key_mass <= max_mass_int)                    
                        proceed = true; 
                }
                else
                    proceed = true; 

                if(proceed)
                {
                    if (!_dGlycopeptides.ContainsKey(key_mass))
                    {
                        if (value_gp_list.Count > 0)
                            value_gp_list.Clear();
                        value_gp_list.Add(g);
                        _dGlycopeptides.Add(key_mass, value_gp_list);
                    }
                    else
                    {
                        value_gp_list.Clear();
                        value_gp_list = _dGlycopeptides[key_mass];
                        value_gp_list.Add(g);
                        _dGlycopeptides[key_mass] = value_gp_list;

                    }
                }

            }
        }


        /// <summary>
        /// Function to write out glycopeptides dictionary to file so that next time
        /// </summary>
        /// <param name="outfile"></param>
        public void WriteOutGlycopeptidesToFile(string filename)
        {
            CSVFileHandler outfile = new CSVFileHandler(filename, CSVFileHandler.WRITE_ONLY);
            outfile.openFile();

            int numglycopeptides = _glycopeptides.Count;
            string[] header = { "GlycoPeptideMonoMass", "GlycoPeptideAvgMass", "Protein", "Peptide", "DECOY_Peptide",
                                  "Site", "PeptideMonoMass", "PeptideAvgMass", "Glycan", "DECOY_Glycan", "GlycanMonoMass", "GlycanAvgMass" };
            outfile.writeLine(header) ;
            for (int i = 0; i < numglycopeptides; i++)
            {
                string[] outline = {Convert.ToString(_glycopeptides[i].GP_Mono_Mass), 
                                       Convert.ToString(_glycopeptides[i].GP_Average_Mass), 
                                       Convert.ToString(_glycopeptides[i].Sequence.proteinName), 
                                       Convert.ToString(_glycopeptides[i].Sequence.sequence),
                                       Convert.ToString(_glycopeptides[i].Sequence.is_decoy),
                                       Convert.ToString(_glycopeptides[i].Sequence.nGlycoSite),
                                       Convert.ToString(_glycopeptides[i].SequenceMonoMass) , 
                                       Convert.ToString(_glycopeptides[i].SequenceAverageMass), 
                                       Convert.ToString(_glycopeptides[i].Glycan.composition) ,
                                       Convert.ToString(_glycopeptides[i].Glycan.is_decoy), 
                                       Convert.ToString(_glycopeptides[i].GlycanMonoMass), 
                                       Convert.ToString(_glycopeptides[i].GlycanAverageMass)};

                outfile.writeLine(outline);
            }

            outfile.closeFile();           
        }





        /// <summary>
        /// Search a particular glycopeptide by mass (the dictionary version, which is faster)
        /// </summary>
        /// <param name="mass">Mass to be searched for</param>
        /// <param name="matched_gps"> list of matched glycopeptides</param>
        /// <param name="check_only_sialylated">Simple filtering criterion</param>
        /// <returns>True/False of a found match</returns>
        public bool SearchGlycoPeptidesDictionary(double mass, ref  List<GlycopeptideRecord> matched_gps, bool check_only_sialylated)
        {
            bool found_match = false;
            matched_gps = new List<GlycopeptideRecord>();
            double min_ppm = _PPM_Diff;

            int gpmass = (int)Math.Floor(mass);
            int[] key_mass = { gpmass - 1, gpmass, gpmass + 1 };
            

            for (int i = 0; i < key_mass.Length; i++)
            {
                List<GlycopeptideRecord> value_gp_list = new List<GlycopeptideRecord>();
                if (_dGlycopeptides.ContainsKey(key_mass[i]))
                {
                    value_gp_list = _dGlycopeptides[key_mass[i]];
                    foreach (GlycopeptideRecord g in value_gp_list)
                    {
                        bool proceed = false;
                        if (check_only_sialylated)
                        {
                            string comp = g.Glycan.composition;
                            if (g.Glycan.GlycanCompositionHasNeuAC(comp))
                                proceed = true;
                        }
                        else
                            proceed = true;

                       
                        if (proceed)
                        {
                            if (_utils.CalculateDelMassPPM(g.GP_Mono_Mass, mass) < min_ppm)
                            {
                                found_match = true;
                                matched_gps.Add(g);

                                
                            }
                        }
                    }
                }
            }

            return found_match;
        }

        /// <summary>
        /// Search all glycoeptide records for a particular mass
        /// </summary>
        /// <param name="mass">Mass to be searched for</param>
        /// <param name="matched_gps">List of GPs that matched</param>
        /// <param name="check_only_sialylated">Simple filtering rule</param>
        /// <returns>TRUE/FALSE found at least one match</returns>
        /* deperecated nov 2012 public bool SearchGlycopeptides(double mass, ref  List <GlycopeptideRecord> matched_gps, bool check_only_sialylated)
        {
            bool found_match = false;
            matched_gps = new List<GlycopeptideRecord>();
            double min_ppm = _PPM_Diff; 
            foreach (GlycopeptideRecord g in _glycopeptides)
            {
                bool proceed = false;
                if (check_only_sialylated)
                {
                    string comp = g.Glycan.composition;
                    if (g.Glycan.GlycanCompositionHasNeuAC(comp))
                        proceed = true;
                }
                else
                    proceed = true;           

                if (proceed)
                {
                    if (_utils.CalculateDelMassPPM(g.GP_Mass, mass) < min_ppm)
                    {
                        found_match = true;
                        matched_gps.Add(g);
                    }
                }
                if (g.GP_Mass > mass + 0.5)
                    break; 

            }

            return found_match; 
        }*/


    

       

        /// <summary>
        ///  Function to Load all fragmentation spectra for glycopeptide ions into the map.
        /// </summary>
        /// <param name="_glycoMap"> The Map where stuff gets loaded into</param>
        /// <param name="filterGlycoPeps">Select only glycopeptides, usually set to true</param>
        /// <param name="glycanListFile">The glycan file</param>
        /// <param name="fastaFile">The FASTA file</param>
        /// <param name="gpFile">A file for glycopeptides</param>
        public void LoadFragmentationSpectraIntoMap(ref Classes.MapRecord _glycoMap, ref Classes.Params _parameters, string glycanListFile, string fastaFile, string gpFile )
        {            
            List<string> DatasetNames = _glycoMap._AssociatedDatasetNames; 
            List <bool> isCID = _glycoMap._IsCID;
            List <bool> isHCD = _glycoMap._IsHCD;
            List <bool> isETD = _glycoMap._IsETD;

            bool filter_glycopeptides = _parameters.ProcessOnlyNGlycopeptides;
            double max_coverage = _parameters.MaxUMCCoverage;

            GlypID.Sequence.clsSequence[] sequences = new GlypID.Sequence.clsSequence[1];
            GlypID.Glycan.clsGlycan[] glycans = new GlypID.Glycan.clsGlycan[1];
            GlypID.Glycopeptide.clsGlycopeptide GlycoPeptide = new GlypID.Glycopeptide.clsGlycopeptide();
            Classes.MapRecord _tMap = new MapRecord(); 

            // -- Get Glycopeptide List --//
            if (_parameters.ProcessOnlyNGlycopeptides)
            {
                if (!_parameters.UseGlycoPeptideFile)
                {
                    Console.WriteLine("Reading Glycan List");
                    try
                    {
                        GlycoPeptide.LoadGlycansFromList(glycanListFile, ref glycans, _parameters.UseDecoyGlycan);
                    }
                    catch (Exception e)
                    {
                        System.Console.WriteLine(e.Message);
                        return; 
                    }
                    Console.WriteLine("Reading Fasta File"); 
                    GlycoPeptide.LoadNGlycopeptidesFromFasta(fastaFile, ref sequences, _parameters.UseDecoyPeptide);
                    Console.WriteLine("Setting Glycopeptides"); 
                    SetGlycoPeptides(ref sequences, ref glycans);  //slightly inflated*/    

                    if (_parameters.CreateGlycoPeptideFile)
                    {
                        int last_pos2 = fastaFile.LastIndexOf(".");
                        string outputgpfile = fastaFile.Substring(0, last_pos2) + "_GPInfo_v2.csv";
                        Console.WriteLine("Writing to GP File");
                        WriteOutGlycopeptidesToFile(outputgpfile);
                    }
                    Console.WriteLine("Creating Hash"); 
                    CreateHashTableGlycopeptides(_glycoMap.MapMinMass - 1, _glycoMap.MapMaxMass + 1, true); 
                }
                else
                {
                    Console.WriteLine("Reading in from GP file");
                    LoadGlycopeptidesFromFile(gpFile , _glycoMap.MapMinMass - 1, _glycoMap.MapMaxMass + 1, true );
                    Console.WriteLine("Creating Hash") ; 
                    //CreateHashTableGlycopeptides(_glycoMap.MapMinMass - 1, _glycoMap.MapMaxMass + 1, false); 

                }
            }

            Console.WriteLine("Loading Frag events and scoring"); 
            _glycopeptides.Clear(); 
            
            _glycoMap._AllMLNRecords.ForEach(delegate(MultiAlignRecord m)
            {
                bool process_mrecord = false; 
                Console.WriteLine("UMC = " + m.ID) ;

                 
                  
                // Only load in thosse ions that match the glycopeptides list

                if (filter_glycopeptides)
                {
                    List<GlycopeptideRecord> candidate_gps = new List<GlycopeptideRecord>();
                    process_mrecord =  SearchGlycoPeptidesDictionary(m.Mass, ref candidate_gps, false);

                    if (process_mrecord)
                    {
                        /*m._CandidatePeptideSeq = new string[candidate_gps.Count];
                        m._CandidateProteinName = new string[candidate_gps.Count];
                        m._CandidateGlycanComposition = new string[candidate_gps.Count];
                        m._CandidateGlycanMass = new double[candidate_gps.Count];
                        m._CandidatePeptideMass = new double[candidate_gps.Count];
                        m._CandidateNGlycoSite = new string[candidate_gps.Count];
                        for (int i = 0; i < candidate_gps.Count; i++)
                        {
                            m._CandidateGlycanComposition[i] = candidate_gps[i].Glycan.composition;
                            m._CandidatePeptideSeq[i] = candidate_gps[i].Sequence.sequence;
                            m._CandidateProteinName[i] = candidate_gps[i].Sequence.proteinName;
                            m._CandidateGlycanMass[i] = candidate_gps[i].GlycanMonoMass;
                            m._CandidatePeptideMass[i] = candidate_gps[i].SequenceMonoMass; ;
                            m._CandidateNGlycoSite[i] = candidate_gps[i].Sequence.nGlycoSite;



                           
                        }*/
                        m._CandidateGlycopeptideRecords = new GlycopeptideRecord[candidate_gps.Count];
                        for (int i = 0; i < candidate_gps.Count; i++)
                        {
                            m._CandidateGlycopeptideRecords[i] = candidate_gps[i]; 
                        }
                    }

                }
                else
                {
                    process_mrecord = true;
                }


                if (process_mrecord)
                {
                    bool store_record = false ; // This determins if the record needs to be stored 
                    for (int i = 0; i < m._AssociatedUMCRecords.Count; i++)
                    {
                        // -- Start processing each umc -- //                       
                        UMCRecord u = m._AssociatedUMCRecords[i];
                        if (u.ScanRep == 0)
                            continue;

                     

                        int datasetid = u.DatasetID;
                       // _RawData = new GlypID.Readers.clsRawData(DatasetNames[datasetid], GlypID.Readers.FileType.FINNIGAN);
                        

                        double scanrange = Convert.ToDouble(u.ScanEnd - u.ScanStart);
                        int numscans = _RawData.GetNumScans();
                        double coverage = scanrange / numscans;
                        if (coverage >= max_coverage)
                        {
                            //_RawData.
                            continue;
                        }

                        // umc has passed filtering conditions so start processing each scan in umc               
                        double min_ppm = 50;
                        

                        // Process each MSn scan from start to stop of UMC
                        List<double> parents_observed = new List<double>();
                        List<FragEvents>allobserved_frag_events = new List<FragEvents>() ; 
                        Classes.FragEvents e = new FragEvents();

                        double most_recent_precursor = 0;
                        int frag_id = 0; 
                        for (int scan = u.ScanStart ; scan <= u.ScanEnd ; scan++)
                        {
                            if (_RawData.IsMSScan(scan))
                                continue;                            
                                                 
                            // Start processing MSn scan                            
                            bool process_scan = false;
                            bool record_scan = false ; 

                            // get parent                               
                            double parent_mz = _RawData.GetParentMz(scan);                            
                            short header_cs = (short) _RawData.GetParentChargeFromHeader(scan);
                            
                            if (header_cs > 0)
                            {
                                double mass = _utils.CalculateMass(parent_mz, header_cs);
                                if (_utils.CalculateDelMassDa(mass, m.Mass) < _DA_Diff)
                                    process_scan = true;
                            }
                            else
                            {
                                // try default charge states
                                double mass1 = _utils.CalculateMass(parent_mz, 2);
                                if (_utils.CalculateDelMassDa(mass1, m.Mass) < _DA_Diff)
                                    process_scan = true;
                                double mass2 = _utils.CalculateMass(parent_mz, 3);
                                if (_utils.CalculateDelMassDa(mass2, m.Mass) < _DA_Diff)
                                    process_scan = true;
                                double mass3 = _utils.CalculateMass(parent_mz, 4);
                                if (_utils.CalculateDelMassDa(mass1, m.Mass) < _DA_Diff)
                                    process_scan = true;
                                double mass4 = _utils.CalculateMass(parent_mz, 5);
                                if (_utils.CalculateDelMassDa(mass2, m.Mass) < _DA_Diff)
                                    process_scan = true;
                                double mass5 = _utils.CalculateMass(parent_mz, 6);
                                if (_utils.CalculateDelMassDa(mass1, m.Mass) < _DA_Diff)
                                    process_scan = true;
                                double mass6 = _utils.CalculateMass(parent_mz, 7);
                                if (_utils.CalculateDelMassDa(mass2, m.Mass) < _DA_Diff)
                                    process_scan = true;
                            }
                            if (!process_scan)
                                continue;

                            // Scan has passed ALL filters                            
                            float[] msms_mzs = new float[1];
                            float[] msms_intensities = new float[1];
                            _RawData.GetRawData(scan, ref msms_mzs, ref msms_intensities);
                            if (_RawData.IsProfileScan(scan))
                                _MSMSPeakProcessor.ProfileType = GlypID.enmProfileType.PROFILE;
                            else
                                _MSMSPeakProcessor.ProfileType = GlypID.enmProfileType.CENTROIDED;

                            if (most_recent_precursor != parent_mz)
                            {
                                // Only deisotope if parent has not been observed

                                // Store the previous Fragmentation event
                                if (e.HCDScore < 1)
                                {
                                    /*FragEvents tEvent = new FragEvents();
                                    tEvent = e; 
                                    m._AssociatedUMCRecords[i]._AssociatedFragEvents.Add(tEvent);*/
                                    allobserved_frag_events.Add(e);
                                }

                                // Clear up everything 
                               /* e.ClearTransformRecord();
                                e.ClearHCD();
                                e.ClearETD();
                                e.ClearCID();
                                e.ClearGPInfo();*/
                                e = new FragEvents();
                                e.ID = frag_id;
                                frag_id++; 

                                // Start over
                                parents_observed.Add(parent_mz);
                                most_recent_precursor = parent_mz; 
                                float[] parent_mzs = new float[1];
                                float[] parent_intensities = new float[1];
                                Reset();

                                // ** ---- deisotope the precursor -- ** //
                                int parent_scan = _RawData.GetParentScan(scan);
                                double parent_scan_time = _RawData.GetScanTime(parent_scan);  //TODO
                                _RawData.GetRawData(parent_scan, ref parent_mzs, ref parent_intensities);

                                // Do peak finding
                                if (_RawData.IsProfileScan(parent_scan))
                                    _ParentPeakProcessor.ProfileType = GlypID.enmProfileType.PROFILE;
                                else
                                    _ParentPeakProcessor.ProfileType = GlypID.enmProfileType.CENTROIDED;
                                double thresh = GlypID.Utils.GetAverage(ref parent_intensities, float.MaxValue);
                                double background_intensity = GlypID.Utils.GetAverage(ref parent_intensities, (float)(5 * thresh));
                                _ParentPeakProcessor.SetPeakIntensityThreshold(background_intensity * _ParentPeakProcessorParams.PeakBackgroundRatio);
                                _ParentPeakProcessor.DiscoverPeaks(ref parent_mzs, ref parent_intensities, ref _ParentPeaks,
                                    Convert.ToSingle(_TransformParameters.MinMZ), Convert.ToSingle(_TransformParameters.MaxMZ), true);

                                // Pick out precursor and transform
                                double pep_intensity;
                                if (_TransformParameters.UseAbsolutePeptideIntensity)
                                    pep_intensity = _TransformParameters.AbsolutePeptideIntensity;
                                else
                                    pep_intensity = background_intensity * _TransformParameters.PeptideMinBackgroundRatio;
                                bool found = false;
                                if (_RawData.IsFTScan(parent_scan))
                                {
                                    found = _HornTransform.FindPrecursorTransform(Convert.ToSingle(background_intensity), Convert.ToSingle(pep_intensity), ref parent_mzs, ref parent_intensities,
                                        ref _ParentPeaks, Convert.ToSingle(parent_mz), ref _ParentTransformResults);
                                }
                                if (!found)
                                {
                                    // Low resolution data or bad high res spectra
                                    if (header_cs > 0)
                                    {
                                        double mono_mz = _RawData.GetParentMonoMzFromHeader(scan);
                                        if (mono_mz == 0)
                                            mono_mz = parent_mz;
                                        /*   GlypID.Peaks.clsPeak monoPeak = new GlypID.Peaks.clsPeak();
                                           _ParentPeakProcessor.GetClosestPeakMz(monoPeak, Convert.ToSingle(mono_mz));*/

                                        short[] charges = new short[1];
                                        charges[0] = header_cs;
                                        _HornTransform.AllocateValuesToTransform(Convert.ToSingle(mono_mz), 0, ref charges, ref _ParentTransformResults);
                                        found = true;
                                        record_scan = true;
                                        e.TransformResult = _ParentTransformResults[0];
                                        e.ParentMz = parent_mz;
                                        e.ParentScan = parent_scan;
                                        e.ParentScanTime = parent_scan_time; 
                                    }
                                    /*else Removind this for now
                                    {
                                        // instrument has no charge just store 2 and 3.      
                                        short[] charges = new short[2];
                                        charges[0] = 2;
                                        charges[1] = 3;
                                        _HornTransform.AllocateValuesToTransform(Convert.ToSingle(parent_mz), 0, ref charges, ref _ParentTransformResults);
                                    }*/
                                }
                                else
                                {
                                    double ppm = _utils.CalculateDelMassPPM(m.Mass, _ParentTransformResults[0].mdbl_mono_mw);
                                    if ((ppm <= min_ppm) || (Math.Abs(ppm - min_ppm) < 1E-05))
                                    {
                                        record_scan = true;
                                        e.TransformResult = _ParentTransformResults[0];
                                        e.ParentMz = parent_mz;
                                        e.ParentScan = parent_scan;
                                        e.ParentScanTime = parent_scan_time; 
                                    }                                    
                                }                                       
                            }
                            else
                            {
                                if (e.TransformResult.mdbl_mono_mw > 0)
                                {
                                    // Has been deisotoped successfulle
                                    record_scan = true ; 
                                }
                            }

                            

                            //  -------- Process scan according to type   -------- //
                            if (isCID[datasetid] && _RawData.IsCIDScan(scan) && record_scan)
                            {
                                // CID scan                                         
                                _MSMSPeakProcessor.SetPeakIntensityThreshold(0); // Since yin takes top 20 on his own.
                                _MSMSPeakProcessor.DiscoverPeaks(ref msms_mzs, ref msms_intensities, ref _MSMSPeaks, Convert.ToSingle(_TransformParameters.MinMZ), Convert.ToSingle(_TransformParameters.MaxMZ), false);                             
                                if (_MSMSPeaks.Length > 0)
                                {                                   
                                    //Scoring                                 
                                    GlypID.CIDScoring.clsCIDScoringScanResults[] tCIDScoreResults = new GlypID.CIDScoring.clsCIDScoringScanResults[1];
                                    bool found_score = _CIDScoring.ScoreCIDSpectra(ref _MSMSPeaks, ref msms_mzs, ref msms_intensities, ref _ParentTransformResults, ref tCIDScoreResults);
                                    if (found_score & tCIDScoreResults[0].mdbl_cid_score > e.CIDScore)
                                    {
                                      
                                        e.CIDScan = scan;
                                        e.CIDMzs = msms_mzs;
                                        e.CIDIntensities = msms_intensities;
                                        e.CIDPeaks = new GlypID.Peaks.clsPeak[_MSMSPeaks.Length];
                                        e.CIDProfileType = _MSMSPeakProcessor.ProfileType;
                                        e.CIDScore = tCIDScoreResults[0].mdbl_cid_score;
                                        Array.Copy(_MSMSPeaks, e.CIDPeaks, _MSMSPeaks.Length);
                                    }
                                }
                            }

                            if (isHCD[datasetid] && _RawData.IsHCDScan(scan) && record_scan)
                            {
                                // HCD Scan
                                double hcd_background_intensity = GlypID.Utils.GetAverage(ref msms_intensities, ref msms_mzs, Convert.ToSingle(_ScoringParameters.MinHCDMz), Convert.ToSingle(_ScoringParameters.MaxHCDMz));
                                _MSMSPeakProcessor.SetPeakIntensityThreshold(hcd_background_intensity);
                                _MSMSPeakProcessor.DiscoverPeaks(ref msms_mzs, ref msms_intensities, ref _MSMSPeaks, Convert.ToSingle(_ScoringParameters.MinHCDMz), Convert.ToSingle(_ScoringParameters.MaxHCDMz), false);                               
                                if (_MSMSPeaks.Length > 0)
                                {
                                    // Score
                                    GlypID.HCDScoring.clsHCDScoringScanResults[] tHCDScoreResults = new GlypID.HCDScoring.clsHCDScoringScanResults[1];
                                    double score = _HCDScoring.ScoreHCDSpectra(ref _MSMSPeaks, ref msms_mzs, ref msms_intensities, ref _ParentTransformResults, ref tHCDScoreResults);
                                    if (score < e.HCDScore) // This makes sure that within frag events of the same parent, the lowest one is chosen
                                    {        
                                        // Store
                                        store_record = true;
                                       
                                        e.HCDScan = scan;                                         
                                        e.HCDMzs = msms_mzs;
                                        e.HCDIntensities = msms_intensities;
                                        e.HCDPeaks = new GlypID.Peaks.clsPeak[_MSMSPeaks.Length];
                                        /*m._AssociatedUMCRecords[i].HCDScanRep = scan;
                                        m._AssociatedUMCRecords[i].HCDScanRepET = _RawData.GetScanTime(scan);*/
                                        e.HCDProfileType = _MSMSPeakProcessor.ProfileType;
                                        Array.Copy(_MSMSPeaks, e.HCDPeaks, _MSMSPeaks.Length);
                                        e.HCDScore = tHCDScoreResults[0].mdbl_hcd_score; 
                                        e.GlycanType = (GlypID.enmGlycanType) tHCDScoreResults[0].menm_glycan_type ;
                                        
                                    }
                                }
                            }

                            if (isETD[datasetid] && _RawData.IsETDScan(scan) && record_scan)
                            {

                                // ETD Scan
                                _MSMSPeakProcessor.SetPeakIntensityThreshold(0); 
                                _MSMSPeakProcessor.ProfileType = GlypID.enmProfileType.CENTROIDED; 
                                _MSMSPeakProcessor.DiscoverPeaks(ref msms_mzs, ref msms_intensities, ref _MSMSPeaks, 
                                    0, Convert.ToSingle(parent_mz), false);
                                                                                               
                                if ((_MSMSPeaks.Length > 1) && (e.HCDScore <1))                                
                                {                                   
                                    // Score
                                    int best_scoring_index = -1;
                                    double max_score = 0;
                                    
                                    for (int k = 0; k < m._CandidateGlycopeptideRecords.Length; k++)
                                    {
                                        // Confirm glycan type matches with HCD prediction
                                        bool process_result = false;
                                        if ((e.GlycanType == GlypID.enmGlycanType.CS) || e.GlycanType == GlypID.enmGlycanType.HY)
                                        {                                            
                                            if (m._CandidateGlycopeptideRecords[k].Glycan.numNeuAc >0)
                                                process_result = true;
                                            else
                                                process_result = false;
                                        }
                                        else
                                            process_result = true;


                                        // Stupidly planned type conversion ; Just doing the bare minimum
                                        GlypID.ETDScoring.clsETDScoringScanResults thisResult = new GlypID.ETDScoring.clsETDScoringScanResults();
                                        thisResult.mdbl_parent_mz = e.ParentMz;
                                        thisResult.mdbl_mono_mw = e.TransformResult.mdbl_mono_mw;
                                        thisResult.mshort_cs = e.TransformResult.mshort_cs;
                                        thisResult.mstr_glyco_site = m._CandidateGlycopeptideRecords[k].Sequence.nGlycoSite;
                                        thisResult.mstr_nglyco_site = m._CandidateGlycopeptideRecords[k].Sequence.nGlycoSite;
                                        thisResult.mstr_pep_seq = m._CandidateGlycopeptideRecords[k].Sequence.sequence;
                                        thisResult.mdbl_glycan_mass = m._CandidateGlycopeptideRecords[k].GlycanMonoMass;
                                       
                                        if (thisResult.mstr_pep_seq != "" && process_result)
                                        {                                            
                                            double etd_score = _ETDScoring.ScoreETDSpectra(ref _MSMSPeaks, thisResult);
                                            if (etd_score >= max_score)
                                            {
                                                thisResult.mdbl_etd_score = etd_score;
                                                best_scoring_index = k;
                                                max_score = etd_score; 
                                            }
                                        }
                                    }

                                    // Store
                                    if ((best_scoring_index != -1) && (max_score > e.ETDScore))
                                    {
                                        // Store     
                                        e.ETDScan = scan; 
                                        e.ETDMzs = msms_mzs;
                                        e.ETDIntensities = msms_intensities;
                                        e.ETDPeaks = new GlypID.Peaks.clsPeak[_MSMSPeaks.Length];
                                        Array.Copy(_MSMSPeaks, e.ETDPeaks, _MSMSPeaks.Length);

                                        e.ETDScore = max_score;
                                        e.GP_Record = m._CandidateGlycopeptideRecords[best_scoring_index];
                                        e.FalseHit = e.GP_Record.IsDecoy; 

                                       /* if (e.FalseHit)
                                            _tMap._AllFalseHitsFDRScore.Add(e.ETDScore);
                                        else
                                            _tMap._AllTrueHitsFDRScore.Add(e.ETDScore) ; */
                                       
                                        // m._AssociatedUMCRecords[i].ETDScanRep = scan;
                                        
                                    }
                                }                                
                            }
                        }
                        // Store the most recent fragmentation event
                        if (e.HCDScore < 1)
                        {
                            /*FragEvents tEvent = new FragEvents(e); 
                            m._AssociatedUMCRecords[i]._AssociatedFragEvents.Add(tEvent);*/
                            allobserved_frag_events.Add(e); 
                        }
                       //  _RawData.Close();


                        m._AssociatedUMCRecords[i]._AssociatedFragEvents.Clear(); 
                        for (int kk = 0; kk < allobserved_frag_events.Count; kk++)
                        {
                            FragEvents tempE = allobserved_frag_events[kk] ; //as FragEvents;
                            m._AssociatedUMCRecords[i]._AssociatedFragEvents.Add(tempE); 
                        }
                        m._AssociatedUMCRecords[i].SetRepScores(); 
                    }
                    if (store_record)                    
                        _tMap.AddRecord(m);
                    
                }
                
            });

            if (_tMap._AllMLNRecords.Count != _glycoMap._AllMLNRecords.Count)
            {
                _glycoMap.ClearRecords();
                for (int k = 0; k < _tMap._AllMLNRecords.Count; k++)
                {
                    _glycoMap.AddRecord(_tMap._AllMLNRecords[k]);
                }
            }
           // _glycoMap._AllTrueHitsFDRScore = _tMap._AllTrueHitsFDRScore;
            //_glycoMap._AllFalseHitsFDRScore = _tMap._AllFalseHitsFDRScore; 

            _tMap.ClearRecords(); 
        }


        public void LoadPrecursorInfoIntoMap(ref Classes.MapRecord _glycoMap, List<int> conflictedIDs)
        {
            /*List<string> DatasetNames = _glycoMap._AssociatedDatasetNames;

            Classes.MapRecord _tMap = new MapRecord();

            _glycoMap._AllMLNRecords.ForEach(delegate(MultiAlignRecord m)
            {
                bool exists = conflictedIDs.Exists(element => element == m.ID);
                if (exists)
                {
                    Console.WriteLine("UMC = " + m.ID);
                    for (int i = 0; i < m._AssociatedUMCRecords.Count; i++)
                    {
                        UMCRecord u = m._AssociatedUMCRecords[i];
                        if (u.ScanRep > 0)
                        {
                            int datasetid = u.DatasetID;
                            _RawData = new GlypID.Readers.clsRawData(DatasetNames[datasetid], GlypID.Readers.FileType.FINNIGAN);

                            double scanrange = Convert.ToDouble(u.ScanEnd - u.ScanStart);
                            int numscans = _RawData.GetNumScans();
                            double min_ppm = _ScoringParameters.PPMTolerance;
                            for (int scan = u.ScanStart; scan <= u.ScanEnd; scan++)
                            {
                                if (!_RawData.IsMSScan(scan))
                                {
                                    float[] parent_mzs = new float[1];
                                    float[] parent_intensities = new float[1];
                                    float[] msms_mzs = new float[1];
                                    float[] msms_intensities = new float[1];

                                    // Not an MS scan, 
                                    bool process_scan = false;
                                    bool record_scan = false;

                                    // get parent                                
                                    double parent_mz = _RawData.GetParentMz(scan);
                                    // Do a quick check if parent mass is somewhat close to multi-align calibrated mass
                                    short header_cs = _RawData.GetMonoChargeFromHeader(scan);
                                    if (header_cs > 0)
                                    {
                                        double mass = _utils.CalculateMass(parent_mz, header_cs);
                                        if (_utils.CalculateDelMassDa(mass, m.Mass) < _DA_Diff)
                                            process_scan = true;
                                    }
                                    else
                                    {
                                        // try default charge states
                                        double mass1 = _utils.CalculateMass(parent_mz, 2);
                                        if (_utils.CalculateDelMassDa(mass1, m.Mass) < _DA_Diff)
                                            process_scan = true;
                                        double mass2 = _utils.CalculateMass(parent_mz, 3);
                                        if (_utils.CalculateDelMassDa(mass2, m.Mass) < _DA_Diff)
                                            process_scan = true;

                                    }

                                    if (process_scan)
                                    {
                                        Reset();

                                        // deisotope the precursor
                                        int parent_scan = _RawData.GetParentScan(scan);
                                        _RawData.GetSpectrum(scan, ref msms_mzs, ref msms_intensities);
                                        _RawData.GetSpectrum(parent_scan, ref parent_mzs, ref parent_intensities);
                                        // Do peak finding
                                        _ParentPeakProcessor.DiscoverPeaks(ref parent_mzs, ref parent_intensities, ref _ParentPeaks,
                                            Convert.ToSingle(_TransformParameters.MinMZ), Convert.ToSingle(_TransformParameters.MaxMZ), true);
                                        double bk_intensity = 0; // _ParentPeakProcessor.GetBackgroundIntensity(ref parent_intensities);
                                        double pep_intensity;
                                        if (_TransformParameters.UseAbsolutePeptideIntensity)
                                            pep_intensity = _TransformParameters.AbsolutePeptideIntensity;
                                        else
                                            pep_intensity = 0; // bk_intensity * _TransformParameters.PeptideMinBackgroundRatio;
                                        // deisotope precursor
                                        bool found = false;
                                        if (_RawData.IsFTScan(parent_scan))
                                        {
                                            found = _HornTransform.FindPrecursorTransform(Convert.ToSingle(bk_intensity), Convert.ToSingle(pep_intensity), ref parent_mzs, ref parent_intensities,
                                                ref _ParentPeaks, Convert.ToSingle(parent_mz), ref _ParentTransformResults);
                                        }
                                        if (!found)
                                        {
                                            // Low resolution data or bad high res spectra
                                            short cs = _RawData.GetMonoChargeFromHeader(scan);
                                            if (cs > 0)
                                            {
                                                short[] charges = new short[1];
                                                charges[0] = cs; //NEED TO FIX
                                               // _HornTransform.AllocateValuesToTransform(Convert.ToSingle(parent_mz), ref charges, ref _ParentTransformResults);
                                            }
                                            else
                                            {
                                                // instrument has no charge just store 2 and 3.      
                                                short[] charges = new short[2];
                                                charges[0] = 2;
                                                charges[1] = 3;
                                               // _HornTransform.AllocateValuesToTransform(Convert.ToSingle(parent_mz), ref charges, ref _ParentTransformResults);
                                            }
                                        }

                                        for (int k = 0; k < _ParentTransformResults.Length; k++)
                                        {
                                            double ppm = _utils.CalculateDelMassPPM(m.Mass, _ParentTransformResults[k].mdbl_mono_mw);
                                            if ((ppm <= min_ppm) || (Math.Abs(ppm - min_ppm) < 1E-05))
                                            {
                                                min_ppm = ppm;
                                                record_scan = true;
                                                m._AssociatedUMCRecords[i].ClearTransformRecord();
                                                m._AssociatedUMCRecords[i].TransformResult = _ParentTransformResults[k];
                                                m._AssociatedUMCRecords[i].ClassRepMz = _ParentTransformResults[k].mdbl_mz;
                                                m._AssociatedUMCRecords[i].ClassRepCS = _ParentTransformResults[k].mshort_cs;
                                                m._AssociatedUMCRecords[i].ClassRepMW = _ParentTransformResults[k].mdbl_mono_mw; 
                                                m._AssociatedUMCRecords[i].ElutionTime = _RawData.GetScanTime(parent_scan);
                                            }
                                        }
                                    }
                                }
                            }
                            _RawData.Close();
                        }
                    }
                    _tMap.AddRecord(m);
                }
            });
            
            _glycoMap.ClearRecords();
            for (int k = 0; k < _tMap._AllMLNRecords.Count; k++)
            {
                _glycoMap.AddRecord(_tMap._AllMLNRecords[k]);
            }*/
             
        }
    }
}
