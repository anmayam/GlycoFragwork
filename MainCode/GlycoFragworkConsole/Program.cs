using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using GlycoFragworkDLL;
using System.Xml; 

namespace GlycoFragworkConsole
{
    class Program
    {
        // Default Glycan File
        private static string glycanFile; // =  @"C:\development\GlycoFragwork\GlycoFragworkDLL\DLLs\Default_Combination_V2.csv";
        private static string fastaFile; 
        private static string mlnfile;
        private static string gpFile; 
        private static bool useGpFile ;
        private static string nonGpFile; 

        private static List<string> datasetnames = new List<string>();
        private static List<bool> cidflags = new List<bool>();
        private static List<bool> hcdflags = new List<bool>();
        private static List<bool> etdflags = new List<bool>();
        private static List<string> datatype = new List<string>();

        

        private static void PrintUsage()
        {
            Console.WriteLine("Usage: GlycoFragworkConsole.exe -i inputfile.xml -p parameters.xml");
        }

        private static void LoadFileNames(XmlReader rdr)
        {

        }        

        /// <summary>
        /// Function to parse the input XML file
        /// </summary>
        /// <param name="file_name"></param>
        private static void ParseInputFile(string file_name)
        {
            XmlTextReader rdr = new XmlTextReader(file_name);
            string dataset = " ";
            bool cidflag = false;
            bool hcdflag = false;
            bool etdflag = false;
            string type = " ";
            while (rdr.Read())
            {             
                switch(rdr.NodeType)
                {
                    case XmlNodeType.Element:
                        if (rdr.Name == "MAPFile")
                        {
                            rdr.Read();
                            while (rdr.NodeType == XmlNodeType.Whitespace || rdr.NodeType == XmlNodeType.SignificantWhitespace)
                            {
                                rdr.Read();
                            }
                            if (rdr.NodeType != XmlNodeType.Text)
                                throw new Exception("Missing MultiAlignfile input");
                            mlnfile = System.Convert.ToString(rdr.Value);
                        }                          
                        if (rdr.Name == "GlycanFile")
                        {
                            rdr.Read();
                            while (rdr.NodeType == XmlNodeType.Whitespace || rdr.NodeType == XmlNodeType.SignificantWhitespace)
                            {
                                rdr.Read();
                            }
                            if (rdr.NodeType != XmlNodeType.Text)
                                throw new Exception("Missing Glycan file input");
                            glycanFile = System.Convert.ToString(rdr.Value);
                        }                       
                        if (rdr.Name == "GlycopeptideFile")
                        {    // Anoop : need to debug this before shipping      
                          
                            rdr.Read(); 
                            while(rdr.NodeType == XmlNodeType.Whitespace || rdr.NodeType == XmlNodeType.SignificantWhitespace)
                            {
                                rdr.Read();
                            }
                            if (rdr.NodeType != XmlNodeType.Text)
                                throw new Exception("Missing infromation on glycopeptide filename");
                            gpFile = Convert.ToString(rdr.Value); 
                        }
                        if (rdr.Name == "NonGlycopeptideFile")
                        {
                            rdr.Read();
                            while (rdr.NodeType == XmlNodeType.Whitespace || rdr.NodeType == XmlNodeType.SignificantWhitespace)
                            {
                                rdr.Read();
                            }
                            if (rdr.NodeType != XmlNodeType.Text)
                                throw new Exception("Missing infromation on non-glycopeptide filename");
                            nonGpFile = Convert.ToString(rdr.Value); 
                        }
                        if (rdr.Name == "FastaFile")
                        {
                            rdr.Read();
                            while (rdr.NodeType == XmlNodeType.Whitespace || rdr.NodeType == XmlNodeType.SignificantWhitespace)
                            {
                                rdr.Read();
                            }
                            if (rdr.NodeType != XmlNodeType.Text)
                                throw new Exception("Missing fasta file input");
                            fastaFile = System.Convert.ToString(rdr.Value);
                        }
                        if (rdr.Name == "RAW")
                        {
                            rdr.Read();
                            dataset = " ";
                            cidflag = false;
                            hcdflag = false;
                            etdflag = false;
                            type = " ";
                        }
                        if (rdr.Name == "File")
                        {
                            rdr.Read();
                            while (rdr.NodeType == XmlNodeType.Whitespace || rdr.NodeType == XmlNodeType.SignificantWhitespace)
                            {
                                rdr.Read();
                            }
                            if (rdr.NodeType != XmlNodeType.Text)
                                throw new Exception("Missing dataset input");
                            dataset = System.Convert.ToString(rdr.Value);
                        }
                        if (rdr.Name == "CID")
                        {
                            rdr.Read();
                            while (rdr.NodeType == XmlNodeType.Whitespace || rdr.NodeType == XmlNodeType.SignificantWhitespace)
                            {
                                rdr.Read();
                            }
                            if (rdr.NodeType != XmlNodeType.Text)
                                throw new Exception("Missing CID (true/false) input");
                            cidflag = System.Convert.ToBoolean(rdr.Value);
                        }                          
                        if (rdr.Name == "HCD")
                        {
                            rdr.Read();
                            while (rdr.NodeType == XmlNodeType.Whitespace || rdr.NodeType == XmlNodeType.SignificantWhitespace)
                            {
                                rdr.Read();
                            }
                            if (rdr.NodeType != XmlNodeType.Text)
                                throw new Exception("Missing HCD(true/false) input");
                            hcdflag = System.Convert.ToBoolean(rdr.Value);
                        }
                        if (rdr.Name == "ETD")
                        {
                            rdr.Read();
                            while (rdr.NodeType == XmlNodeType.Whitespace || rdr.NodeType == XmlNodeType.SignificantWhitespace)
                            {
                                rdr.Read();
                            }
                            if (rdr.NodeType != XmlNodeType.Text)
                                throw new Exception("Missing ETD (true/false) input");
                            etdflag = System.Convert.ToBoolean(rdr.Value);
                        }
                        if (rdr.Name == "Type")
                        {
                            rdr.Read();
                            while (rdr.NodeType == XmlNodeType.Whitespace || rdr.NodeType == XmlNodeType.SignificantWhitespace)
                            {
                                rdr.Read();
                            }
                            if (rdr.NodeType != XmlNodeType.Text)
                                throw new Exception("Missing Type (disease/control) input");
                            type = Convert.ToString(rdr.Value);
                        }
                        break;
                    case XmlNodeType.EndElement:
                        if (rdr.Name == "RAW")
                        {
                            datasetnames.Add(dataset);
                            cidflags.Add(cidflag);
                            hcdflags.Add(hcdflag);
                            etdflags.Add(etdflag);
                            datatype.Add(type);
                        }
                        break; 
                    default : 
                        break ; 
                }
            }

        }

   


        static void Main(string[] args)
        {
            if (args.Length == 1)
            {
               PrintUsage();
            }

            GlycoFragworkDLL.Classes.Params _params = new GlycoFragworkDLL.Classes.Params(); 

            string filename ; 
            string paramname ; 
            for (int argIndex = 0; argIndex < args.Length; argIndex++)
            {
                string commandLine = args[argIndex];

                switch (commandLine)
                {
                    case "-i":
                        filename = args[argIndex + 1];
                        ParseInputFile(filename);
                        break;
                    case "-p":
                        paramname = args[argIndex + 1];
                        _params.LoadParamsFromFile(paramname);
                        break;
                    default:
                        PrintUsage();
                        break;
                }
            }


          
            int last_pos = mlnfile.LastIndexOf(".");
            string outputfilecsv = null;
            string outputfilexml = null; 

            _params.PrintOutputAsCSV = true ;
            if (_params.RunInGlycoMode)
            {
                if (_params.PrintOutputAsCSV)
                    outputfilecsv = mlnfile.Substring(0, last_pos) +"_map_glycopeptides.csv";
                if (_params.PrintOutputAsXML)
                    outputfilexml = mlnfile.Substring(0, last_pos) + "_map_glycopeptides.xml";
            }
            else
            {
                if (_params.PrintOutputAsCSV)
                    outputfilecsv = mlnfile.Substring(0, last_pos) + "_map_nonglycopeptides.csv" ; 
                else
                    outputfilexml = mlnfile.Substring(0, last_pos) + "_map_nonglycopeptides.xml" ; 
            }


            int folder_pos = mlnfile.LastIndexOf("\\") ;
            string folder = mlnfile.Substring(0, folder_pos) + "\\CIDSequencing_output";

          
            //string deGlycoFile = @"D:\data\GlycoMapSera\MascotAnalysis\ForwardReverseSearch\BuildSummary\TestScore25EValue05.glycopeptides.txt" ; 

            // -- Initialize -- //
            GlycoFragworkDLL.Tasks.MapCreator tMapCreator = new GlycoFragworkDLL.Tasks.MapCreator();
            GlycoFragworkDLL.Tasks.MapParser tMapParser = new GlycoFragworkDLL.Tasks.MapParser();
            //GlycoFragworkDLL.Tasks.MapScorer tMapScorer = new GlycoFragworkDLL.Tasks.MapScorer();
            GlycoFragworkDLL.Tasks.MapWriter tMapWriter = new GlycoFragworkDLL.Tasks.MapWriter(); 
            GlycoFragworkDLL.Classes.MapRecord cMap = new GlycoFragworkDLL.Classes.MapRecord() ;
            //GlycoFragworkDLL.Tasks.MapMergeIDs tMapMergeIDs = new GlycoFragworkDLL.Tasks.MapMergeIDs(); 

 
            Console.WriteLine("Loading Parameters.."); 
            // -- Load in Data -- //
            tMapCreator.LoadParameters(_params);
           
            Console.WriteLine("Loading Map.."); 
            tMapCreator.LoadMap(ref cMap, mlnfile, datasetnames, cidflags, hcdflags, etdflags, datatype);
            
            /*if (!_params.RunInGlycoMode)
            {
                tMapMergeIDs.LoadIdentificationsFromFile(nonGpFile, cMap._AssociatedDatasetNames);  
            }
         //   tMapCreator.LoadDeglycoInfo(ref deGlycoMap, deGlycoMap); */
           
            
           
 
            // Load Map            
            //_params.ProcessOnlyNGlycopeptides = true;                         
            tMapCreator.LoadFragmentationSpectraIntoMap(ref cMap, ref _params, glycanFile, fastaFile, gpFile);
            
            tMapCreator.CleanUp();

            if (cidflags.Contains(true))
            {
                Console.WriteLine("Clustering records..");
                // Cluster CIDs together, break apart into seperate IDs
                tMapParser.ClusterRecordsOnCID(ref cMap, ref _params);
                tMapParser.GetRepresentatives(ref cMap, ref _params);
            }
            
           
            if(etdflags.Contains(true))
            {

                Console.WriteLine("Sequencing Clustered CID representatives..");                 
                tMapParser.SequenceCIDPeaks(ref cMap, ref _params, folder);
                Console.WriteLine("Assigning FDR.."); 
                tMapParser.AssignFDR(ref cMap, ref _params); 
            }

            if (_params.ProcessOnlyNGlycopeptides)
            {
                Console.WriteLine("Resolving..");
                tMapParser.ResolveIdentificationsAcrossClusters(ref cMap, ref _params);
            }
            else
            {
                Console.WriteLine("Choosing Representatives with lowest HCD..");
                tMapParser.GetRepresentativesAcrossClusters(ref cMap, ref _params);
            }

            
            // Following two are for inclusion list creation
           //string infofile = mlnfile.Substring(0, last_pos) + "_precursor_details.csv";
            //tMapWriter.WriteOutPrecursorInfoToFileV2(ref cMap, infofile); 

        /*    string inclusionlist = mlnfile.Substring(0, last_pos) + "_inclusion_list.csv" ; 
            tMapWriter.WriteOutInclusionListCandidates(ref cMap, inclusionlist) ; 
            

            // Score
            //tMapScorer.ScoreMap(ref cMap);            
            //tMapScorer.SetNglycopeptides(nglycopeptides); 



            // Scoring
         /*   tMapScorer.SetOptions(_params);
            tMapScorer.output_folder = folder; 
            tMapScorer.ScoreMapIndividual(ref cMap, ref _params);


            // Handle HCD and the sequence, glycan identifications
            tMapParser.CalculateHCDRepresentativeFragmentationSpectra(ref cMap, ref _params);

            /* Y1 inclusion List 
            string y1inclusionList = mlnfile.Substring(0, last_pos) + "_y1_inclusion_list.csv";
            tMapWriter.WriteOutY1InclusionListCandidates(ref cMap, y1inclusionList); 

            tMapParser.GetRepresentativeIdentificationInformation(ref cMap, ref _params);
            */


            

            // Write out Map to file 
            Console.WriteLine("Writing output..") ; 
            if (_params.PrintOutputAsCSV)
                tMapWriter.WriteOutMapToCSV(ref cMap, outputfilecsv);
            if (_params.PrintOutputAsXML)                
                tMapWriter.WriteOutMapToXML(ref cMap, outputfilexml, _params.ProcessOnlyNGlycopeptides); 

            Console.WriteLine("Done") ; 
        }


       
       

       
    }
}
