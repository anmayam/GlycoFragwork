/* Code to create an XML format input
 * 
 * bool use_only_n_glycopeps = true;
          List<string> datasetnames = new List<string>();
          List<bool> cidflags = new List<bool>();
          List<bool> hcdflags = new List<bool>();
          List<bool> etdflags = new List<bool>();

          string mlnfile = @"D:\data\GlycoFragwork\Combined\Map\combinedmnl2_map.csv"; 
          string rawfile1 = @"D:\data\GlycoFragwork\Combined\Raw\combined_101111_111012161157.raw" ;
          string fastaFile = @"D:\data\Fasta\CombinedMixture.fasta";
          datasetnames.Add(rawfile1);
          cidflags.Add(true);
          hcdflags.Add(true);
          etdflags.Add(true);
 
          string ttest = @"C:\development\GlycoFragwork\Test\testCombined.xml" ;
          XmlTextWriter twriter = new XmlTextWriter(ttest, System.Text.Encoding.UTF8);
          twriter.Formatting = Formatting.None;
          twriter.IndentChar = '\t';
          twriter.Indentation = 1;
          twriter.WriteStartDocument(true);
          twriter.WriteWhitespace("\n");
          twriter.WriteStartElement("parameters");
          twriter.WriteWhitespace("\n\t");
          twriter.WriteElementString("version", "1.0");
          twriter.WriteWhitespace("\n\t");
          twriter.WriteStartElement("Input");
          twriter.WriteWhitespace("\n\t\t");
          twriter.WriteElementString("MAPFile", mlnfile);
          twriter.WriteWhitespace("\n\t\t");
          twriter.WriteElementString("GlycanFile", glycanFile);
          twriter.WriteWhitespace("\n\t\t");
          twriter.WriteElementString("FastaFile", fastaFile);
          twriter.WriteWhitespace("\n\t\t");
          twriter.WriteEndElement();
          twriter.WriteWhitespace("\n\t\t");
            
          twriter.WriteStartElement("RAW") ;
          twriter.WriteElementString("File", rawfile1);
          twriter.WriteWhitespace("\n\t\t");
          twriter.WriteEndElement();
          twriter.WriteWhitespace("\n\t\t");
            
          twriter.WriteStartElement("CIDFlags");
          twriter.WriteElementString("CID", true.ToString());
          twriter.WriteWhitespace("\n\t\t");
          twriter.WriteEndElement();
          twriter.WriteWhitespace("\n\t\t");

          twriter.WriteStartElement("HCDFlags");
          twriter.WriteElementString("HCD", true.ToString());
          twriter.WriteWhitespace("\n\t\t");
          twriter.WriteEndElement();
          twriter.WriteWhitespace("\n\t\t");

          twriter.WriteStartElement("ETDFlags");
          twriter.WriteElementString("ETD", true.ToString());
          twriter.WriteWhitespace("\n\t\t");
          twriter.WriteEndElement();
          twriter.WriteWhitespace("\n\t\t");
            
          twriter.WriteEndElement();
          twriter.WriteEndDocument();
          twriter.WriteWhitespace("\n");
          twriter.Close(); */


/** // Temp code to create inclusion lists (to be inserted in the main program
            /*string conflictedUMCs = @"D:\data\GlycoMapSera\AllDatasests\Map\Map1\InclusionListBuilding\QueryConflictedUMCs.csv";
            GlycoFragworkDLL.Utils.CSVFileHandler conflictCSV = new GlycoFragworkDLL.Utils.CSVFileHandler(conflictedUMCs, GlycoFragworkDLL.Utils.CSVFileHandler.READ_ONLY);
            List<int> conflictedUMCIDs = new List<int>(); 
            conflictCSV.openFile();
            String[] Attributes2;
            conflictCSV.readLine();
            while ((Attributes2 = conflictCSV.readLine()) != null)
            {
                int umc = int.Parse(Attributes2[0]);
                conflictedUMCIDs.Add(umc); 
            }
            conflictCSV.closeFile();
            tMapCreator.LoadPrecursorInfoIntoMap(ref cMap, conflictedUMCIDs);
            tMapWriter.WriteOutPrecursorInfoToFile(ref cMap, outputfile); */
