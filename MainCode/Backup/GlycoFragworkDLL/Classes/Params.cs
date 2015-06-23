using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Xml; 

namespace GlycoFragworkDLL.Classes
{
    public class Params
    {
        GlypID.Scoring.clsScoringParameters _scoringParams;
        GlypID.HornTransform.clsHornTransformParameters _transformParams;
        GlypID.Peaks.clsPeakProcessorParameters _peakParams; 
        
        // UMC parameters
       double _maxUMCCoverage;
       bool _processOnlyNGlycopeptides;
       bool _useGlycopeptideFile;
       bool _createGlycopeptideFile; 
       bool _printAsCSV;
       bool _printAsXML; 
       bool _inGlycoMode;
       bool _use_decoy_glycan;
       bool _use_decoy_peptide;   
        


        public  Params()
        {
            _scoringParams = new GlypID.Scoring.clsScoringParameters();
            _transformParams = new GlypID.HornTransform.clsHornTransformParameters();
            _peakParams = new GlypID.Peaks.clsPeakProcessorParameters();
            _maxUMCCoverage = 0.2;
            _processOnlyNGlycopeptides = true;
            _useGlycopeptideFile = false;
            _printAsCSV = false;
            _inGlycoMode = true; 
        }

        public bool RunInGlycoMode
        {
            get { return _inGlycoMode; }
            set { _inGlycoMode = value; }
        }

        public GlypID.Scoring.clsScoringParameters ScoringParams
        {
            get { return _scoringParams; }
            set { _scoringParams = value; }
        }
        
        public GlypID.Peaks.clsPeakProcessorParameters PeakParams
        {
            get { return _peakParams; }
            set { _peakParams = value; }
        }

        public GlypID.HornTransform.clsHornTransformParameters TransformParams
        {
            get { return _transformParams; }
            set { _transformParams = value; }        
        }

        public double MaxUMCCoverage
        {
            get { return _maxUMCCoverage; }
            set { _maxUMCCoverage = value; }
        }
        public bool ProcessOnlyNGlycopeptides
        {
            get { return _processOnlyNGlycopeptides; }
            set { _processOnlyNGlycopeptides = value; }
        }

        public bool UseGlycoPeptideFile
        {
            get { return _useGlycopeptideFile; }
            set { _useGlycopeptideFile = value; }
        }
        
        public bool CreateGlycoPeptideFile
        {
            get { return _createGlycopeptideFile; }
            set { _createGlycopeptideFile = value; }
        }


        public bool PrintOutputAsCSV
        {
            get { return _printAsCSV; }
            set { _printAsCSV = value; }
        }
        public bool PrintOutputAsXML
        {
            get { return _printAsXML; }
            set { _printAsXML = value; }
        }
        public bool UseDecoyGlycan
        {
            get { return _use_decoy_glycan; }
            set { _use_decoy_glycan = value; }
        }
        public bool UseDecoyPeptide
        {
            get { return _use_decoy_peptide; }
            set { _use_decoy_peptide = value; }
        }

        /// <summary>
        /// Function to check the validity of parameters
        /// </summary>
        /// <returns></returns>
        public bool ParamsCheck()
        {
            bool safe = true;

            if (UseDecoyPeptide & UseDecoyGlycan)
            {
                safe = false;
                throw new Exception("Only one of peptide and glycan are allowed to be decoy");
            }

            return safe; 
        }

      
        private void LoadUMCParameters(XmlReader rdr)
        {
            while (rdr.Read())
            {
                switch (rdr.NodeType)
                {
                    case XmlNodeType.Element:
                        if (rdr.Name == "MaxUMCCoverage")
                        {
                            rdr.Read();
                            while (rdr.NodeType == XmlNodeType.Whitespace || rdr.NodeType == XmlNodeType.SignificantWhitespace)
                            {
                                rdr.Read();
                            }
                            if (rdr.NodeType != XmlNodeType.Text)
                                throw new Exception("Missing information for MaxUMCCoverage in parameter file");
                            this.MaxUMCCoverage = System.Convert.ToDouble(rdr.Value);
                        }
                        if (rdr.Name == "ProcessOnlyNGlycopeptides")
                        {
                             rdr.Read();
                            while (rdr.NodeType == XmlNodeType.Whitespace || rdr.NodeType == XmlNodeType.SignificantWhitespace)
                            {
                                rdr.Read();
                            }
                            if (rdr.NodeType != XmlNodeType.Text)
                                throw new Exception("Missing information for ProcessOnlyNglycopeptides in parameter file");
                            this.ProcessOnlyNGlycopeptides= System.Convert.ToBoolean(rdr.Value);
                        }
                        if (rdr.Name == "UseGlycopeptideFile")
                        {                            
                            rdr.Read();
                            while (rdr.NodeType == XmlNodeType.Whitespace || rdr.NodeType == XmlNodeType.SignificantWhitespace)
                            {
                                rdr.Read();
                            }
                            if (rdr.NodeType != XmlNodeType.Text)
                                throw new Exception("Missing information on whether to use glycopeptide file");
                            this.UseGlycoPeptideFile = Convert.ToBoolean(rdr.Value);
                            rdr.Read();
                        }
                        if (rdr.Name == "CreateGlycopeptideFile")
                        {
                            rdr.Read();
                            while (rdr.NodeType == XmlNodeType.Whitespace || rdr.NodeType == XmlNodeType.SignificantWhitespace)
                            {
                                rdr.Read();
                            }
                            if (rdr.NodeType != XmlNodeType.Text)
                                throw new Exception("Missing information on whether to use glycopeptide file");
                            this.CreateGlycoPeptideFile = Convert.ToBoolean(rdr.Value);
                            rdr.Read();
                        }
                        if (rdr.Name == "GlycoMode")
                        {
                            rdr.Read() ; 
                            while (rdr.NodeType == XmlNodeType.Whitespace || rdr.NodeType == XmlNodeType.SignificantWhitespace)
                            {
                                rdr.Read();
                            }
                            if (rdr.NodeType != XmlNodeType.Text)
                                throw new Exception("Missing information on whether to run in GlycoMode");
                            this.RunInGlycoMode = Convert.ToBoolean(rdr.Value);
                            rdr.Read();
                        }
                        if (rdr.Name == "OutputFormat")
                        {
                            rdr.Read();
                            while (rdr.NodeType == XmlNodeType.Whitespace || rdr.NodeType == XmlNodeType.SignificantWhitespace)
                            {
                                rdr.Read();
                            }
                            if (rdr.NodeType != XmlNodeType.Text)
                                throw new Exception("Missing information on output format");
                            if (rdr.Value == "CSV")
                            {
                                PrintOutputAsCSV = true;
                                PrintOutputAsXML = false;
                            }
                            else if (rdr.Value == "XML")
                            {
                                PrintOutputAsXML = true;
                                PrintOutputAsCSV = false;
                            }
                            else if (rdr.Value == "BOTH")
                            {
                                PrintOutputAsCSV = true;
                                PrintOutputAsXML = true;
                            }
                            else
                                throw new Exception("Invalid output format. Only XML, CSV or BOTH are possible"); 

                        }
                        if (rdr.Name == "UseDecoyGlycan")
                        {
                            rdr.Read();
                            while (rdr.NodeType == XmlNodeType.Whitespace || rdr.NodeType == XmlNodeType.SignificantWhitespace)
                            {
                                rdr.Read();
                            }
                            if (rdr.NodeType != XmlNodeType.Text)
                                throw new Exception("Missing information on using a decoy glycan approach");
                            this.UseDecoyGlycan = Convert.ToBoolean(rdr.Value);
                            rdr.Read();
                        }
                        if (rdr.Name == "UseDecoyPeptide")
                        {
                            rdr.Read();
                            while (rdr.NodeType == XmlNodeType.Whitespace || rdr.NodeType == XmlNodeType.SignificantWhitespace)
                            {
                                rdr.Read();
                            }
                            if (rdr.NodeType != XmlNodeType.Text)
                                throw new Exception("Missing information on using a decoy peptide approach");
                            this.UseDecoyPeptide = Convert.ToBoolean(rdr.Value);
                            rdr.Read();
                        }                       
                        break;
                    case XmlNodeType.EndElement:
                        if (rdr.Name == "UMCParameters")
                            return;
                        break; 
                    default:
                        break;
                }
            }
        }

        private void SaveUMCProcessingParameters(XmlTextWriter writer)
        {
            writer.WriteStartElement("UMCParameters");
            writer.WriteWhitespace("\n\t\t");
            writer.WriteElementString("MinUMCCoverage", this.MaxUMCCoverage.ToString());
            writer.WriteWhitespace("\n\t\t");
            writer.WriteElementString("ProcessOnlyNGlycopeptides", this.ProcessOnlyNGlycopeptides.ToString());
            writer.WriteWhitespace("\n\t\t");            
            writer.WriteEndElement(); 
            writer.WriteWhitespace("\n\t\t") ; 
        }

        public void LoadParamsFromFile(string param_file)
        {
            XmlTextReader rdr = new XmlTextReader(param_file);
            while (rdr.Read())
            {
                switch (rdr.NodeType)
                {
                    case XmlNodeType.Element:
                        if (rdr.Name == "PeakParameters")
                            _peakParams.LoadPeakParameters(rdr);
                        else if (rdr.Name == "HornTransformParameters")
                            _transformParams.LoadHornTransformParameters(rdr);
                        else if (rdr.Name == "ScoringParameters")
                            _scoringParams.LoadScoringParameters(rdr);
                        else if (rdr.Name == "Miscellaneous")
                            _transformParams.LoadMiscellaneousParameters(rdr);
                        else if (rdr.Name == "UMCParameters")
                            LoadUMCParameters(rdr);
                        break; 
                    default:
                        break; 
                }
            }
            rdr.Close();             
        }

        public void WriteParamsToFile(string param_file)
        {
            try
            {
                XmlTextWriter xwriter = new XmlTextWriter(param_file, System.Text.Encoding.UTF8);
                xwriter.Formatting = Formatting.None;
                xwriter.IndentChar = '\t';
                xwriter.Indentation = 1;
                xwriter.WriteStartDocument(true);


                xwriter.WriteWhitespace("\n");
                xwriter.WriteStartElement("parameters");
                xwriter.WriteWhitespace("\n\t");
                xwriter.WriteElementString("version", "1.0");
                xwriter.WriteWhitespace("\n\t");

                _peakParams.SavePeakParameters(xwriter);
                _transformParams.SaveHornTransformParameters(xwriter);
                _transformParams.SaveMiscellaneousParameters(xwriter);
                _scoringParams.SaveScoringParameters(xwriter);
                SaveUMCProcessingParameters(xwriter); 

                _transformParams.ElementIsotopeComposition.SaveElementIsotopes(xwriter);
                xwriter.WriteEndElement();
                xwriter.WriteEndDocument();
                xwriter.WriteWhitespace("\n");
                xwriter.Close(); 
            }
            catch (Exception ex)
            {
                Console.WriteLine(ex.Message); 
            }

        }

    }
}
