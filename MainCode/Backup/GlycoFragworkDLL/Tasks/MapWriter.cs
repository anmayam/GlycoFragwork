using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using GlycoFragworkDLL.Utils;
using System.Xml;


namespace GlycoFragworkDLL.Tasks
{
    public class MapWriter
    {
        public void WriteOutY1InclusionListCandidates(ref Classes.MapRecord _glycoMap, string filename)
        {
           /* CSVFileHandler outfile = new CSVFileHandler(filename, CSVFileHandler.WRITE_ONLY);
            outfile.openFile();
            string[] header = { "ClusterID", "Peptide", "Protein", "PrecursorMZ", "PrecursorRT", "Y1Mz", "Y1RT", "Y1Intensity" };
            outfile.writeLine(header);

            GlypID.Peaks.clsPeakProcessor _HCDPeakProcessor = new GlypID.Peaks.clsPeakProcessor();






            _glycoMap._AllMLNRecords.ForEach(delegate(Classes.MultiAlignRecord m)
           {
               if (m.ID == 6)
               {
                   bool debug = true;
               }
               for (int i = 0; i < m._AssociatedUMCRecords.Count; i++)
               {
                   Classes.UMCRecord u = m._AssociatedUMCRecords[i];
                   if ((u.DatasetID == m._RepresentativeDatasetID_HCD) && (u.HCDPeaks != null)) // temp work around a bug
                   {

                       if (u.HCDPeaks.Length > 5)
                       {
                           GlypID.Peaks.clsPeak[] HCDPeaks = new GlypID.Peaks.clsPeak[1]; ;
                           float[] msms_mzs = new float[1];
                           float[] msms_intensities = new float[1];

                           msms_intensities = u.HCDIntensities;
                           msms_mzs = u.HCDMzs;

                           _HCDPeakProcessor.ProfileType = u.HCDProfileType;
                           double bk_intensity = GlypID.Utils.GetAverage(ref msms_intensities, ref msms_mzs, 700, 2000);
                           _HCDPeakProcessor.SetPeakIntensityThreshold(bk_intensity);
                           _HCDPeakProcessor.DiscoverPeaks(ref msms_mzs, ref msms_intensities, ref HCDPeaks, Convert.ToSingle(0), Convert.ToSingle(2000), false);
                           GlypID.HornTransform.clsHornTransformResults[] TransformResults = new GlypID.HornTransform.clsHornTransformResults[1];

                           TransformResults[0] = u.TransformResult;


                           // Get out all candidates                            
                           GlypID.ETDScoring.clsETDScoringScanResults[] ETDScoreResults = new GlypID.ETDScoring.clsETDScoringScanResults[m._PeptideSeq.Length];

                           List<string> peptides_as_string = new List<string>();

                           for (int k = 0; k < m._PeptideSeq.Length; k++)
                           {
                               ETDScoreResults[k] = new GlypID.ETDScoring.clsETDScoringScanResults();
                               ETDScoreResults[k].mstr_pep_seq = m._PeptideSeq[k];
                               ETDScoreResults[k].mstr_pro_seq_name = m._ProteinName[k];
                               ETDScoreResults[k].mdbl_seq_mass = m._PeptideMass[k];
                               ETDScoreResults[k].mdbl_glycan_mass = m._GlycanMass[k];
                               ETDScoreResults[k].mstr_glycan_composition = m._GlycanComposition[k];
                               ETDScoreResults[k].mstr_glyco_site = m._NGlycoSite[k];
                               peptides_as_string.Add(ETDScoreResults[k].mstr_pep_seq);
                           }

                           GlypID.HCDScoring.clsHCDScoring HCDScoring = new GlypID.HCDScoring.clsHCDScoring();
                           HCDScoring.SearchForY1Ion(ref HCDPeaks, ref ETDScoreResults, ref TransformResults);
                           for (int k = 0; k < ETDScoreResults.Length; k++)
                           {
                               double y1_mz = ETDScoreResults[k].mdbl_y1_mz;
                               if (y1_mz > 0)
                               {
                                   GlypID.Peaks.clsPeak y1Peak = new GlypID.Peaks.clsPeak();
                                   _HCDPeakProcessor.GetClosestPeakMz(y1Peak, Convert.ToSingle(y1_mz));

                                   if (y1Peak.mdbl_intensity > 0)
                                   {
                                       string[] line = { m.ID.ToString(), ETDScoreResults[k].mstr_pep_seq,ETDScoreResults[k].mstr_pro_seq_name ,
                                                        TransformResults[0].mdbl_mz.ToString(), u.ElutionTime.ToString(),
                                                        y1_mz.ToString(), u.HCDScanRepET.ToString(), y1Peak.mdbl_intensity.ToString() };
                                       outfile.writeLine(line);
                                   }
                                   else
                                   {
                                       bool debug = true;
                                   }
                               }
                           }

                       }
                   }


               }

           });

            outfile.closeFile();*/
        }


        public void WriteOutYxInclusionListCandidates(ref Classes.MapRecord _glycoMap, string filename)
        {
            /* Anoop November 2012 : Will have to recode this to use new sequencing information 
             CSVFileHandler outfile = new CSVFileHandler(filename, CSVFileHandler.WRITE_ONLY);
             outfile.openFile(); 
             string[] header = {"ClusterID","Peptide","Protein", "PrecursorMZ", "PrecursorRT", "YxMz", "YxRT"} ; 
             outfile.writeLine(header) ; 

             _glycoMap._AllMLNRecords.ForEach(delegate(Classes.MultiAlignRecord m)
             {
                 if (m.ID == 6)
                 {
                     bool debug = true;
                 }
                 for (int i = 0; i < m._AssociatedUMCRecords.Count; i++)
                 {
                     Classes.UMCRecord u = m._AssociatedUMCRecords[i];
                     if ((u.DatasetID == m._RepresentativeDatasetID_CID) && (u.CIDPeaks != null) ) // temp work around a bug
                     {
                         if (u.CIDPeaks.Length > 5)
                         {
                             // Get Spectra 
                             GlypID.Peaks.clsPeak[] CIDPeaks = u.CIDPeaks;
                             GlypID.HornTransform.clsHornTransformResults[] TransformResults = new GlypID.HornTransform.clsHornTransformResults[1];

                             TransformResults[0] = u.TransformResult;

                             float[] mzs = new float[CIDPeaks.Length];
                             float[] intensities = new float[CIDPeaks.Length];
                             int index = 0;
                             foreach (GlypID.Peaks.clsPeak peak in CIDPeaks)
                             {
                                 mzs[index] = Convert.ToSingle(CIDPeaks[index].mdbl_mz);
                                 intensities[index] = Convert.ToSingle(CIDPeaks[index].mdbl_intensity);
                                 index++;
                             }

                             // Get out all candidates
                             GlypID.ETDScoring.clsETDScoringScanResults[] ETDScoreResults = new GlypID.ETDScoring.clsETDScoringScanResults[m._PeptideSeq.Length];
                             List<string> peptides_as_string = new List<string>();
                             List<GlycoSeq.ProteinInfo> peptides_as_objects = new List<GlycoSeq.ProteinInfo>();

                             //List<GlycoSeq.ProteinInfo> ProteinList = GlycoSeq.ProteinParser(@"C:\development\GlycoFragwork\Test\Test.fasta");


                             for (int k = 0; k < m._PeptideSeq.Length; k++)
                             {
                                 GlycoSeq.ProteinInfo this_peptide = new GlycoSeq.ProteinInfo(m._ProteinName[k], m._PeptideSeq[k]);
                                 peptides_as_objects.Add(this_peptide);
                             }

                             // Create MS scan object
                             GlycoSeq.MSScan _msScan = new GlycoSeq.MSScan(CIDPeaks, mzs, intensities, Convert.ToSingle(TransformResults[0].mdbl_mz), Convert.ToInt32(TransformResults[0].mshort_cs), 1);


                             GlycoSeq.InclusionList IL = new GlycoSeq.InclusionList(_msScan, peptides_as_objects, 500.0f);

                             List<string> result = IL.Process(-1);

                             if (result.Count > 0)
                             {
                                 for (int zz = 0; zz < result.Count; zz++)
                                 {
                                     string[] yx_entry = result[zz].Split(',');
                                     string yx_peptide = yx_entry[2];
                                     double yx_mz = Convert.ToDouble(yx_entry[3]);
                                     string yx_protein = null;

                                     for (int k = 0; k < m._PeptideSeq.Length; k++)
                                     {
                                         if (m._PeptideSeq[k] == yx_peptide)
                                         {
                                             yx_protein = m._ProteinName[k];
                                             break;
                                         }
                                     }


                                     string[] line = { m.ID.ToString(), yx_peptide, yx_protein,
                                                     TransformResults[0].mdbl_mz.ToString(), u.ElutionTime.ToString(),
                                                     yx_mz.ToString(), u.CIDScanRepET.ToString() };
                                     outfile.writeLine(line);
                                 }
                             }
                             break;
                         }
                     }
                 }
             });
             outfile.closeFile(); 


             */


        }

        public void WriteOutPrecursorInfoToFileV2(ref Classes.MapRecord _glycoMap, string filename)
        {
                CSVFileHandler outfile = new CSVFileHandler(filename, CSVFileHandler.WRITE_ONLY);
            outfile.openFile();
            int numdatasets = _glycoMap._AssociatedDatasetNames.Count;

            string[] header = { "ClusterID", "Mass", "NET", "Protein", "Peptide", "Site", "Glycan", "PeptideMass", "GlycanMass", 
                              "TypeID", "TRUE_FALSE", "RepHCDScore", "RepCIDScore","RepETDScore", "RepCIDSeqScore","ParentMZ", "ParentScanTime"};
            for (int i = 0; i < numdatasets; i++)
            {
            }
            outfile.writeLine(header);

            _glycoMap._AllMLNRecords.ForEach(delegate(Classes.MultiAlignRecord m)
           {
               if (m.BestMatchProteinName != "")
               {
                   string true_false = "TRUE";
                   if (m.BestMatchFalseHit)
                       true_false = "FALSE";



                   string[] outline = { Convert.ToString(m.ID), Convert.ToString(m.Mass), Convert.ToString(m.NET), m.BestMatchProteinName, 
                                       m.BestMatchPeptideSeq, m.BestMatchNGlycoSite, m.BestMatchGlycanComposition, 
                                   m.BestMatchPeptideMass.ToString(), m.BestMatchGlycanMass.ToString(), 
                                   m.IDLabel.ToString(), true_false ,  m.RepHCDScore.ToString(), m.RepCIDScore.ToString(), m.RepETDScore.ToString(), m.RepCIDSequencingScore.ToString(), 
                                  m.BestMatchParentMz.ToString(), m.BestMatchParentScanTime.ToString()};
                   outfile.writeLine(outline);
               }
               
           });

           outfile.closeFile();


        }

        public void WriteOutPrecursorInfoToFile(ref Classes.MapRecord _glycoMap, string filename)
        {
            /*CSVFileHandler outfile = new CSVFileHandler(filename, CSVFileHandler.WRITE_ONLY);
            outfile.openFile();
            int numdatasets = _glycoMap._AssociatedDatasetNames.Count;

            string[] header = { "ClusterID", "Mass", "NET"};
            for (int i = 0; i < numdatasets; i++)
            {
                string[] ar2;
                ar2 = new string[header.Length + 5];
                header.CopyTo(ar2, 0);
                ar2.SetValue("DatasetID" + "." + i.ToString(), header.Length);
                ar2.SetValue("ClassRepMZ" + "." + i.ToString(), header.Length + 1);
                ar2.SetValue("ClassRepCS" + "." + i.ToString(), header.Length + 2);
                ar2.SetValue("Mass" + "." + i.ToString(), header.Length + 3);
                ar2.SetValue("ElutionTime" + "." + i.ToString(), header.Length + 4);
                header = ar2;
            }
            outfile.writeLine(header);

            _glycoMap._AllMLNRecords.ForEach(delegate(Classes.MultiAlignRecord m)
            {
                string[] outline = { Convert.ToString(m.ID), Convert.ToString(m.Mass), Convert.ToString(m.NET), m._RepresentativeTransformResult };
                m._AssociatedUMCRecords.Sort(delegate(Classes.UMCRecord u1, Classes.UMCRecord u2) { return u1.DatasetID.CompareTo(u2.DatasetID); });

                int count_index = 0;
                for (int i = 0; i < numdatasets; i++)
                {
                    string[] ar3 = new string[outline.Length + 5];
                    outline.CopyTo(ar3, 0);
                    Classes.UMCRecord u = m._AssociatedUMCRecords[count_index];
                    if (u.DatasetID == i)
                    {
                        ar3.SetValue(Convert.ToString(u.DatasetID), outline.Length);
                        ar3.SetValue(Convert.ToString(u.UMCRepMz), outline.Length + 1);
                        ar3.SetValue(Convert.ToString(u.UMCRepCS), outline.Length + 2);
                        ar3.SetValue(Convert.ToString(u.UMCRepMW), outline.Length + 3);
                        ar3.SetValue(Convert.ToString(u.UMCElutionTime), outline.Length + 4);
                        count_index++;
                        if (count_index >= m._AssociatedUMCRecords.Count)
                            count_index--; // To keep it within index.                        
                    }
                    else
                    {

                        ar3.SetValue(Convert.ToString(i), outline.Length);
                        ar3.SetValue(Convert.ToString(0), outline.Length + 1);
                        ar3.SetValue(Convert.ToString(0), outline.Length + 2);
                        ar3.SetValue(Convert.ToString(0), outline.Length + 3);
                        ar3.SetValue(Convert.ToString(0), outline.Length + 4);
                    }
                    outline = ar3;
                }
                outfile.writeLine(outline);
            });
            outfile.closeFile();*/

        }

        /*public void WriteOutMapToCsv(ref Classes.MapRecord _glycoMap, string filename)
        {
            CSVFileHandler outfile = new CSVFileHandler(filename, CSVFileHandler.WRITE_ONLY);
            outfile.openFile();
            int numdatasets = _glycoMap._AssociatedDatasetNames.Count;

            string[] header = { "ClusterID", "Mass", "NET" };
            for (int i = 0; i < numdatasets; i++)
            {
                string[] ar2;
                ar2 = new string[header.Length + 9];
                header.CopyTo(ar2, 0);
                ar2.SetValue("DatasetID" + "." + i.ToString(), header.Length);
                ar2.SetValue("HCDScore" + "." + i.ToString(), header.Length + 1);
                ar2.SetValue("CIDScore" + "." + i.ToString(), header.Length + 2);
                ar2.SetValue("ETDScore" + "." + i.ToString(), header.Length + 3);
                ar2.SetValue("Protein" + "." + i.ToString(), header.Length + 4);
                ar2.SetValue("Peptide" + "." + i.ToString(), header.Length + 5);
                ar2.SetValue("Site" + "." + i.ToString(), header.Length + 6);
                ar2.SetValue("Glycan" + "." + i.ToString(), header.Length + 7);
                ar2.SetValue("Abundance" + "." + i.ToString(), header.Length + 8);
                header = ar2;
            }
            outfile.writeLine(header);

            _glycoMap._AllMLNRecords.ForEach(delegate(Classes.MultiAlignRecord m)
            {
                string[] outline = { Convert.ToString(m.ID), Convert.ToString(m.Mass), Convert.ToString(m.NET) };

                m._AssociatedUMCRecords.Sort(delegate(Classes.UMCRecord u1, Classes.UMCRecord u2) { return u1.DatasetID.CompareTo(u2.DatasetID); });

                int count_index = 0;
                for (int i = 0; i < numdatasets; i++)
                {
                    string[] ar3 = new string[outline.Length + 9];
                    outline.CopyTo(ar3, 0);
                    Classes.UMCRecord u = m._AssociatedUMCRecords[count_index];
                    if (u.DatasetID == i)
                    {
                        ar3.SetValue(Convert.ToString(u.DatasetID), outline.Length);
                        ar3.SetValue(Convert.ToString(u.HCDScore), outline.Length + 1);
                        ar3.SetValue(Convert.ToString(u.CIDScore), outline.Length + 2);
                        ar3.SetValue(Convert.ToString(u.ETDScore), outline.Length + 3);
                        ar3.SetValue(u.ProteinName, outline.Length + 4);
                        ar3.SetValue(u.PeptideSeq, outline.Length + 5);
                        ar3.SetValue(u.NGlycoSite, outline.Length + 6);
                        ar3.SetValue(u.GlycanComposition, outline.Length + 7);
                        ar3.SetValue(Convert.ToString(u.Abundance), outline.Length + 8);
                        count_index++;
                        if (count_index >= m._AssociatedUMCRecords.Count)
                            count_index--; // To keep it within index.
                    }
                    else
                    {
                        ar3.SetValue(Convert.ToString(i), outline.Length);
                        ar3.SetValue(Convert.ToString(1), outline.Length + 1);
                        ar3.SetValue(Convert.ToString(0), outline.Length + 2);
                        ar3.SetValue(Convert.ToString(0), outline.Length + 3);
                        ar3.SetValue("", outline.Length + 4);
                        ar3.SetValue("", outline.Length + 5);
                        ar3.SetValue("", outline.Length + 6);
                        ar3.SetValue("", outline.Length + 7);
                        ar3.SetValue(Convert.ToString(0), outline.Length + 8);

                    }
                    outline = ar3;
                }
                outfile.writeLine(outline);




            });

            outfile.closeFile();
        }*/

        public void WriteOutMapToCSV(ref Classes.MapRecord _glycoMap, string filename)
        {
            CSVFileHandler outfile = new CSVFileHandler(filename, CSVFileHandler.WRITE_ONLY);
            outfile.openFile();
            int numdatasets = _glycoMap._AssociatedDatasetNames.Count;

            string[] header = { "ClusterID", "Mass", "NET", "Protein", "Peptide", "Site", "Glycan", "PeptideMass", "GlycanMass", 
                              "TypeID", "TRUE_FALSE", "RepHCDScore", "RepCIDScore","RepETDScore", "RepCIDSeqScore" };
            for (int i = 0; i < numdatasets; i++)
            {
                string[] ar2;
                ar2 = new string[header.Length + 5];
                header.CopyTo(ar2, 0);
                ar2.SetValue("DatasetID" + "." + i.ToString(), header.Length);
                ar2.SetValue("HCDScore" + "." + i.ToString(), header.Length + 1);
                ar2.SetValue("CIDScore" + "." + i.ToString(), header.Length + 2);
                ar2.SetValue("ETDScore" + "." + i.ToString(), header.Length + 3);              
                ar2.SetValue("Abundance" + "." + i.ToString(), header.Length + 4);
                header = ar2;
            }
            outfile.writeLine(header);

            _glycoMap._AllMLNRecords.ForEach(delegate(Classes.MultiAlignRecord m)
          {

              if (m.ID == 2290)
              {

                  bool test = true;
              }
              if (m.BestMatchProteinName != "")
              {
                  string true_false= "TRUE" ;
                  if(m.BestMatchFalseHit)
                      true_false = "FALSE"; 
                     
                  string[] outline = { Convert.ToString(m.ID), Convert.ToString(m.Mass), Convert.ToString(m.NET), m.BestMatchProteinName, 
                                       m.BestMatchPeptideSeq, m.BestMatchNGlycoSite, m.BestMatchGlycanComposition, 
                                   m.BestMatchPeptideMass.ToString(), m.BestMatchGlycanMass.ToString(), 
                                   m.IDLabel.ToString(), true_false ,  m.RepHCDScore.ToString(), m.RepCIDScore.ToString(), m.RepETDScore.ToString(), m.RepCIDSequencingScore.ToString()};

                  m._AssociatedUMCRecords.Sort(delegate(Classes.UMCRecord u1, Classes.UMCRecord u2) { return u1.DatasetID.CompareTo(u2.DatasetID); });

                  int count_index = 0;
                  for (int i = 0; i < numdatasets; i++)
                  {
                      string[] ar3 = new string[outline.Length + 5];
                      outline.CopyTo(ar3, 0);
                      Classes.UMCRecord u = m._AssociatedUMCRecords[count_index];
                      if (u.DatasetID == i)
                      {
                          ar3.SetValue(Convert.ToString(u.DatasetID), outline.Length);
                          ar3.SetValue(Convert.ToString(u.UMCRepHCDScore), outline.Length + 1);
                          ar3.SetValue(Convert.ToString(u.UMCRepCIDScore), outline.Length + 2);
                          ar3.SetValue(Convert.ToString(u.UMCRepETDScore), outline.Length + 3);
                          ar3.SetValue(Convert.ToString(u.Abundance), outline.Length + 4);
                          count_index++;
                          if (count_index >= m._AssociatedUMCRecords.Count)
                              count_index--; // To keep it within index.
                      }
                      else
                      {
                          ar3.SetValue(Convert.ToString(i), outline.Length);
                          ar3.SetValue(Convert.ToString(1), outline.Length + 1);
                          ar3.SetValue(Convert.ToString(0), outline.Length + 2);
                          ar3.SetValue(Convert.ToString(0), outline.Length + 3);
                          ar3.SetValue(Convert.ToString(0), outline.Length + 4);

                      }
                      outline = ar3;
                  }
                  outfile.writeLine(outline);
              }
          });

          outfile.closeFile();
        }

        public void WriteOutMapToXML(ref Classes.MapRecord _glycoMap, string filename, bool printall)
        {
            try
            {
                XmlTextWriter xwriter = new XmlTextWriter(filename, System.Text.Encoding.UTF8);
                xwriter.Formatting = Formatting.None;
                xwriter.IndentChar = '\t';
                xwriter.Indentation = 1;
                xwriter.WriteStartDocument(true);

                xwriter.WriteWhitespace("\n");
                xwriter.WriteStartElement("GlycoMap");
                xwriter.WriteWhitespace("\n\t");

                int numdatasets = _glycoMap._AssociatedDatasetNames.Count;
                List<string> names = _glycoMap._AssociatedDatasetNames;
                List<string> types = _glycoMap._AssociatedDatasetTypes;

                SpectralUtilities _spectralUtils = new SpectralUtilities();


                _glycoMap._AllMLNRecords.ForEach(delegate(Classes.MultiAlignRecord m)
                 {
                     if ((m.BestMatchProteinName != "") || (!printall))
                     {
                         string gsequence = " "; 
                         if ((m.BestMatchGlycanSequence !="") && (m.BestMatchGlycanSequence!= null))
                             gsequence = m.BestMatchGlycanSequence ; 

                         xwriter.WriteStartElement("GlycoRecord");
                         xwriter.WriteWhitespace("\n\t\t");

                         //Main Part
                         xwriter.WriteElementString("ID", m.ID.ToString());
                         xwriter.WriteWhitespace("\n\t\t");
                         xwriter.WriteElementString("Mass", m.Mass.ToString());
                         xwriter.WriteWhitespace("\n\t\t");
                         xwriter.WriteElementString("NET", m.NET.ToString());
                         xwriter.WriteWhitespace("\n\t\t");
                         xwriter.WriteElementString("Protein", m.BestMatchProteinName);
                         xwriter.WriteWhitespace("\n\t\t");
                         xwriter.WriteElementString("Peptide", m.BestMatchPeptideSeq);
                         xwriter.WriteWhitespace("\n\t\t");
                         xwriter.WriteElementString("PeptideMass", m.BestMatchPeptideMass.ToString());
                         xwriter.WriteWhitespace("\n\t\t");
                         xwriter.WriteElementString("Site", m.BestMatchNGlycoSite);
                         xwriter.WriteWhitespace("\n\t\t");
                         xwriter.WriteElementString("Glycan", m.BestMatchGlycanComposition);
                         xwriter.WriteWhitespace("\n\t\t");                         
                         xwriter.WriteElementString("GlycanSequence", gsequence);
                         xwriter.WriteWhitespace("\n\t\t");
                         xwriter.WriteElementString("GlycanMass", m.BestMatchGlycanMass.ToString());
                         xwriter.WriteWhitespace("\n\t\t");
                         string true_false = "TRUE";
                         if (m.BestMatchFalseHit)
                             true_false = "FALSE";
                         xwriter.WriteElementString("TRUE_FALSE", true_false);
                         xwriter.WriteWhitespace("\n\t\t");
                         xwriter.WriteElementString("IDType", m.IDLabel.ToString());
                         xwriter.WriteWhitespace("\n\t\t");

                         // The spectra
                         string repCID = " ";
                         int repCIDLength = 0;
                         if (m._RepresentativeCIDPeaks != null)
                         {
                             repCID = _spectralUtils.EncodePeaks(ref m._RepresentativeCIDPeaks);
                             repCIDLength = m._RepresentativeCIDPeaks.Length;
                         }
                         xwriter.WriteElementString("RepresentativeCIDLength", repCIDLength.ToString());
                         xwriter.WriteWhitespace("\n\t\t");
                         xwriter.WriteElementString("RepresentativeCIDSpectra", repCID);
                         xwriter.WriteWhitespace("\n\t\t");
                         string repHCD = " ";
                         int repHCDLength = 0;
                         if (m._RepresentativeHCDPeaks != null)
                         {
                             repHCD = _spectralUtils.EncodePeaks(ref m._RepresentativeHCDPeaks);
                             repHCDLength = m._RepresentativeHCDPeaks.Length;
                         }
                         xwriter.WriteElementString("RepresentativeHCDLength", repHCDLength.ToString());
                         xwriter.WriteWhitespace("\n\t\t");
                         xwriter.WriteElementString("RepresentativeHCDSpectra", repHCD);
                         xwriter.WriteWhitespace("\n\t\t");
                         string repETD = " ";
                         int repETDLength = 0;
                         if (m._RepresentativeETDPeaks != null)
                         {
                             repETD = _spectralUtils.EncodePeaks(ref m._RepresentativeETDPeaks);
                             repETDLength = m._RepresentativeETDPeaks.Length;
                         }
                         xwriter.WriteElementString("RepresentativeETDLength", repETDLength.ToString());
                         xwriter.WriteWhitespace("\n\t\t");
                         xwriter.WriteElementString("RepresentativeETDSpectra", repETD);
                         xwriter.WriteWhitespace("\n\t\t");


                         // The scores
                         xwriter.WriteElementString("RepCIDScore", m.RepCIDScore.ToString());
                         xwriter.WriteWhitespace("\n\t\t");
                         xwriter.WriteElementString("RepHCDScore", m.RepHCDScore.ToString());
                         xwriter.WriteWhitespace("\n\t\t");
                         xwriter.WriteElementString("RepETDScore", m.RepETDScore.ToString());
                         xwriter.WriteWhitespace("\n\t\t");
                         xwriter.WriteElementString("RepCIDSequencing", m.RepCIDSequencingScore.ToString());
                         xwriter.WriteWhitespace("\n\t\t");


                         // Print out all datasets but in order of dataset id
                         m._AssociatedUMCRecords.Sort(delegate(Classes.UMCRecord u1, Classes.UMCRecord u2) { return u1.DatasetID.CompareTo(u2.DatasetID); });
                         int count_index = 0;
                         xwriter.WriteStartElement("DatasetInfo");

                         for (int i = 0; i < numdatasets; i++)
                         {
                             xwriter.WriteWhitespace("\n\t\t\t");
                             xwriter.WriteStartElement("DatasetRecord");
                             Classes.UMCRecord u = m._AssociatedUMCRecords[count_index];

                             if (u.DatasetID == i)
                             {
                                 // There is a UMC present, so write it out.
                                 int last_pos = names[i].LastIndexOf("\\");
                                 string dataset_name = names[i].Substring(last_pos + 1); //, names[i].Length);

                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("DatasetID", u.DatasetID.ToString());
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("DatasetName", dataset_name);
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("DatasetType", types[i]);
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("Abundance", u.Abundance.ToString());
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("PrecursorMz", u.UMCRepMz.ToString());
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("PrecursorMass", u.UMCRepMW.ToString());
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("PrecursorCS", u.UMCRepCS.ToString());
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("PrecursorRT", u.UMCElutionTime.ToString());
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("HCDScan", u.HCDScanRep.ToString());
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("HCDScore", u.UMCRepHCDScore.ToString());
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("CIDScan", u.CIDScanRep.ToString());
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("CIDScore", u.UMCRepCIDScore.ToString());
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("ETDScan", u.ETDScanRep.ToString());
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("ETDScore", u.UMCRepETDScore.ToString());
                                 count_index++;
                                 if (count_index >= m._AssociatedUMCRecords.Count)
                                     count_index--; // To keep it within index.
                             }
                             else
                             {
                                 int last_pos = names[i].LastIndexOf("\\");
                                 string dataset_name = names[i].Substring(last_pos + 1); //, names[i].Length);

                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("DatasetID", i.ToString());
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("DatasetName", dataset_name);
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("DatasetType", types[i]);
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("Abundance", Convert.ToString(0));
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("PrecursorMz", Convert.ToString(0));
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("PrecursorMass", Convert.ToString(0));
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("PrecursorCS", Convert.ToString(0));
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("PrecursorRT", Convert.ToString(0));
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("HCDScan", Convert.ToString(0));
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("HCDScore", Convert.ToString(1));
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("CIDScan", Convert.ToString(0));
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("CIDScore", Convert.ToString(0));
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("ETDScan", Convert.ToString(0));
                                 xwriter.WriteWhitespace("\n\t\t\t\t");
                                 xwriter.WriteElementString("ETDScore", Convert.ToString(0));
                             }
                             xwriter.WriteWhitespace("\n\t\t\t");
                             xwriter.WriteEndElement(); // Dataset record

                         }
                         xwriter.WriteWhitespace("\n\t\t");
                         xwriter.WriteEndElement(); //Datasetinfo
                         xwriter.WriteWhitespace("\n\t");
                         xwriter.WriteEndElement(); // glycorecord
                         xwriter.WriteWhitespace("\n\t");
                     }


                 });

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