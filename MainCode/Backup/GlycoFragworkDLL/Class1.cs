using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Xml; 

namespace GlycoFragworkDLL
{
    public class Class1
    {

        void TempCode()
        {
            string glycanFile = @"C:\development\GlycoFragwork\GlycoFragworkDLL\DLLs\Default_Combination.csv";
                      

            string mlnfile = @"D:\data\GlycoFragwork\Combined\Map\combinedmnl2_map.csv";
            string rawfile1 = @"D:\data\GlycoFragwork\Combined\Raw\combined_101111_111012161157.raw";
            string fastaFile = @"D:\data\Fasta\CombinedMixture.fasta";
          

            string ttest = @"C:\development\GlycoFragwork\Test\testCombined.xml";
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

            twriter.WriteStartElement("RAW");
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
            twriter.Close(); 


        }
    }
}
