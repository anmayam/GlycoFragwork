using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using GlycoFragworkDLL.Utils; 

namespace GlycoFragworkDLL.Classes
{
    public class MapRecord
    {
        public List<MultiAlignRecord> _AllMLNRecords;
        public List<string> _AssociatedDatasetNames;
        public List<bool> _IsCID;
        public List<bool> _IsHCD;
        public List<bool> _IsETD;
        public List<string> _AssociatedDatasetTypes ;

        public List<double> _AllFalseHitsFDRScore;
        public List<double> _AllTrueHitsFDRScore; 

        double _MAP_MIN_MASS = 10000;
        double _MAP_MAX_MASS = 0; 


        public MapRecord()
        {
            _AllMLNRecords = new List<MultiAlignRecord>();
            _AssociatedDatasetNames = new List<string>();
            _IsCID = new List<bool>();
            _IsHCD = new List<bool>();
            _IsETD = new List<bool>();
            _AssociatedDatasetTypes = new List<string>();
            _AllTrueHitsFDRScore = new List<double>();
            _AllFalseHitsFDRScore = new List<double>(); 
        }       

        public void AddRecord(MultiAlignRecord m_record)
        {
            _AllMLNRecords.Add(m_record); 
        }

        public void ClearRecords()
        {
            _AllMLNRecords.Clear(); 
        }

        public void SetDatasetNames(List<string> files)
        {
            _AssociatedDatasetNames = files; 
        }
        public void SetCID(List<bool> cid)
        {
            _IsCID = cid;
        }
        public void SetHCD(List<bool> hcd)
        {
            _IsHCD = hcd; 
        }
        public void SetETD(List<bool> etd)
        {
            _IsETD = etd; 
        }
        public void SetDatasetTypes(List<string> types)
        {
            _AssociatedDatasetTypes = types; 
        }
        public double MapMinMass
        {
            get { return _MAP_MIN_MASS; }
            set { _MAP_MIN_MASS = value; }
        }
        public double MapMaxMass
        {
            get { return _MAP_MAX_MASS; }
            set { _MAP_MAX_MASS = value; }
        }


        public void ReanInMLNV2File(string filename, int _numdatasets)
        {
            // Read in new version of MLN file
            CSVFileHandler CsvFileHandler2 = new CSVFileHandler(filename, CSVFileHandler.READ_ONLY);
            CsvFileHandler2.openFile();
            String[] Attributes2;
            CsvFileHandler2.readLine();
            while ((Attributes2 = CsvFileHandler2.readLine()) != null)
            {
                MultiAlignRecord mrecord = new MultiAlignRecord();
                mrecord.ID = int.Parse(Attributes2[0]);
                
                mrecord.Mass = double.Parse(Attributes2[1]);
                mrecord.NET = double.Parse(Attributes2[2]);
                int numDatasetsToAllocate = 0; 

                if (mrecord.Mass < _MAP_MIN_MASS)
                    _MAP_MIN_MASS = mrecord.Mass;
                if (mrecord.Mass > _MAP_MAX_MASS)
                    _MAP_MAX_MASS = mrecord.Mass;

                if ((Attributes2.Length-5) / 8 != _numdatasets)
                    numDatasetsToAllocate = (Attributes2.Length -5) / 8;
                else
                    numDatasetsToAllocate = _numdatasets; 

                mrecord.AllocateNumberDatasets(numDatasetsToAllocate);
                for (int i = 5; i < Attributes2.Length; i += 8)
                {
                    UMCRecord u = new UMCRecord();
                    if (Attributes2[i] != "")
                    {
                        u.ID = int.Parse(Attributes2[i]);
                        u.DatasetID = int.Parse(Attributes2[i + 1]);
                        u.MW = double.Parse(Attributes2[i + 2]);                        
                        u.Abundance = double.Parse(Attributes2[i + 4]);
                        u.ScanRep = int.Parse(Attributes2[i + 5]);
                        u.ScanStart = int.Parse(Attributes2[i + 6]);
                        u.ScanEnd = int.Parse(Attributes2[i + 7]);
                        mrecord._AssociatedUMCRecords.Add(u);
                    }
                }               
                AddRecord(mrecord);
                if (mrecord._AssociatedUMCRecords.Count > 10)
                {
                    bool debig = true;
                }
            }
        }

        public void ReadInMLNFile(string filename, int _numdatasets)
        {
            // Read in MultiAlign file
            // Note : sep 26, 2011 - This will need to be updated when the new MLN format comes out
            CSVFileHandler CsvFileHandler2 = new CSVFileHandler(filename, CSVFileHandler.READ_ONLY);
            CsvFileHandler2.openFile();
            String[] Attributes2;
            CsvFileHandler2.readLine();
            while ((Attributes2 = CsvFileHandler2.readLine()) != null)
            {
                MultiAlignRecord mrecord = new MultiAlignRecord();
                mrecord.ID = int.Parse(Attributes2[0]);
                mrecord.Size = int.Parse(Attributes2[1]);
                mrecord.Mass = double.Parse(Attributes2[2]);
                mrecord.NET = double.Parse(Attributes2[3]);
                mrecord.AllocateNumberDatasets(_numdatasets);
                int n_dataset = 0;
                for (int i = 4; i < Attributes2.Length; i += 4)
                {
                    UMCRecord u = new UMCRecord();
                    if (Attributes2[i] != "")
                    {
                        u.ScanRep = int.Parse(Attributes2[i]);
                        u.Abundance = double.Parse(Attributes2[i + 1]);
                        u.ScanStart = int.Parse(Attributes2[i + 2]);
                        u.ScanEnd = int.Parse(Attributes2[i + 3]);
                    }
                    mrecord._AssociatedUMCRecords.Add(u);
                    n_dataset++;
                }
                AddRecord(mrecord);
            }
        }

       

    }
}
