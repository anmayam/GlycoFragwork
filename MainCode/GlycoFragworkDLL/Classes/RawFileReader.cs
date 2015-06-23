using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ThermoRawFileReaderDLL.FinniganFileIO;


namespace GlycoFragworkDLL.Classes
{
    class RawFileReader
    {
        XRawFileIO _RawData ;
        FinniganFileReaderBaseClass.udtScanHeaderInfoType _header = new FinniganFileReaderBaseClass.udtScanHeaderInfoType();


        void ReadHeaderInfo(int scan)
        {
            bool success = _RawData.GetScanInfo(scan,  out _header); 
        }

        public RawFileReader()
        {
            _RawData = new XRawFileIO();
        }

        public RawFileReader(string inpFile)
        {
            _RawData = new XRawFileIO();
            _RawData.OpenRawFile(inpFile);
        }

        public void GetRawData(int scan, ref float[] mzs, ref float[] intensities)
        {
            double[] dmzs = new double[1];
            double[] dintensities = new double[1];
            int numpoints = _RawData.GetScanData(scan, ref dmzs, ref dintensities, ref _header) ;
            mzs = new float[numpoints];
            intensities = new float[numpoints];
            for (int i = 0; i < numpoints; i++)
            {
                mzs[i] = Convert.ToSingle(dmzs[i]);
                intensities[i] = Convert.ToSingle(dintensities[i]); 
            }
        }


        public bool IsMSScan(int scan)
        {
            ReadHeaderInfo(scan);         
            if(_header.MSLevel >1)
                return false ; 
            else
                return true ; 
        }

        public bool IsCentroidScan(int scan)
        {
            ReadHeaderInfo(scan);
            return _header.IsCentroidScan;

        }

        public int GetMSLevel(int scan)
        {
            ReadHeaderInfo(scan);
            return _header.MSLevel; 
        }


        public int GetParentScan(int scan)
        {
            int parent_scan = -1;
            int msn_level = GetMSLevel(scan); 
            bool found = false ;
            scan--;
            int temp_level = GetMSLevel(scan); 
            while (temp_level==msn_level && scan>1)
            {
                scan--;
                temp_level = GetMSLevel(scan); 
            }
            parent_scan = scan;
            return parent_scan; 
        }

        public double GetParentMz(int scan)
        {
            ReadHeaderInfo(scan); 
            return _header.ParentIonMZ;
        }

        public int GetNumScans()
        {
            return _RawData.GetNumScans(); 
        }

        public double GetParentMonoMzFromHeader(int scan)
        {
            ReadHeaderInfo(scan);
            return Convert.ToDouble(_header.ScanEventValues[10]);
        }

        public double GetParentChargeFromHeader(int scan)
        {
            ReadHeaderInfo(scan);
            return Convert.ToDouble(_header.ScanEventValues[9]);
        }

        public double GetScanTime(int scan)
        {
            // To do
            return 0;
        }

        public bool IsProfileScan(int scan)
        {
            // To Do
            return true;
        }

        public bool IsCIDScan(int scan)
        {
            // To do
            return false;
        }

        public bool IsHCDScan(int scan)
        {
            // to do
            return false;
        }

        public bool IsETDScan(int scan)
        {
            // to do
            return  false; 
        }

        public bool IsFTScan(int scan)
        {
            // to do
            return true; 
        }

    }
}
