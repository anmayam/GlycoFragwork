using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GlycoFragworkDLL.Classes
{
    public class UMCRecord
    {
        int _ID;
        int _Count ; 
        double _MW ; 
        double _MWStd ;
        int _DatasetID;  

      
        int _ScanStart;
        int _ScanEnd;
        double _Abundance;
        double _AvgFit ; 

        int _ScanRep ;
        int _CIDScanRep;
        int _HCDScanRep;
        int _ETDScanRep;

        double _CIDScanRepET;
        double _HCDScanRepET;         


        double _UMCRepMW ; 
        double _UMCRepMz;        
        short _UMCRepCS;
        double _UMCRepAbun ; 
        double _UMCElutionTime ;

        double _UMCRepHCDScore;
        double _UMCRepCIDScore;
        double _UMCRepETDScore;

        /*string _UMCRepProteinName;
        string _UMCRepPeptideSeq;
        string _UMCRepNGlycoSite;
        string _UMCRepPeptideMass;
        string _UMCRepGlycanComposition;
        double _UMCRepGlycanMass;
        double _UMCRepY1;
        double _UMCRepPPM;*/

        public List<FragEvents> _AssociatedFragEvents;
        public List<FragEvents> _ClusteredFragEvents; 


        public UMCRecord()
        {

            _ID = -1;
            _DatasetID = -1; 
            _Count = 0;
            _MW = 0;
            _MWStd = 0;             
            _ScanStart = 0;
            _ScanEnd = 0;
            _Abundance = 0;
            _AvgFit = 0;

            _ScanRep = 0;
            _CIDScanRep = 0;
            _HCDScanRep = 0;
            _ETDScanRep = 0;

            _CIDScanRepET = 0;
            _HCDScanRepET = 0;

            _UMCRepCS = 0;
            _UMCRepMz = 0;
            _UMCRepMW = 0;
            _UMCRepAbun = 0;
            _UMCElutionTime = 0; 

            

            _AssociatedFragEvents = new List<FragEvents>();
            _ClusteredFragEvents = new List<FragEvents>(); 
            
        }
        public int ID
        {
            get { return _ID; }
            set { _ID = value; }
        }
        public int DatasetID
        {
            get { return _DatasetID; }
            set { _DatasetID = value; }
        }
        public int Count
        {
            get { return _Count; }
            set { _Count = value; }
        }
        public double MW
        {
            get { return _MW; }
            set { _MW = value; }
        }
        public double MWStd
        {
            get { return _MWStd; }
            set { _MWStd = value; }
        }

        public int ScanStart
        {
            get { return _ScanStart; }
            set { _ScanStart= value; }
        }
        public int ScanEnd
        {
            get { return _ScanEnd; }
            set { _ScanEnd = value; }
        }       
        public double Abundance
        {
            get { return _Abundance; }
            set { _Abundance = value; }
        }
        public double AvgFit
        {
            get { return _AvgFit; }
            set { _AvgFit= value; }
        }
        public int ScanRep
        {
            get { return _ScanRep; }
            set { _ScanRep = value; }
        }  
        public int CIDScanRep
        {
            get { return _CIDScanRep; }
            set { _CIDScanRep = value; }
        }
        public double CIDScanRepET
        {
            get { return _CIDScanRepET; }
            set { _CIDScanRepET = value; }
        }
        public double HCDScanRepET
        {
            get { return _HCDScanRepET; }
            set { _HCDScanRepET = value; }
        }

        public int HCDScanRep
        {
            get { return _HCDScanRep; }
            set { _HCDScanRep = value; }
        }
        public int ETDScanRep
        {
            get { return _ETDScanRep; }
            set { _ETDScanRep = value; }
        }

        public double UMCRepMz
        {
            get { return _UMCRepMz; }
            set { _UMCRepMz = value; }
        }
        public short UMCRepCS
        {
            get { return _UMCRepCS; }
            set { _UMCRepCS = value; }
        }
        public double UMCRepMW
        {
            get { return _UMCRepMW; }
            set { _UMCRepMW = value; }
        }
        public double UMCRepAbun
        {
            get { return _UMCRepAbun; }
            set { _UMCRepAbun = value; }
        }
        public double UMCElutionTime
        {
            get { return _UMCElutionTime; }
            set { _UMCElutionTime = value; }
        }
        public double UMCRepHCDScore
        {
            get { return _UMCRepHCDScore; }
            set { _UMCRepHCDScore = value; }
        }
        public double UMCRepCIDScore
        {
            get { return _UMCRepCIDScore; }
            set { _UMCRepCIDScore = value; }
        }
        public double UMCRepETDScore
        {
            get { return _UMCRepETDScore; }
            set { _UMCRepETDScore = value; }
        }
            
        public void SetRepScores()
        {
            double min_hcd = 1 ; 
            double max_cid = 0 ; 
            double max_etd = 0 ;
            int hcd_scan = -1;
            int cid_scan = -1;
            int etd_scan = -1;

            for (int i = 0; i < _AssociatedFragEvents.Count; i++)
            {
                if (_AssociatedFragEvents[i].HCDScore < min_hcd)
                {
                    min_hcd = _AssociatedFragEvents[i].HCDScore;
                    hcd_scan = _AssociatedFragEvents[i].HCDScan;
                }
                if (_AssociatedFragEvents[i].CIDScore > max_cid)
                {
                    max_cid = _AssociatedFragEvents[i].CIDScore;
                    cid_scan = _AssociatedFragEvents[i].CIDScan;
                }
                if (_AssociatedFragEvents[i].ETDScore > max_etd)
                {
                    max_etd = _AssociatedFragEvents[i].ETDScore;
                    etd_scan = _AssociatedFragEvents[i].ETDScan;
                }
            }
            _UMCRepETDScore = max_etd;
            _UMCRepHCDScore = min_hcd;
            _UMCRepCIDScore = max_cid;
            _HCDScanRep = hcd_scan;
            _ETDScanRep = etd_scan;
            _CIDScanRep = cid_scan; 
        }

        public FragEvents FragEventWithLowestHCDScore( List<int> IdsToLook)
        {
            FragEvents t = new FragEvents();
            double min_hcd_score = 1;
            double prev_cid_score = 0 ; 
            int min_hcd_score_index = -1;

            for (int i = 0; i < _AssociatedFragEvents.Count; i++)
            {
                if ((_AssociatedFragEvents[i].HCDScore <= min_hcd_score) && (IdsToLook.Contains(_AssociatedFragEvents[i].ID)))
                {
                    if (_AssociatedFragEvents[i].HCDScore < min_hcd_score)
                    {
                        min_hcd_score_index = i;
                        min_hcd_score = _AssociatedFragEvents[i].HCDScore;
                        prev_cid_score = _AssociatedFragEvents[i].CIDScore; 
                    }
                    else
                    {
                        if (_AssociatedFragEvents[i].CIDScore > prev_cid_score)
                        {
                            min_hcd_score_index = i;
                            min_hcd_score = _AssociatedFragEvents[i].HCDScore;
                            prev_cid_score = _AssociatedFragEvents[i].CIDScore;
                        }
                    }
                }
            }
            if (min_hcd_score_index > -1)
                t = _AssociatedFragEvents[min_hcd_score_index];
            else
            {
                bool debug ;
                debug = true;
            }

            return t; 
        }
        public FragEvents FragEventWithHighestETDScore(bool use_associated, List<int> IdsToLook)
        {
            FragEvents t = new FragEvents();
            double max_etd_score = 0;
            int max_etd_score_index = -1;            
            if (use_associated)
            {
                for (int i = 0; i < _AssociatedFragEvents.Count; i++)
                {
                    if ((_AssociatedFragEvents[i].ETDScore > max_etd_score) && (IdsToLook.Contains(_AssociatedFragEvents[i].ID)))
                    {
                        max_etd_score_index = i;
                        max_etd_score = _AssociatedFragEvents[i].ETDScore;
                    }
                }
            }
            else
            {
                for (int i = 0; i < _ClusteredFragEvents.Count; i++)
                {
                    if ((_ClusteredFragEvents[i].ETDScore > max_etd_score) && (IdsToLook.Contains(_AssociatedFragEvents[i].ID)))
                    {
                        max_etd_score_index = i;
                        max_etd_score = _ClusteredFragEvents[i].ETDScore;
                    }
                }

            }
            if (max_etd_score > 0)
                t = _AssociatedFragEvents[max_etd_score_index];
            return t; 

        }
        public FragEvents FragEventWithHighestETDScore(bool use_associated)
        {
            FragEvents t = new FragEvents();
            double max_etd_score = 0; 
            int max_etd_score_index = -1;
            if (use_associated)
            {
                for (int i = 0; i < _AssociatedFragEvents.Count; i++)
                {
                    if (_AssociatedFragEvents[i].ETDScore > max_etd_score)
                    {
                        max_etd_score_index = i;
                        max_etd_score = _AssociatedFragEvents[i].ETDScore;
                    }
                }
            }
            else
            {
                for (int i = 0; i < _ClusteredFragEvents.Count; i++)
                {
                    if (_ClusteredFragEvents[i].ETDScore > max_etd_score)
                    {
                        max_etd_score_index = i;
                        max_etd_score = _ClusteredFragEvents[i].ETDScore;
                    }
                }

            }
            if (max_etd_score > 0)
                t = _AssociatedFragEvents[max_etd_score_index]; 
            return t; 

        }

        
       
        
    }
}
