using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GlycoFragworkDLL.Classes
{
    public class NonGlycopeptideRecord
    {
        GlypID.Sequence.clsSequence _sequence;
        int _datasetID;
        int _datasetScan;
        double _mass;
        short _cs; 
        

        public NonGlycopeptideRecord()
        {
            _sequence = new GlypID.Sequence.clsSequence(); 
        }

        public GlypID.Sequence.clsSequence Sequence
        {
            get { return _sequence; }
            set { _sequence = value; }
        }

        public int DatasetID
        {
            get { return _datasetID; }
            set { _datasetID = value; }
        }

        public int DatasetScan
        {
            get { return _datasetScan; }
            set { _datasetScan = value; }
        }

        public double Mass
        {
            get { return _mass; }
            set { _mass = value; }
        }

        public short CS
        {
            get { return _cs; }
            set { _cs = value; }
        }
        public bool unassigned = true; 
 
       
        
    }
}
