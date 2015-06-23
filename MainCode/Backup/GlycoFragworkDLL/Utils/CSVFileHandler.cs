using System;
using System.IO;
using System.Collections;   
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GlycoFragworkDLL.Utils
{
    public class CSVFileHandler
    {
        public static readonly int READ_ONLY = 0;
        public static readonly int WRITE_ONLY = 1;
        int _Mode;
        string _FileName;
        StreamReader _ReadStream = null;
        StreamWriter _WriteStream = null;


        //constructor
        public CSVFileHandler(String FileName, int Mode)
        {
            _FileName = FileName;
            _Mode = Mode;
        }


        public String getFileName
        {
            get
            {
                return _FileName;
            }
        }

        public int getMode
        {
            get
            {
                return _Mode;
            }
        }


        //open file
        public void openFile()
        {
           if (_Mode == CSVFileHandler.READ_ONLY)
                    _ReadStream = File.OpenText(_FileName);
                else if (_Mode == CSVFileHandler.WRITE_ONLY)
                    _WriteStream = File.CreateText(_FileName);
                else
                    throw new IOException("Unkown file mode.");
        }

        //close file
        public void closeFile()
        {
            if ((_Mode == CSVFileHandler.READ_ONLY) && (_ReadStream != null))
                _ReadStream.Close();
            else if (_WriteStream != null)
                _WriteStream.Close();
        }

        //read a line from the file and return the attributes from the line.
        public String[] readLine()
        {
            if (_ReadStream == null)
                throw new IOException("Reading file failed because the stream is null.");
            if (_ReadStream.Peek() != -1)
            {
                String Str = _ReadStream.ReadLine();
                return Str.Split(new char[] { ',' });
            }
            else
                return null;

        }

        //write a list of attributes to a line in the file.
        public void writeLine(String[] Attributes)
        {
            if (_WriteStream == null)
                throw new IOException("Writing file failed because the stream is empty.");
            for (int i = 0; i < Attributes.Length - 1; i++)
                _WriteStream.Write(Attributes[i] + ",");
            _WriteStream.Write(Attributes[Attributes.Length - 1] + "\n");
        }

        public void writeLine(String line)
        {
            if (_WriteStream == null)
                throw new IOException("Writing file failed because the stream is empty.");
            _WriteStream.Write(line + "\n");
        }
    }
}
