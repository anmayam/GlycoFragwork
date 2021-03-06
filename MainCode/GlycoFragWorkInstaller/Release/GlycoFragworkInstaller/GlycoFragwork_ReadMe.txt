GlycoFragwork version 0.1
---------------------------
Developed by Anoop Mayampurath, SoIC, Indiana University
Email Contacts:
Anoop Mayampurath: anmmayam@indiana.edu
Haixu Tang : hatang@indiana.edu

Main Notes
-----------
This package contains 3 folders

- GlycoFragwork_ReadMe.txt : this text

- Installer : containing the latest release of GlycoFragwork

- Sample : contains
	- Sample_input - sample of  input and parameter XML files for a Fetuin study. Also contains a sample input map acquired from MultiAlign analysis of three datasets.
	- Default_Combination_V2.csv : sample glycan composition file. This file was created based on synthetic N-glycosylation pathway. For details, see paper [Cite paer]
	- Sample_output - sample outputs corresponding to Fetuin study.  Both CSV and XML formats are reported.

- Extra Code
	- LDA_Analysis/CIDETDTrainMaps.R : R code for performing LDA on ETD and CID glycan sequencing in order to reclassify TRUE-FALSE entries based on a comined score. 


GlycoFragwork can identify glycopeptide ions within an aligned map of datasets, given a fasta file and a glycan list.

Installation
-------------
- Use the main installer file GlycoFragworkInstaller


Dependencies
-------------
- Windows XP/Vista with .NET framework 3.5.  Windows 7 users are encouraged to install XP mode on their systems (http://www.microsoft.com/windows/virtual-pc/download.aspx) and to install/run GlycoFragwork within that framework.
- Thermo Excalibur libraries must be installed, since the .RAW reader dlls are not distributable.  Support for mzXML has been incorporated. 
- The R code requires R version 2.10.1 or higher and the following packages : verification, ROCR, MASS, klaR.

Usage
------

GlycoFragwork is a command-line tool. Change dirrectory to where GlycoFragworkConsole.exe is.

C:\Installation_Directory>GlycoFragworkConsole.exe -i [input-file] -p [params=file]

Example Command Line

C:\Dir\GlycoFragworkConsole.ext -i input_sample_FetuinStudy.xml -p param_sample_FetuinStudy.xml


Input file
----------

The format of the input file is XML which can be edited in any text editor (e.g. gVim). Here is a detailed explanation of all the necessary tags

<MAPFile> - This is the output of MultiAlign v. 5.0.2 (http://omics.pnl.gov/software/MultiAlign.php : Console 32-bit).  Details of input are given in paper [CITE PAPER] (or see sample input file)

<GlycanFile> - Path of the glycan file to be used. An example file is given (Default_Combination_V2.csv) to be used as the default file.  Any other file can be used that keeps the same format and header information as Default_Combination_V2.csv.

<FastaFile> - Path of the fasta file to be used for glycopeptide identification

<GlycopeptideFile> - GlycoFragwork creates a glycopeptide file that combines tryptically digested peptides with sequons from the fasta file with each glycan in the glycan list. Thus, the fasta and the glycan file need only to be read in once. Reading the Glycopeptide file is much faster for small fastas ( < 100 glycoproteins). The glycopeptide file is read only if the UseGlycopeptideFile tag in parameters is set to TRUE. For larger fastas (>= 100 glycoproteins), there is no comparable speed up, so set UseGlycopeptideFile to FALSE.

<RAW> - Tag that represents the datasets used in the MultiAlign analysis
- <File> - Name of LC-MS/MS dataset (.RAW/ .mzXML)
- <CID> - Does dataset contain CID?
- <HCD> - Does dataset contain HCD?
- <ETD> - Does dataset contain ETD?


Parameter file
--------------

The format of the parameter file is XMl which can be edited in any text editor (e.g. gVim). Here is a detailed explanation of all necessary tags

<Peak Parameters> - These are peak parameters for peak picking of isotopic distributions in the MS scan. These are similar to the ones used in Decon2LS. For most cases, the default works very well.

<UMC Parameters>  - These are parameters for processing each LC-MS record cluster (or UMC)
- <MaxUMCCoverage> - An upper bound of elution time coverage for an LC-MS record to be considered. 0.5 would indicated GlycoFragwork to ignore all clusters that elute over 50% of the total dataset run (i.e. stop_scantime - start_scantime should < 50% of ttal LC time).
- <ProcessOnlyNglycopeptides> - Should be set to TRUE.
- <UseGlycopeptideFile> - If TRUE, the glycopeptide file mentioned in the input XML is read. If False, the fasta and glycan file are read in seperately.  Set to TRUE only for small databases (< 100 glycoproteins)
- <CreateGlycopeptideFile> - If True, a glycopeptide file is created which can be used in future runs to speed up processing.  Set to TRUE only for small databases ( < 100 glycoproteins). Set FALSE for large databases.
- <OutputFormat> - Format of output. Users can select - CSV, XML or BOTH.
- <UseDecoyPeptide> - Set to TRUE for FDR assignment.


<HornTransformParameters> - Same parameters as Decon2LS software tools.  

