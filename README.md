# TATA box prediction
Predict TATA box with HMM data downloaded from EPD Databese (https://epd.epfl.ch/promoter_elements.php) in Python3.

@ EPD Database Citation: 
  1. _The eukaryotic promoter database in its 30th year: focus on non-vertebrate organisms. Dreos, R., Ambrosini, G., Groux, R., Périer, R., Bucher, P. Nucleic Acids Res. (2017)._

  2. _The Eukaryotic Promoter Database: expansion of EPDnew and new promoter analysis tools. Dreos, R., Ambrosini, G., Périer, R., Bucher, P. Nucleic Acids Res. (2014)._


Steps to do:
  1. Prepare a sequence file in .fna(.fasta) format and a genome annotation file in .gff(.gff3) format.
  2. Prepare a HMM matrix. (A default HMM TATA box related matrix downloaded from EPD has been uploaded in the repository.)
  3. Check the file names in the Python script.
  4. Run HMM.sh script.
  5. Label each predicted TATA box in a distance threshold with + or -.
  6. Develop statistic analysis with R plotting. (Distribution Plot, scatter plot, ROC plot)
  
Example results has also been uploaded.
