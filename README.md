# TATA-box-prediction
Predict TATA box with HMM data downloaded from EPD Databese(https://epd.epfl.ch/promoter_elements.php)in Python3.

@ EPD Database Citation: 
  1. The eukaryotic promoter database in its 30th year: focus on non-vertebrate organisms. Dreos, R., Ambrosini, G., Groux, R., Périer, R., Bucher, P. Nucleic Acids Res. (2017). doi: 10.1093/nar/gkw1069; PMID: 27899657.

  2. The Eukaryotic Promoter Database: expansion of EPDnew and new promoter analysis tools. Dreos, R., Ambrosini, G., Périer, R., Bucher, P. Nucleic Acids Res. (2014) doi: 10.1093/nar/gku1111; PMID: 25378343


Steps to do:
  1. Prepare a sequence file in fna(fasta) format and a genome annotation file in gff(gff3) format.
  2. Prepare a HMM matrix.(A default HMM TATA box related matrix downloaded from EPD has been uploaded in the repository.)
  3. Check the file names in the Python script.
  4. Run HMM.sh script.
  5. Label each predicted TATA box in a distance threshold with + or -.
  6. Develop statistic analysis with R plotting.(Distribution Plot, scatter plot, ROC plot)
  
Example results has also been uploaded.
