# TATA-box-prediction
Predict TATA box with HMM data downloaded from EPD in Python3

Steps to do:
  1. Prepare a sequence file in fna(fasta) format and a genome annotation file in gff(gff3) format.
  2. Prepare a HMM matrix.(A default HMM TATA box related matrix downloaded from EPD has been uploaded in the repository.)
  3. Check the file names in the Python script.
  4. Run HMM.sh script.
  5. Label each predicted TATA box in a distance threshold with + or -.
  6. Develop statistic analysis with R plotting.(Distribution Plot, scatter plot, ROC plot)
  
Example results has been uploaded.
