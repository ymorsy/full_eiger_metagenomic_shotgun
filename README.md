# gdc_kraken
Snake make workflow to set the project id in the config file (for example BRCA or COAD..): <br>
- download the genomic files. <br>
- convert to fastq. <br>
- run kraken2 against Bacteria database. <br>
- keep the results of kraken2 and only succesful aligned reads. <br> 
