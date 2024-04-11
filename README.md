# Project Plan

I plan to build a genetic analysis toolkit that processes raw sequencing data to identify genetic missense variants and annotate them with potential impact on health or disease. 

Step 1: Call Raw Sequencing Data  
-Input raw sequencing data of potential pathogenic missense variants  
-Perform necessary alignment of the raw sequencing data to the reference sequencing data at hand

Step 2: Call Variant Annotation Module  
-Use online databases (such as ClinVar) to annotate variants that are pathogenic  
-Use pathogenicity scores predicted by AI models (such as AlphaMissense and ESM1b) to annotate variants that are pathogenic  
-Output flagged missense variants that could be potentially pathogenic 

Step 3: Population Genetics Analysis  
-Calculate relevant population genetics statistics (such as minor allele frequency)  
-Generate relevant plots for visual aid
