# Quantifying_Contribution_TranscriptionFactors_EColi
Code to quantitatively assess the contribution of transcription factors to biological processes in *Escherichia coli*

Description
This pipeline integrates coexpression analysis to distinguish between transcription factors (TFs) that drive the expression of a process and those that only cooperate with the leading TF.

If you have questions, problems or suggestions, contact Edna Karen Rivera Zagal (see contact info below)

## Features
-v 0.0.1
- Replicates the whole analysis

## Installation
Download pipeline from Github repository:

```bash 
git clone https://github.com/Ednakrz/Quantifying_Contribution_TranscriptionFactors_EColi.git
```
```bash
cd Quantifying_Contribution_TranscriptionFactors_EColi
``` 
## Test our pipeline

Estimated test time: 1 minute or less

- Go version
  ```bash
  bash Generator_analysis_ContriMatrix.sh GOs
  ```
- Pathway version
  ```bash
  bash Generator_analysis_ContriMatrix.sh Pathways
  ```
Pipeline results for test data should be in the following directory:
  ```bash
  ./results/
  ```

## Cite us
If you use our results or pipeline, please, cite us!

## Contact
If you have questions, requests, or bugs to report, open an issue in github, or email dledezma@ccg.unam.mx and ednakrz@lcg.unam.mx 

