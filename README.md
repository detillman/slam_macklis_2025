# Distinct regulatory elements control in vivo mRNA stability and subcellular localization in developing neurons.
Scripts for subtype-specific, in vivo SLAM-seq: quantification of T -> C conversions, elucidation of conserved sequence features, and identification of UTR-mediated axonal targeting.

All raw sequencing files, counts tables, and T -> C conversion quantifications are avabile on the Harvard Dataverse (https://doi.org/10.7910/DVN/QTKVVC) and NCBI GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309877).

Information for each script is provided below - refer to our preprint for additional information: https://doi.org/10.1101/2025.04.23.650263.

helper_script.R: Miscellaneous utilities and file locations.

annotate_mismatch_reads.py: Annoates mismatched reads for identification of T -> C conversions.

conversion_counter.py: Quantifies T -> C conversions.

defining_utrs.R: Delineates UTR sequence boundaries based on QuantSeq reads.

soma_conversions.R: Identifies T -> C conversions in CPN somata.

soma_labeling_validation.R: Validates T -> C conversions in CPN somata.

injection_strategies.ipynb: Compares labeling efficiency for subcutaneous and intraperitoneal injection strategies. 

soma_kinetics: Quantifies transcript half-lives in CPN somata.

soma_features.R: Elulcidates conserved sequence features across transcript classes in CPN somata.

are_scoring.R: Identifies AU-rich elements in CPN somata.

luciferase_reporterAssay.R: Determines regulatory roles of specific UTRs over transcript turnover and transport in CPN.

axon_conversions.R: Identifies T -> C coversions in CPN axons.

axon_analysis.R: Elucidates conserved sequence features across transcript classes in CPN axons.
