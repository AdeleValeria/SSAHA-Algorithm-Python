The coding and testing of the program were done using a notebook with the following 
specifications: AMD Ryzen 7 4800H processor (8 cores), 8 GB RAM and 512 GB SSD. The 
software used are Ubuntu 18.04 LTS and Spyder with built-in Python 3.8.5. Imported packages 
include tabulate v0.8.7, requests v2.24.0 and optparse v1.5.3. The program also interacts with 
Ensembl’s REST API. I tested my program on two reference genome files: 
Arabidopsis_thaliana.TAIR10.dna.toplevel.fa(ftp://ftp.ensemblgenomes.org/pub/plants/rele
ase-49/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz) and 
Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa(ftp://ftp.ensemblgenomes.org/pub/me
tazoa/release-49/fasta/drosophila_melanogaster/dna/). I used AT1G04645 sequence 
(https://www.arabidopsis.org/servlets/TairObject?type=sequence&id=1002493798) for A. 
thaliana and simulated sequence for D. melanogaster (generated by a Python function 
provided in source code).
