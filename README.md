# pathway_vari.py

## hangelog and bug information
Version 2.0         2020/7/15

## Purpose of the tool
This tool evaluates and visualizes which metabolic pathways the specified mutated genes accumulate in multigenerational mutants. In order to test whether or not the mutant genes have accumulated in each pathway, this tool performs statistical hypothesis testing by calculating the binomial probability. The null hypothesis is that mutated genes do not accumulate in the set of genes corresponding to the nodes in the pathway of interest. The alternative hypothesis is that mutated genes accumulate in the set of genes corresponding to the nodes in the pathway of interest. Instead of using p-value, we used Q-value for statistical hypothesis testing. 

## Requirements 
* Python3.7
* nsdm 


## Installation 
* nsdm

```
        pip install -U git+https://github.com/CompBio-TDU-Japan/nsdm
```


## Usage

pathway_vari.py works with the following command.

```
python payhway_vari.py [geneIDfile1]  [geneIDfile2] ・・・ [geneIDfileN][gfffile] 
```

[geneIDfile]   A list of gene IDs that are used in gfffile.

[gfffile]  The gff file of the target organism.


### option

#### -p , --ppi
Add `-p` or `--ppi` to the command to analyze protein-protein interaction data in an integrated metabolic pathway.
Enter the commands in the following order.

```
python payhway_vari.py [geneIDfile1] [geneIDfile2] ・・・[geneIDfileN] [gfffile] -p [ppifile]
```

The [ppifile] should use the full version of the target organism from the STRING database(url:https://string-db.org/cgi/download.pl?sessionId=l3jvlSNOTnAs) ).

With this option, you can specify the evidence type on the command.  Specify the name of evidence type or the following numbers after `-p` or `--ppi`.  If you use only `-p` or `--ppi`, information on all evidence types is given.

```

 evidence_type
                0.  neighborhood
                1.  neighborhood_transferred
                2.  fusion
                3.  cooccurence
                4.  homology
                5.  coexpression
                6.  coexpression_transferred
                7.  experiments
                8.  experiments_transferred
                9.  database
                10. database_transferred
                11. textmining
                12. textmining_transferred
```


The details of the evidence type follow STRING.(url: http://version10.string-db.org/help/getting_started/)

ex )

 ```
 python pathway_vari.py [geneIDfile1] [gfffile] -p [ppifile] 1  database
 ```
 
 When you run this command, information on protein-protein interactions that have an evidence type of  `neighborhood_transferred` and  `database` is added to the pathways.
 
 
## output
 
The following files and directories are output.

* q-value.html
* q-value.tsv
* pathway_images
* gene_result

### q-value.html
A graph showing the Q value of each pathway by generation.　where Q value is represented by storey method[1].  If we do this with test data, it looks like the following image.

<img src='https://github.com/CompBio-TDU-Japan/pathway_vari/blob/images/testgraph.png' width='520px'>


[1]John D. Storey .et al:Statistical significance for genomewide studies. Proc Natl Acad Sci U S A. 100: 9440–9445,2003

### q-value.tsv
A text summary of the Q values in q-value.html.

### pathway_images
This contains a colored diagram of mutation accumulation in the KEGGpathway(url:https://www.genome.jp/kegg/pathway.html) .


### On the terminal
You will see the following.

* pathway number
* pathway name
* number of variant proteins
* The number of proteins in each pathway
* p-value

