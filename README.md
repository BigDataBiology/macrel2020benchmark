# MACREL Benchmarking (2019/20)

This repository includes code for benchmarking *MACREL*.

This is a companion repository to:

>   MACREL: antimicrobial peptide screening in genomes and metagenomes
>   Celio Dias Santos-Junior, Shaojun Pan, Xing-Ming Zhao, Luis Pedro Coelho
>   bioRxiv 2019.12.17.880385; doi:
>   [https://doi.org/10.1101/2019.12.17.880385](https://doi.org/10.1101/2019.12.17.880385)

## Contents

It contains the rules to rebuild the benchmarks in the paper.

However, instead just running the code, we strongly recommend you read it, as some steps depended on inputs obtained from manual curation

- To evaluate benchmarking results over tested AMP and hemolytic peptides prediction models, please refer to the *"train"* folder in [Macrel](https://github.com/BigDataBiology/macrel).

The other results showed in the MACREL benchmarking can be reproduced using the scripts in the following order:

(1) Benchmark.sh

(2) Macrel_in_real_metagenomes.sh

(3) Annotation_rules.sh

-- To generate Figure 3, please run:

```
$ python3 Figure_3_rendering.py
```

-- To generate Figure 4, please run:

```
$ ./python3 Figure_4_rendering.py
```

### Homology effect

In order to check homology in the training and testing data sets, please go to *"homology effects"* folder and follow the command:

```
$ ./retrain_complete.sh
```

This will retrain all models from MACREL, iAMP-2L and AMP Scanner v.2 with the non-redundant data sets, previously clustered with cd-hit at 80% of identity. The measures of accuracy, precision, and the confusion matrices will also be available. Be aware some of them can be generated in different time and will be printed in the screen.

## Third party softwares

In order to run all the codes, you will need besides MACREL:

- [Spurio](https://bitbucket.org/bateman-group/spurio/src/master/)
- [ArtMountRainier](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)
- BlastAll+
- [pigz](https://zlib.net/pigz/)
- R v3.5+
- [samtools](http://samtools.sourceforge.net/)
- [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)
- [Macrel](https://github.com/BigDataBiology/macrel)
- Python 3+
