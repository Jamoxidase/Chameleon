# Chameleon

UCSC TABI presents project Chameleon, a pipeline and software package that optimizes plasmids for bacterial transformation efficiency in non-model organisms.



## Description

Chameleon is a software project written for UCSC TABI as a submission for the Best Software award for iGEM 2023. For more information on anything here, visit our [team wiki](https://2023.igem.wiki/ucsc/software)

`pstealth` is a software pipeline included in Chameleon built around Stealth. Stealth is a bioinformatics tool written by our PI David L. Bernick here a UCSC. Stealth statistically finds and reports underrepresented kmer motifs within a genome. The original version of Stealth can be found [here](https://git.ucsc.edu/dbernick/stealth).

Chameleon is published as a python package to the Python Package Index (PyPI) under the package `chameleontools`. Find the PyPI page for the project [here](https://pypi.org/project/chameleontools/).


## Installation
Installation steps are tailored to Mac/Unix-based operating systems. For instruction for Windows, skip to **Windows OS**

##### MacOS/Linux (Ubuntu 20.04)

Chameleon can be installed under the Python package `chameleontools`. Installation is made simple with the pip python package installer. 

The pip installer should be automatically installed with most Python distributions. Check if you have properly installed python and pip installer, the output should look similar
```bash
usr:~$ python --version
Python 3.x.x
usr:~$ python -m pip --version
pip X.Y.Z from /<path>/<to>/<your>/pip (python 3.x.x)  
```
If your output looks similar, you can skip to **Installing Chameleon**.

##### Installing Python / pip
Follow this section if your output did not look the same as above.

If you do not have Python, get started by installing Python 3.10 or above from [python.org](https://www.python.org/downloads/), or through a distriubtion such as [anaconda](https://www.anaconda.com/download).

If you do not have a working pip installer, follow steps found [here](https://pip.pypa.io/en/stable/installation/).

Validate you have a working Python / pip installer by running 
```bash
usr:~$ python --version
Python 3.x.x
usr:~$ python -m pip --version
pip X.Y.Z from /<path>/<to>/<your>/pip (python 3.x.x)  
```
and making sure the output looks similar.

##### Installing Chameleon

Installing Chameleon is simple as running the command
```bash
usr:~$ pip install chameleontools
# OR
usr:~$ python -m pip install chameleontools
```
with a valid pip installer to install the chameleontools package from the Python Package Index (PyPI)

##### Windows Installation
Installation on Windows requires a```Path``` environment variable be set for ```pstealth``` to be run as a CLI tool. Firstly, make sure Python 3.10+ is installed from [python.org](https://www.python.org/downloads/) or other distribution.

Validate that Python and a valid pip installer are installed using the command in the command prompt or Windows powershell
```bat
C:\usr>python --version
Python 3.x.x
C:\usr>python -m pip --version 
pip X.Y.Z from C:\<path>\<to>\<your>\pip (python 3.x.x)
```

Install chameleontools with the command
```bat
C:\usr>python -m pip install chameleontools
```
If running `pstealth` gives an output that looks like the following
```bat
C:\usr>pstealth -h
'pstealth' is not recognized as an internal or external command,
operable program or batch file.
C:\usr>
```
continue to the following section. If pstealth runs normally (prints usage information), chameleontools is properly installed and you can ignore the following section.

###### Adding to Path 
You must add the Python `Scripts` directory to your `Path` environment variables. This can be done by navigating to your `Scripts` directory where `pstealth.exe` is located. This is most likely under the path `C:\Users\<USER>\AppData\Local\Programs\Python\Python3<x>\Scripts` where `<USER>` is the local Windows user profile and `Python3<x>` represents the version of your Python installation (eg. `Python312`). Once naviagted to your `Scripts` directory, save the filepath.

In the Windows search (hitting the windows button), type and enter `run` and then enter `sysdm.cpl` into the dialog box that pops up. This opens up a new window titled `System Properties`. From there navigate to the `Advance` tab and select `Environment Variables`. From there, scroll until you see a `Path` variable, select, and click `Edit`. In the `Edit Environment Variable` window, select `New` and paste the filepath you saved into the new text box.

If you do not see a `Path` variable, simply select `New` in the `Environment Variables` window and create a new variable named `Path` with the value of your saved file path.

Once done, `pstealth` should be able to be used from the commandline as intended.



## Usage

#### Plasmid-Stealth (pstealth) CLI tool
Once installed, the pipeline can be easily run with the command `pstealth`
```bash
# usage
pstealth --genome (-g) <genome infile> --plasmid (-p) <plasmid infile> --outfile -o [outfile | default: stdout] -[zPMmrs]
            Optional Args:
                --zScore (-z) -> [zscore cutoff value | default: -4]
                --pseudo (-P) [pseudo-count value | default: 0]
                --max (-M) [maximum motif size | default: 8]
                --min (-m) [minimum motif size | default: 1]
                --palindrome (-r) [Remove RC palindromes only | default: off]
                --silent (-s) [Hide report message | default: show]
```
The `pstealth` command takes two required arguments `--plasmid (-p)` and `--genome (-g)`. 
`--plasmid` is the annotated plasmid which you wish to Stealth optimize in GenBank format (.gb/.gbk)
`--genome` is the genomic sequence of the species you wish transform your plasmid into in GenBank or FastA format (.gb/.gbk || .fasta/.fa)

An example usage is as follows
```bash
usr:~$ pstealth -g M_aeruginosa.gb -p pSPDY.gb -o pSPDY_stealth.gb -r 
...
usr:~$ cat pSPDY_stealth.gb | head
LOCUS       pSPDY                   8825 bp    DNA     circular UNA 18-AUG-2023
DEFINITION  synthetic circular DNA.
ACCESSION
VERSION     .
KEYWORDS    .
SOURCE      synthetic DNA construct
  ORGANISM  synthetic DNA construct.
FEATURES             Location/Qualifiers
     misc_feature    1
```

The above example optimizes plasmid `pSPDY.gb` for transformation into species `M_aeruginosa.gb`, removing only reverse-complement palindrome motifs identified with StealthV0. The modifed plasmid is saved to a new file `pSPDY_stealth.gb`, keeping all the same annotations, only changing the plasmid sequence to remove underrepresented motifs in coding regions.

#### chameleontools Module

Development of chameleontools started when we wrote software to process the Stealth output of <em>M.aeruginosa</em> and realized we could write a pipeline to automate what was a manual and time consuming process. This software ultimately became what is now the `StealthParser` module.

We recognized the usefulness of having each submodule perform their own specific tasks and have written the chameleontools package in a way that reflects that ideal.

To import chameleontools as a Python module, simplly import it as follows
```python
import chameleontools.<submodule>
```
The submodules supported for importing include 
```py
import chameleontools.ChromatoSeq # Motif removal from sequence
import chameleontools.CodonAnalyzer # Codon Analysis of CDS
import chameleontools.FastAreader 
import chameleontools.ORFfinder # Very basic ORF finder
import chameleontools.SeqParser # Plasmid and Genome input handling
import chameleontools.Stealth # Stealth analysis
import chameleontools.StealthParser # Stealth output parser
```

Visit each submodule's folder for detailed README documentation. 


## Contributing

Chameleon was created for the 2023 iGEM competition and will not be maintained on this repository after the conclusion of this year's iGEM cycle. 

This software is published under the MIT license. Feel free to use any and all code provided by the project in any way and for any purpose.

## Authors and acknowledgment
Stealth was written by our PI David L. Bernick at UC Santa Cruz who advised the development of the `pstealth` pipeline (contact: dbernick@soe.ucsc.edu) 

Chameleon was written and contributed to by 
* Tyler Gaw (contact: tagaw@ucsc.edu)
* Allison Jaballas (contact: acjaball@ucsc.edu)
* James Larbaleister (contact: jlarbale@ucsc.edu)

Special thanks to Reto Stamm (contact: rstamm@ucsc.edu | [github](https://github.com/retospect)) for guidance in developing and publishing a packge to the Python Package Index