# Chameleon

UCSC TABI presents Chameleon, a pipeline and collection of software based around Stealth that optimizes plasmids for bacterial transformation efficiency in non-model organisms.



## Description
Let people know what your project can do specifically. Provide context and add a link to any reference visitors might
be unfamiliar with (for example your team wiki). A list of Features or a Background subsection can also be added here.
If there are alternatives to your project, this is a good place to list differentiating factors.

## Installation
Installing chameleontools is made simple with the pip python package installer. 

The pip installer should be automatically installed with most python distributions. Check if you have properly installed python and pip installer, the output should look similar
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



## Contributing

Chameleon was created for the 2023 iGEM competition and will not be maintained on this repository after the conclusion of this year's iGEM cycle. The project may be maintained [here](https://google.com) in the future. Chameleon is licenced under the [license] so you are free to modify and distribute any code found. Feel free to use any and all code provided by the module to create a pipeline tailored to any specific use case.

## Authors and acknowledgment
Stealth was written by David L. Bernick, UC Santa Cruz who advised the development of the software pipeline built around Stealth
(contact: dbernick@soe.ucsc.edu) 

Chameleon and the `chameleontools` package was written by 
* Tyler Gaw (contact: tagaw@ucsc.edu)
* James Larbaleister (contact: jlarbale@ucsc.edu)

Special thanks to Reto Stamm (contact: rstamm@ucsc.edu | [github](https://github.com/retospect)) for guidance in developing and publishing a packge to the Python Package Index (PyPI)