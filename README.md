# CrisprDiag

## Background
The nucleic acid detection method based on the CRISPR/Cas system mainly consists of three steps, including nucleic acid amplification, specific nucleic acid sequence identification, and detection result reading. In the nucleic acid amplification stage, the small amount of DNA or RNA in the sample is amplified by polymerase chain reaction (PCR), recombinase polymerase amplification (RPA) or loop-mediated isothermal amplification (LAMP) for Cas /crRNA complex recognition.


## Retuirement
python >=3.5

## Usage
gene.fa & genome.fa  test file
```python
python RAVI.py -i gene.fa -g genome.fa -o output
```
```
Options:
    -h         Display help information
    -s         PAM sequence                                 <string>         [default:NGG]
    -r         pamori.enter 5 or 3                           <int>           [default:3]
    -l         Length of protospacer                         <int>           [default:20]
    -m         The num of mismatch                           <int>           [default:3]
    --amax     The maximum value of GC content              <int>           [default:80]
    --amin     The minimum value of GC content              <int>           [default:20]
    -1         PCR
    -2         RPA
    -3         LAMP
```
