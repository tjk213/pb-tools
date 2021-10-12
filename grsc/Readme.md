# Guide-RNA Sample Compiler

Script for compiling GRNA sample files & filtering by RNA sequence.

Run with:
```
./grsc.py > out.csv
```

Sample File Format:

  - No quotes
  - Header line: `gene,sgRNA,gRNASeq,Count` (Note: no leading comma)
  - Data lines: `<id>,<gene>,<gene>_<seq>,<seq>,<count>`

See `S2wk_B2.csv` for an example.
