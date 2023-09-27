# Skimming functionality

This has been deprecated.
Although it technically works (at least until first-order testing), it turns out to be almost unusable in practice (for large datasets read from DAS) because of:
- very slow remote input file reading
- very slow output file writing
- many random xrootd errors.

Instead of the code in this folder, use a NanoAOD-Tools based solution,
that can be submitted via CRAB:
https://github.com/GhentAnalysis/nanoSkimming
