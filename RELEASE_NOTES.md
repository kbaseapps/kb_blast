### Version 1.7.0
__Changes__
- updated BLAST to 2.13.0
- fixed bug with long feature IDs
- added SpeciesTree as possible target
- add option to write off-genetic-code protein translations
- tightened up client instantiation

### Version 1.6.0
__Changes__
- updated BLAST to 2.12.0
- tidied up output FeatureSet obj names
- made all test data be uploaded

### Version 1.5.2
__Changes__
- tidied up App Docs

### Version 1.5.1
__Changes__
- allow multiple targets as input

### Version 1.5.0
__Changes__
- updated BLAST to 2.11.0

### Version 1.4.1
__Changes__
- Added genome obj name as output display option in addition to genome sci name

### Version 1.4.0
__Changes__
- Added support for AnnotatedMetagenomeAssembly as target and output hits as FeatureSet

### Version 1.3.5
__Changes__
- fixed 3*pro_len bug that was overfiltering BLASTx hits

### Version 1.3.4
__Changes__
- updated BLAST+ to 2.10.1

### Version 1.3.3
__Changes__
- set KBaseDataObjectToFileUtils to beta so can access AnnotatedMetagenomeAssemblyToFASTA() method for testing

### Version 1.3.1
__Changes__
- removed output_one_name arg.  set name automatically now

### Version 1.3.0
__Changes__
- updated BLAST to 2.10.0
- disabled psiBLAST until debugged
- disabled tBLASTn (only makes sense if searching scaffold directly)
- disabled tBLASTx (only makes sense if searching scaffold directly)
- added AnnotatedMetagenomeAssembly as valid target type
- added more unit tests

### Version 1.2.0
__Changes__
- updated BLAST to 2.9.0
- updated base Docker image to sdkbase2:latest
- moved repeated methods to Utils/BlastUtil.py

### Version 1.1.0
__Changes__
- patched CheckJob() "Bad Status Line" error
- updated BLAST to 2.8.1
- Update to Python 3

### Version 1.0.7
__Changes__
- changed citations to be PLOS format

### Version 1.0.6
__Changes__
- made protein sequence checking stricter

### Version 1.0.5
__Changes__
- bugfix for psiBLAST with MSA index
- changed base docker image to sdkbase2

### Version 1.0.4
__Changes__
- updated BLAST+ to v2.7.1
- changed kbase help contact

### Version 1.0.0
- Initial release (BLAST+ v2.6.0)
