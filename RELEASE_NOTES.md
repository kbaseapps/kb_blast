### Version 1.3.5
__Changes__
- added HTML report to BLASTx

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
