# 1 - Added --pipeline argument to the parser to ingest skim2mito or skim2ribo, so that neither need to be the final base directory.
# 2 - Replaced the basename inference with a direct assignment from the argument so pipeline outputs can be located elsewhere than jsut the working directory.
# 3 - Added --trans_table argument so that the translation table needed for EMBLmyGFF3 can be changed.
# 4 - Fixed average_coverage error with contig-level mitogenome assemblies, which either resulted in 'NameError: name 'average_coverage' is not defined' or average_coverage getting populated from the previous sample. The COVERAGE line in the manifest is now only written for mitochondrial complete and mitochondrial partial assemblies.
