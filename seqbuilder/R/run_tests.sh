

# Create a chain of SDs
./sample_bed_sd_files.R -t 10000 -s 500 -i 100 -d 500 -m 1 -o ../res/test.txt

# Convert to tsv
./sd_bed_format_and_qc.R ../res/test.txt ../res/test.tsv

# Run seqbuilder
Rscript seqbuilder_wrapper.R -l 10000 -s ../res/test.tsv -o test -c 500

# Test output
./sd_compare_input_output.R -p ../res/paf/test_chunked.paf -q ../res/test.tsv -n autotest -o ../resfile.txt -c 100
