A suite of test cases for mpileup variant calling

### Adding a new test
Create a bam and fasta slice
```
./misc/create-bam-test -f ref.fa -b file.bam -r 12:34567 -o new-puzzzling-case
```
and copy the slices to ./dat. Then either create a new file or add to an existing file in ./tests the command and expected result, for example
```
# A brief description what is this test about, ideally with a snapshot if it should help understanding
# IGV snapshot: new-puzzling-case.png
cmd: bcftools mpileup -t 1:234 -f new-puzzling-case.fa new-puzzling-case.bam
fmt: %CHROM:%POS %REF %ALT
exp: chr1:129 CCCTTGCTGTG C
```

### Running the tests
```
# run all tests
./run-tests.pl

# run only tests that have been curated (have comments) and print only the failed ones
./run-tests.pl -c -q

# run only tests that use new-puzzling-case.bam
./run-tests.pl -r new-puzzling-case.bam
```
