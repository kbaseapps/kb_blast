FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# -----------------------------------------

# BLAST 2.13.0+ needs libgomp1
RUN apt-get update -y && \
  apt-get install -y libgomp1


ENV BLAST_VERSION='2.13.0'

WORKDIR /kb/module
RUN \
  curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz > ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz && \
  tar xfz ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz && \
  ln -s ncbi-blast-${BLAST_VERSION}+ blast && \
  rm -f ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz && \
  rm -f blast/bin/blastdb_aliastool && \
  rm -f blast/bin/blastdbcheck && \
  rm -f blast/bin/blastdbcmd && \
  rm -f blast/bin/blast_formatter && \
  rm -f blast/bin/convert2blastmask && \
  rm -f blast/bin/deltablast && \
  rm -f blast/bin/dustmasker && \
  rm -f blast/bin/legacy_blast.pl && \
  rm -f blast/bin/makembindex && \
  rm -f blast/bin/makeprofiledb && \
  rm -f blast/bin/rpsblast && \
  rm -f blast/bin/rpstblastn && \
  rm -f blast/bin/segmasker && \
  rm -f blast/bin/update_blastdb.pl && \
  rm -f blast/bin/windowmasker

COPY ./ /kb/module
RUN mkdir -p /kb/module/work && \
  chmod -R a+rw /kb/module && \
  make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
