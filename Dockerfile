FROM kbase/kbase:sdkbase.latest
MAINTAINER KBase Developer
# -----------------------------------------

# Insert apt-get instructions here to install
# any required dependencies for your module.



# Update Transform (should go away eventually)
RUN \
  . /kb/dev_container/user-env.sh && \
  cd /kb/dev_container/modules && \
  rm -rf transform && \ 
  git clone https://github.com/kbase/transform -b develop

# setup the transform, but ignore errors because sample data cannot be found!
RUN \
  . /kb/dev_container/user-env.sh; \
  cd /kb/dev_container/modules/transform/t/demo; \
  python setup.py; \
  exit 0;



# RUN apt-get update

# -----------------------------------------

# Install SDK Module
#
RUN mkdir -p /kb/module
COPY ./ /kb/module
RUN mkdir -p /kb/module/work
WORKDIR /kb/module
RUN make


# Install BLAST+
#
WORKDIR /kb/module
RUN curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.3.0+-x64-linux.tar.gz > ncbi-blast-2.3.0+-x64-linux.tar.gz
RUN tar xfz ncbi-blast-2.3.0+-x64-linux.tar.gz
RUN ln -s ncbi-blast-2.3.0+ blast
RUN rm -f ncbi-blast-2.3.0+-x64-linux.tar.gz
RUN rm -f blast/bin/blastdb_aliastool
RUN rm -f blast/bin/blastdbcheck
RUN rm -f blast/bin/blastdbcmd
RUN rm -f blast/bin/blast_formatter
RUN rm -f blast/bin/convert2blastmask
RUN rm -f blast/bin/deltablast
RUN rm -f blast/bin/dustmasker
RUN rm -f blast/bin/legacy_blast.pl
RUN rm -f blast/bin/makembindex
RUN rm -f blast/bin/makeprofiledb
RUN rm -f blast/bin/rpsblast
RUN rm -f blast/bin/rpstblastn
RUN rm -f blast/bin/segmasker
RUN rm -f blast/bin/update_blastdb.pl
RUN rm -f blast/bin/windowmasker


WORKDIR /kb/module
ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
