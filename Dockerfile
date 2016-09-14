FROM kbase/kbase:sdkbase.latest
MAINTAINER KBase Developer
# -----------------------------------------

# Insert apt-get instructions here to install
# any required dependencies for your module.



## Update Transform (should go away eventually)
#RUN \
#  . /kb/dev_container/user-env.sh && \
#  cd /kb/dev_container/modules && \
#  rm -rf transform && \ 
#  git clone https://github.com/kbase/transform -b develop
#
## setup the transform, but ignore errors because sample data cannot be found!
#RUN \
#  . /kb/dev_container/user-env.sh; \
#  cd /kb/dev_container/modules/transform/t/demo; \
#  python setup.py; \
#  exit 0;


# Make sure SSL certs are properly installed  (NOT WORKING FOR ME ON APPDEV)
#RUN apt-get install python-dev libffi-dev libssl-dev \
#    && pip install pyopenssl ndg-httpsclient pyasn1 \
#    && pip install requests --upgrade \
#    && pip install 'requests[security]' --upgrade

# Install KBase Data API Library + dependencies
RUN mkdir -p /kb/module && cd /kb/module && git clone https://github.com/kbase/data_api && \
    mkdir -p lib/ && cp -a data_api/lib/doekbase lib/ && \
    pip install -r /kb/module/data_api/requirements.txt


# Install SDK utils
#RUN mkdir -p /kb/module && cd /kb/module && git clone https://github.com/kbaseapps/KBaseDataObjectToFileUtils && \
#    mkdir -p lib/ && cp -a KBaseDataObjectToFileUtils/lib/* lib/


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
RUN \
  curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.4.0+-x64-linux.tar.gz > ncbi-blast-2.4.0+-x64-linux.tar.gz && \
  tar xfz ncbi-blast-2.4.0+-x64-linux.tar.gz && \
  ln -s ncbi-blast-2.4.0+ blast && \
  rm -f ncbi-blast-2.4.0+-x64-linux.tar.gz && \
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


WORKDIR /kb/module
ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
