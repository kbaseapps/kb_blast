After you generate the GTDB species representative genome objects for
GTDB-Tk and load them into public workspaces, do the following for kb_blast
to access their proteomes:


Make BLAST dbs (for BLASTp App in kb_blast module)
--------------------------------------------------

  A. get mapping from Genome ID to UPA and save as Genome_UPAs-GTDB.tsv
  
     can use the one from kb_gtdbtk: kb_gtdbtk/data/Genome_UPAs-GTDB_r207_r214.tsv

  B. run scripts/add_UPA_to_genome_faa_header.pl for each group
     (Archaea-RS, Archaea-GB, Bacteria-RS, Bacteria-GB)

    $ ./scripts/add_UPA_to_genome_faa_header.pl -genomeidsfile lists/reps/IDs-Archaea-RS.list -idtoupafile tables/Genome_UPAs-GTDB_r207.tsv -faaindir gtdb_downloads_extracted/reps/proteome/archaea -stagingdir proteomes/reps/Archaea-RS

  C. merge upa-decorated faa files into one per group
    $ cat proteomes/reps/Archaea-RS/*faa > blast_dbs/reps/Archaea-RS.faa
#    $ gzip proteomes/reps/Archaea-RS.faa
    $ rm -rf proteomes/reps/Archaea-RS

  D. make blast formatted dbs (use kb_blast to run bash shell with data mounted)
    copy test_local/run_bash.sh to test_local/run_bash_with_mnt.sh and add the argument "-v /home/ac.dchivian/proj/genome_ref_dbs/GTDB/GTDB_r214.1/blast_dbs/reps:/src_data"

    RUN
    $ ./test_local/run_bash_with_mnt.sh
    % cd /src_data
    % makeblastdb -in Archaea-RS.faa -parse_seqids -dbtype prot -out Archaea-RS
    % makeblastdb -in Archaea-GB.faa -parse_seqids -dbtype prot -out Archaea-GB
    % makeblastdb -in Bacteria-RS.faa -parse_seqids -dbtype prot -out Bacteria-RS
    % makeblastdb -in Bacteria-GB.faa -parse_seqids -dbtype prot -out Bacteria-GB
    % chown 21610 ./*
    % chgrp 20001 ./*
    % exit

    foreach Archaea-RS, Archaea-GB, Bacteria-RS, Bacteria-GB
    $ tar cfz GTDB_r214-sprep-Archaea-RS-blast_2.14.1_dbs.tgz Archaea-RS.*p??
    $ rm Archaea-RS.*p??
    $ gzip Archaea-RS.faa
    $ mv Archaea-RS.faa.gz proteomes/reps

  E. copy blast db tarballs to public download site at NERSC (make sure dirs are world r-x)
    $ scp blast_dbs/reps/GTDB_r214-sprep-*-blast_2.14.1_dbs.tgz dylan@perlmutter.nersc.gov:/global/cfs/cdirs/kbase/www/GTDB/r214.1/blast_dbs/v1.1.0/

    log into perlmutter and
    $ chmod 644 /global/cfs/cdirs/kbase/www/GTDB/r214.1/blast_dbs/v1.1.0/*
