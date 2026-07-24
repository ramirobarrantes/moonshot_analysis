from /Users/rbarrant/Desktop/CurrentWork/moonshot_pipeline_analysis/testData

rsync -avz rbarrant@login.vacc.uvm.edu:/netfiles/davidkragseqdata/moonshot_run_output/MS008/sarek_output/annotation/ .
rsync -avz --exclude bam/ rbarrant@login.vacc.uvm.edu:/netfiles/davidkragseqdata/moonshot_run_output/MS008/moonshotAnalysis/ .
rsync -avz --exclude bam/ rbarrant@login.vacc.uvm.edu:/netfiles/davidkragseqdata/moonshot_run_output/MS008/seqinspector/multiqc .

and from /Users/rbarrant/Desktop/CurrentWork/moonshot_pipeline_analysis



docker run --rm \
  -e ANTHROPIC_API_KEY="$ANTHROPIC_API_KEY" \
  -v "$PWD/testData":/data/testData \
  moonshot_pipeline ...

docker run --rm \
  -e ANTHROPIC_API_KEY="$ANTHROPIC_API_KEY" \
  -v /Users/rbarrant/Desktop/CurrentWork/moonshot_pipeline_analysis/testData:/data/testData \
  -v /Users/rbarrant/Desktop/CurrentWork/moonshot_pipeline_analysis/DBs:/data/DBs \
  moonshot-analysis:latest \
  /data/testData/samplesheet.csv \
  /data/DBs/ensembl_peptides/ensembl_peptides.csv \
  /data/DBs/uniprot_database/HUMAN_9606_Database.RDS \
  /data/DBs/Cosmic_MutantCensus_Tsv_v104_GRCh37/Cosmic_MutantCensus_v104_GRCh37.tsv \
  MS008 \
  /data/testData/output

docker run --rm \
  -v /Users/rbarrant/Desktop/CurrentWork/moonshot_pipeline_analysis/testData:/data/testData \
  -v /Users/rbarrant/Desktop/CurrentWork/moonshot_pipeline_analysis/DBs:/data/DBs \
  moonshot-analysis:latest \
  /data/testData/samplesheet.csv \
  /data/DBs/ensembl_peptides/ensembl_peptides.csv \
  /data/DBs/uniprot_database/HUMAN_9606_Database.RDS \
  /data/DBs/Cosmic_MutantCensus_Tsv_v104_GRCh37/Cosmic_MutantCensus_v104_GRCh37.tsv \
  MS008 \
  /data/testData/output
