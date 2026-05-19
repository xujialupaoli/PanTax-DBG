#themis
cd ${workpath}/themis
/usr/bin/time -v -o query_time.log themis profile  -r ./reads/read1.fq -r ./reads/read2.fq  \
--db-prefix  ./themis_db \
--ref-info RefDB_13404_genomes_info.txt \
--out query_themis --threads 64 -k 31 2>&1 |tee log_thmis

python evaluate_species_abundance.py \
  --real $ref_spyID_abundance \
  --predict predict_spy.ID.abundance \
  --output spy_abundance_evaluation.tsv


#ganon
cd ${workpath}/ganon
db="${workpath}/ganon/ganon_db"

/usr/bin/time -v -o query_time.log_ganonclassify  ganon classify --db-prefix $db --paired-reads $read1 $read2 --output-prefix results --report-type abundance -t $threads 2>&1 |tee log_ganonclassify
/usr/bin/time -v -o query_time.log_ganonreport ganon report -i results.rep --db-prefix $db --output-prefix tax_profile --report-type abundance -r all  2>&1 |tee log_ganonreport
python ganon_species_process.py tax_profile.tre
python /home/work/wenhai/metaprofiling/bacteria_refgenome_NCBIdata/alternative_methods_0208/evaluation_scripts/time_process.py query_time.log_ganonclassify > time_evaluation.txt


python evaluate_species_abundance.py \
  --real $ref_spyID_abundance \
  --predict predict_spy.ID.abundance \
  --output spy_abundance_evaluation.tsv



#centrifuger
cd ${workpath}/centrifuger
db=./centrifuger_db/centrifugerDB
/usr/bin/time -v -o query_time.log centrifuger -k 1 -x $db -1 $read1 -2 $read2 -t 64 > cls.tsv
python time_process.py query_time.log > time_evaluation.txt
"centrifuger-quant" -x $db -c cls.tsv > centrifuger_report.tsv


python evaluate_species_abundance.py \
        --real $ref_spyID_abundance \
        --predict predict_spy.abundance \
        --output spy_abundance_evaluation.tsv


#sylph
cd ${workpath}/sylph
db=./database.syldb
seq2tax=./seqid2taxid.map
/usr/bin/time -v -o query_time.log sylph profile $db -1 $read1 -2 $read2 -o result.tsv
python time_process.py query_time.log > time_evaluation.txt.query_time.log
python sylph_convert.py result.tsv $seq2tax species

#丰度
tail -n +2 sylph_abundance.txt|awk 'BEGIN{FS="\t";OFS="\t"}{print $0}' > predict_spy.ID.abundance
python evaluate_species_abundance.py \
  --real $ref_spyID_abundance \
  --predict predict_spy.ID.abundance \
  --output spy_abundance_evaluation.tsv







#kraken2
cd ${workpath}/Kraken2
Kraken2_DB="./kraken2_db"
/usr/bin/time -v -o query_time.log ${tools} --db $Kraken2_DB --output kraken2_query_reads --report kraken2_query_report --threads 64 --paired $read1 $read2
python time_process.py query_time.log > time_evaluation.txt

python evaluate_species_abundance.py \
        --real $ref_spyID_abundance \
        --predict predict_spy.abundance \
        --output spy_abundance_evaluation.tsv









#KMCP
cd ${workpath}/KMCP
scripts_dir=/home/work/wenhai/metaprofiling/bacteria_refgenome_NCBIdata/scripts
database_genomes_info=13404_strain_genomes_info.txt
db=./kmcp/kmcpDB/kmcp_refs_k21.kmcp
tax2genome=kraken2_strain_taxid.tsv
strain_taxonomy=taxonomy
/usr/bin/time -v -o query_time.log kmcp search --db-dir $db -1 $read1 -2 $read2 --out-file result.kmcp.gz --log result.kmcp.gz.log -j $threads
python time_process.py query_time.log > time_evaluation.txt
awk -F'\t' 'NR>1 {print $3,$1}' OFS='\t' $tax2genome > taxid.map
kmcp profile --taxid-map taxid.map --taxdump $strain_taxonomy result.kmcp.gz --out-file result.kmcp.profile --metaphlan-report result.metaphlan.profile --sample-id 0 --cami-report result.cami.profile --binning-result result.binning.gz --log result.kmcp.profile.log --level strain
python kmcp_strain_process.py result.kmcp.profile
python KMCP_add_taxid_abundance.py \
    --abundance strain_abundance.txt \
    --genome $database_genomes_info \
    --output strain_spy_taxid.txt
tail -n +2  strain_spy_taxid.txt |cut -f 4,5 |sort|uniq |tail -n +2 > predict_spy.ID.abundance
cut -f 1 predict_spy.ID.abundance >  predict_spy_number.ID
#评估
python evaluate_species_abundance.py \
  --real $ref_spyID_abundance \
  --predict predict_spy.ID.abundance \
  --output spy_abundance_evaluation.tsv




#bracken
cd ${workpath}/bracken

db="kraken2_db"
/usr/bin/time -v -o query_time.log bracken -d $db -i ../Kraken2/kraken2_query_report -o bracken_query_report -r 150 2>&1 |tee log_braken_run
python time_process.py query_time.log > time_evaluation.txt


python evaluate_species_abundance.py \
        --real $ref_spyID_abundance \
        --predict predict_spy.abundance \
        --output spy_abundance_evaluation.tsv


