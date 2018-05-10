#!/bin/bash
for i in `seq 119 200`;
do
   celldb-etl-rnaseqer 10.50.101.122 http://www.ebi.ac.uk/fg/rnaseq/api homo_sapiens --limit 1 --offset $i
done    
        
