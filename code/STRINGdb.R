library(STRINGdb)
library(UniProt.ws)

Specie<-(UniProt.ws(taxId=9606))
Specie@taxIdUniprots<- as.vector(NoEC_OL.OX$V1)

# Retrieving the UNIPROT EC numbers
Data <- select(Specie, keys(Specie,"ENSEMBL"), "GENES", "ENSEMBL")


get_STRING_species(version="10", species_name="Homo sapiens")
noEC_OL <- STRINGdb$new()
Map = noEC_OL$map(Data, "ENSEMBL", removeUnmappedRows = TRUE )
Enrichment<-noEC_OL$get_enrichment(string_ids = Map$STRING_id)
Interactions<-noEC_OL$get_interactions(string_ids = Map$STRING_id)
Add<-noEC_OL$get_term_proteins(string_ids = Map$STRING_id,term_ids = Enrichment$term_id)
write.table(Interactions,file = "Interactions.tsv",sep = "\t")
