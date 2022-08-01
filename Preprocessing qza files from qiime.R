##### PREPROCESAMIENTO DE LOS FILES .QZA DE QIIME2 ####

########## IMPORTAR LOS OTUs IDENTIFICADOS ##########

#for import
read_qza("table-deblur.qza")
otus_table = read_qza("table-deblur.qza")
otus_table$data
head(otus_table$data,3) 
#VER LAS PRIMERAS 3 LINEAS DE CADA MUESTRA
#write.table(otus_table$data,"otus_table_abundance.tsv")
read.table("otus_table_abundance.tsv", header = T)
table_abundance = read.table("otus_table_abundance.tsv", header = T)
table_abundance$OTUID = row.names(table_abundance) 
#Para que con los nombres de los renglones hagas una columna adicional con el nombre "OTUID"
head(table_abundance, 3)
#Esta columna OTUID nos va a servir como punto de union entre 2 archivos
taxonomy = read_qza("taxonomy_Shigella.qza")
taxonomy$data
taxonomy = taxonomy$data
dim(taxonomy)
colnames(taxonomy)=c("OTUID", "Taxonomy", "Confidence")
#c De concatenar
head(taxonomy, 3)
taz=cSplit(taxonomy, "Taxonomy", ";")
head(taz, 3)
colnames(taz)= (c("OTUID", "Confidence", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Specie"))
head(taz, 3)
OTUS=merge(table_abundance, taz, by="OTUID")
#Mezlco ambas tablas por una columna que tenian en comun
head(OTUS, 3)
#write.table(OTUS, "merge_otu_tax.tsv")
#Crear la base de datos con la que vamos a trabajar

########## COMIENZA A TRABAJAR CON phyloseq ##########

# otu and taxonomy table
merge = read.table("merge_otu_tax.tsv", header = T)
head(merge, 3)
menus = merge
dim(menus)
row.names(menus)=menus$OTUID
#Agregamos la columna OTUID como nombre de filas
taxo = menus[,43:49]
head(taxo, 3)
head(menus, 3)
abundance=menus[,2:41 ]

#metadata file
metadata=read.table("Documento Metadata.txt", header=T, sep="\t", row.names = 1)
head(metadata, 3)

##### phyloseq file
taxo=as.matrix(taxo)
abundance=as.matrix(abundance)
abundance=otu_table(abundance, taxa_are_rows = T)
taxo=tax_table(taxo)
uni = phyloseq(abundance, taxo, sample_data(metadata))
uni

## filtrated data only bacteria D_0
bacteria = subset_taxa(uni, Kingdom==c("D_0__Bacteria"))
bacteria
