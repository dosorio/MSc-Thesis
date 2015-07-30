#!/bin/bash
# USE: bash Proteomics.sh mgf_folder outputfolder
# Daniel Camilo Osorio
# Maestría en Bioinformática - Universidad Nacional de Colombia (Bogotá)
# Laboratorio de Bioquímica Teórica y Bioinformática - Pontificia Universidad Javeriana (Bogotá)

# Configuración de Idioma
export LC_ALL=C

# Genera la base de datos decoy a partir de la base de datos descargada de UniProt.
java -cp ~/SearchGUI/SearchGUI-1.26.4.jar eu.isas.searchgui.cmd.FastaCLI -in DB/HSA/HSA-210715.fasta -decoy

# Genera el archivo de parametros.
java -cp ~/SearchGUI/SearchGUI-1.26.4.jar eu.isas.searchgui.cmd.IdentificationParametersCLI -db ~/DB/HSA/HSA-210715_concatenated_target_decoy.fasta -out p.parameters -frag_tol 0.02 -prec_tol 6 -prec_ppm 1 -enzyme "Trypsin" -fixed_mods "carbamidomethyl c" -variable_mods "oxidation of m" -min_charge 2 -max_charge 4 -mc 2 -fi b -ri y

# Genera la carpeta resultados.
name=`basename $2`
mkdir ~/$2

# Crea enlaces a los archivos mgf de fragmentos m/z
ln -s ~/$1/*.mgf ~/$2/

# Genera la comparación de peptidos usando SearchGUI
java -Xmx8000m -cp ~/SearchGUI/SearchGUI-1.26.4.jar eu.isas.searchgui.cmd.SearchCLI -spectrum_files ~/$2/ -output_folder ~/$2/ -id_params p.parameters -xtandem 0 -msgf 0 -omssa 0 -ms_amanda 1 -myrimatch 0 -comet 0 -tide 0 -species "Homo sapiens" -species_type "Vertebrates" -output_option 3

# Remueve el log de la base de datos derby
rm ~/derby.log

# Identifica los peptidos en la muestra, los caracteriza y asigna un valor de probabilidad
java -Xmx8000m -cp ~/PeptideShaker/PeptideShaker-0.38.3.jar eu.isas.peptideshaker.cmd.PeptideShakerCLI -replicate 3 -experiment "${name}" -sample "${name}" -identification_files ~/$2/ -out ~/$2/${name}.cps -species 'Human (Homo sapiens)' -species_type 'Vertebrates' -min_peptide_length 6 -species_update 1 -id_params p.parameters

# Reporta los peptidos identificados y los datos asociados.
java -Xmx8000m -cp ~/PeptideShaker/PeptideShaker-0.38.3.jar eu.isas.peptideshaker.cmd.ReportCLI -in ~/$2/${name}.cps -out_reports ~/$2/ -reports "0,1,2,3,4" -documentation "0,1,2,3,4"

# Reporta las proteinas identificadas
grep -E '^[0-9]{1,10}\s' ~/$2/*${2}_3_Default_Hierarchical_Report.txt > ${2}_Proteins.txt

# Extrae los códigos EC de UNIPROT
Rscript UniProt.R ${2}_Proteins.txt

exit
