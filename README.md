# Generación de filogenias con UShER a partir de .annoF (sin rescate) y multifasta (con rescate)

## Cómo ejecutar la pipeline

> [!NOTE]  
> Dataset de prueba disponible en `salas:/data/fcordero/usher_example/data/`

1. Instalación de `snakemake`:

```shell
# Instalación via conda
conda create -c conda-forge -c bioconda -c nodefaults -n snakemake snakemake
conda activate snakemake

# Instalación via pip
pip install snakemake
```

2. Establecer ajustes necesarios en `config/config.yaml`:
   - Para **modo rescue** (con rescate): establecer `pipeline_mode: "rescue"` y proporcionar `multifasta` y `snp_table`
   - Para **modo no-rescue** (sin rescate): establecer `pipeline_mode: "no-rescue"` y proporcionar `annof_dir` (directorio con archivos `.annoF`)

3. Ejecutar `snakemake` desde directorio parental (el que incluye el archivo `snakefile`):

```shell
snakemake --use-conda --cores <cores>
```

<br>

## Descripción de la pipeline

**UShER** toma como input, entre otros, un archivo `.diff` que resume todas las variatnes de cada una de las muestras en la filogenia. Para llegar a este archivo, se pueden seguir 2 procesos distintos, dependiendo de si queremos hacer la filogenia con o sin rescate de variantes.

### .diff a partir de multifasta (con rescate)

Archivos necesarios:

- Multi-FASTA
- SNP table
- Ref FASTA
- Bed mask (`.bed` con las regiones enmascaradas del genoma) 

---

#### Paso 1

A partir del multifasta de interés (generado con `ThePipeline_vX consensus`) separamos el archivo Multi-FASTA en FASTAs individuales:

```shell
./split_multifasta.sh <multifasta> <output_dir> 
```

#### Paso 2

Generamos un *pseudoalineamiento* de la referencia contra nuestro *pseudogenoma* con las variantes incluidas. Además de los Multi-FASTA generados en el paso anterior para cada muestra, utilizamos la SNP table (generada con `ThePipeline_vX consensus`).

El script toma la referencia, el FASTA individual y la SNP table y sustituye las posiciones de la referencia por el nucleótido correspondiente en el FASTA individual ya rescatado:

```shell
python align_ref_and_rescued_fasta.py <ref_fasta> <sample_fasta> <snp_table>

# Ejemplo paralelizando con xargs
ls *.fasta | xargs -I {} -P 10 python align_ref_and_rescued_fasta.py /data/Databases/MTB_ancestor/MTB_ancestor_reference.fasta {} SNP_table.txt
```

#### Paso 3

A partir de cada uno de los *pseudoalineamientos* generados en el paso anterior se genera un archivo `.vcf` usando la herramienta `faToVcf` (incluida con `UShER`):

```shell
ls *_aligned.fasta | cut -f 1 -d '.' | xargs -I {} faToVcf {}\.fasta {}\.vcf
```

#### Paso 4

Se convierten los archivos `.vcf` en `.diff`, que contienen las variantes del `.vcf` pero manteniendo solo la información de la posición y el cambio:

```shell
python vcf_to_diff.py -v <vcf> -d <vcf_dir> -smf <bed_mask>

# Ejemplo paralelizando con xargs
ls *.vcf | xargs -I {} -P 10 python vcf_to_diff.py -v {} -d diffs/ -smf /nas01/Fran/UCSC/cryptic_tb_callable_mask/R00000039_repregions.bed
```

#### Paso 5 

Se concatenan los `.diff` en un único archivo:

```shell
cat *.diff > final_concatenated.diff
```

<br>

### .diff a partir de .annoF (sin rescate)

Archivos necesarios:

- `.annoF` de cada muestra
- Ref FASTA

---

#### Paso 1

Para cada uno de los `.annoF`, se genera un alineamiento de la referencia contra un *pseudogenoma* de la muestra, obtenido al copiar la escuencia de referencia y modificar los nucleótidos que aparecen como variantes en cada `.annoF`:

```shell
python individual_multifasta.py <ref> <annoF>

# Ejemplo paralelizando con xargs
ls *annoF | xargs -I {} -P 10 python individual_multifasta.py /data/Databases/MTB_ancestor/MTB_ancestor_reference.fasta {}
```

#### Paso 2

A partir de los archivos generados en el paso anterior, se obtienen los archivos en formato `.vcf` usando la herramienta `faToVcf`:

```shell
ls *_and_ref_multifasta.fas | cut -f 1 -d '.' | xargs -I {} -P 5 faToVcf {}\.fas {}\.vcf
```

#### Paso 3

A partir de los `.vcf` se generan los `.diff`:

```shell
python vcf_to_diff.py -v <vcf> -d <vcf_dir> -smf <bed_mask>

# Ejemplo paralelizando con xargs
ls *.vcf | xargs -I {} -P 10 python vcf_to_diff.py -v {} -d diffs/ -smf /nas01/Fran/UCSC/cryptic_tb_callable_mask/R00000039_repregions.bed
```

#### Paso 4

Se concatenan los `.diff` en un único archivo:

```shell
cat *.diff > final_concatenated.diff
```

<br>

### Reconstrucción de filogenia usando UShER

Tras obtener el archivo `.diff` con todas las muestras concatenadas se sigue el mismo proceso, con y sin rescate:

#### Paso 1

Es necesario crear un *mock tree* y un *mock diff*. Estos archivos se anotarán con la referencia, actuando como una base/plantilla sobre la que `UShER` añadirá el resto de muestras:

```shell
echo "();" > base_tree.nwk
echo ">ref" > ref.diff
```

El árbol *plantilla* se convierte luego a formato `.pb` (Protobuf), más eficiente: 

```shell
usher-sampled -t base_tree.nwk --diff ref.diff --ref /data/Databases/MTB_ancestor/MTB_ancestor_reference.fasta -o ref.pb
```

#### Paso 2

`UShER` añadirá sobre la base `ref.pb` el resto de muestras, incluidas en el `.diff` concatenado:

```shell
usher-sampled -i ref.pb --diff final_concatenated.diff --ref <ref_fasta> -o final_tree.pb -u -d usher_output/ -T 10
```

En el directorio indicado con `-d`, obtendremos la filogenia final generada en formato Newick: `uncondensed-final-tree.nh`.

> [!NOTE]  
> Con `snakemake` el árbol final tendrá el nombre establecido en `config/config.yaml`.
