from pathlib import Path

# Load configuration file
configfile: "config/config.yaml"

# Define pipeline mode
pipeline_mode = config["pipeline_mode"]  # "rescue" or "no-rescue"

# Define input and output paths - convert to absolute paths
ref_fasta = str(Path(config["ref_fasta"]).resolve())
bed_mask = str(Path(config["bed_mask"]).resolve())
output_dir = str(Path(config["output_dir"]).resolve())
usher_output_dir = str(Path(config["usher_output_dir"]).resolve())

# Mode-specific inputs
if pipeline_mode == "rescue":
    multifasta = str(Path(config["multifasta"]).resolve())
    snp_table = str(Path(config["snp_table"]).resolve())
elif pipeline_mode == "no-rescue":
    annof_dir = str(Path(config.get("annof_dir", "annof")).resolve())  # Directory containing .annoF files
else:
    raise ValueError(f"Invalid pipeline_mode: {pipeline_mode}. Must be 'rescue' or 'no-rescue'")

# Define subdirectories for outputs
split_dir = f"{output_dir}/split"
vcf_dir = f"{output_dir}/vcf"
diff_dir = f"{output_dir}/diff"

# Define final output names
final_diff_file = str(Path(config["final_diff"]).resolve())
final_pb_tree = str(Path(config["final_pb_tree"]).resolve())
final_nwk_tree = str(Path(config["final_nwk_tree"]).resolve())

# Helper function to get sample names based on pipeline mode
def get_samples():
    """Extract sample names from input files"""
    if pipeline_mode == "rescue":
        # Parse multifasta to get sample names
        samples = []
        with open(multifasta, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    sample_name = line[1:].strip().split()[0]  # Get first word after '>'
                    samples.append(sample_name)
        return samples
    else:  # no-rescue
        # Get sample names from .annoF files
        samples = []
        for f in Path(annof_dir).glob("*.annoF"):
            sample = f.name.split(".")[0]
            samples.append(sample)
        return samples

# Get samples list at the start
SAMPLES = get_samples()


# Final target
rule all:
    input:
        final_nwk_tree


# ===== RESCUE MODE RULES =====

if pipeline_mode == "rescue":
    # Step 1 (rescue): Split multifasta into individual fasta files
    rule split_multifasta:
        input:
            multifasta
        output:
            temp(expand(split_dir + "/{sample}.fasta", sample=SAMPLES))
        params:
            split_dir=split_dir
        shell:
            """
            bash code/split_multifasta.sh {input} {params.split_dir}
            """

    # Step 2 (rescue): Align each individual rescued fasta to reference genome
    rule align_ref_and_rescued_fasta:
        input:
            ref_fasta=ref_fasta,
            sample_fasta=split_dir + "/{sample}.fasta",
            snp_table=snp_table
        output:
            aligned_fasta=temp(split_dir + "/{sample}_aligned.fasta")
        threads: 1
        shell:
            """
            python code/align_ref_and_rescued_fasta.py {input.ref_fasta} {input.sample_fasta} {input.snp_table}
            """


# ===== NO-RESCUE MODE RULES =====

if pipeline_mode == "no-rescue":
    # Step 1 (no-rescue): Process .annoF files to generate multifasta files
    rule process_annof:
        input:
            annof = annof_dir + "/{sample}.EPI.snp.final.annoF"
        output:
            multifasta=temp(annof_dir + "/{sample}_and_ref_multifasta.fas")
        params:
            ref_fasta=ref_fasta
        threads: 1
        shell:
            """
            python code/individual_multifasta.py {params.ref_fasta} {input.annof}
            """


# ===== SHARED RULES (mode-dependent inputs) =====

# Step 3: Convert fastas to VCF files (handles both rescue and no-rescue modes)
rule fasta_to_vcf:
    input:
        fasta=lambda wildcards: (
            split_dir + f"/{wildcards.sample}_aligned.fasta" if pipeline_mode == "rescue"
            else annof_dir + f"/{wildcards.sample}_and_ref_multifasta.fas"
        )
    output:
        vcf=vcf_dir + "/{sample}.vcf"
    conda:
        "envs/usher-env.yaml"
    params:
        vcf_dir=vcf_dir
    threads: 1
    shell:
        """
        faToVcf {input.fasta} {output.vcf}
        """


# Step 4: Convert VCF files to diff files
rule vcf_to_diff:
    input:
        vcf=vcf_dir + "/{sample}.vcf"
    output:
        diff=diff_dir + "/{sample}.diff"
    params:
        bed_mask=bed_mask,
        diff_dir=diff_dir
    threads: 1
    shell:
        """
        mkdir -p {params.diff_dir}
        python code/vcf_to_diff.py \
            -v {input.vcf} -d {params.diff_dir} -smf {params.bed_mask}
        """


# Step 5: Merge all .diff files into a single file
rule merge_diffs:
    input:
        diffs=expand(diff_dir + "/{sample}.diff", sample=SAMPLES)
    output:
        final_diff_file
    params:
        diff_dir=diff_dir
    shell:
        """
        find {params.diff_dir} -name "*.diff" -type f | sort | xargs cat > {output}
        """

# Step 6: Create tree template files
rule create_tree_template:
    output:
        tree=temp(usher_output_dir + "/base_tree.nwk"),
        diff=temp(usher_output_dir + "/ref.diff")
    shell:
        """
        echo "();" > {output.tree}
        echo ">ref" > {output.diff}
        """

# Step 7: Convert template to protobuf format
rule convert_template_to_pb:
    input:
        tree=usher_output_dir + "/base_tree.nwk",
        diff=usher_output_dir + "/ref.diff",
        ref_fasta=ref_fasta
    output:
        pb=temp(usher_output_dir + "/ref.pb")
    conda:
        "envs/usher-env.yaml"
    shell:
        """
        cd {usher_output_dir} && usher-sampled -t base_tree.nwk --diff ref.diff --ref {input.ref_fasta} -o ref.pb
        """

# Step 8: Build final phylogenetic tree
rule build_final_tree:
    input:
        pb=usher_output_dir + "/ref.pb",
        diff_file=final_diff_file,
        ref_fasta=ref_fasta
    output:
        pb_tree=final_pb_tree,
        nwk_tree=final_nwk_tree
    params:
        usher_dir=usher_output_dir
    threads: workflow.cores
    conda:
        "envs/usher-env.yaml"
    shell:
        """
        cd {params.usher_dir} && usher-sampled -i ref.pb --diff {input.diff_file} --ref {input.ref_fasta} \
            -o {output.pb_tree} -u -d . -T {threads}
        mv {params.usher_dir}/uncondensed-final-tree.nh {output.nwk_tree}
        """