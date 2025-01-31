rule all:
    input:
        auspice_json = "auspice/WNV-Foco_nextstrain2.json",

rule options:
	params:
		threads = 8

options = rules.options.params

input_fasta = "data/sequences.fasta",
input_metadata = "data/metadata.csv",
dropped_strains = "config/dropped_strains.txt",
mask_index = "config/masked_sites.txt",
reference = "config/reference.gb",
colors = "config/colors.tsv",
lat_longs = "config/lat_longs.tsv",
auspice_config = "config/auspice_config.json",
input_clades = "config/clades.tsv"



rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = input_fasta,
        reference = reference
    output:
        alignment = "results/aligned.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps
        """

rule mask:
    message:
        """
        masking sites in {input.mask}
        """
    input:
        fasta = rules.align.output.alignment,
        mask = mask_index
    output:
        masked_alignment = "results/masked.fasta"
    shell:
        """
        augur mask \
            --sequences {input.fasta} \
            --mask {input.mask} \
            --output {output.masked_alignment}
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.mask.output.masked_alignment
    output:
        tree = "results/tree_raw.nwk"
    shell:
        """
       augur tree \
           --alignment {input.alignment} \
           --method iqtree \
           --override-default-args \
           --tree-builder-args="-ninit 10 -n 4"\
           --nthreads auto \
           --output {output.tree}
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock
expectation
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = input_metadata
    output:
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 100
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --stochastic-resolve \
            --coalescent {params.coalescent} \
            --date-confidence True\
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd} \
            --seed 1 \
            --root AF481864/1998/Israel/ISR/D/Eilat/pre-NY/unknown/Bird-other
        """
#remove/add ## to --seed to ensure same/stochastic results for refine

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "results/nt_muts.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = reference
    output:
        node_data = "results/aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
        """

rule clades:
    message: " Labeling clades as specified in config/clades.tsv"
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = input_clades
    output:
        clade_data = "results/clades.json"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output {output.clade_data}
        """

rule traits:
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = input_metadata
    output:
        node_data = "results/traits.json",
    params:
        columns = "country state CO division zone trap"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output-node-data {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = input_metadata,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        colors = colors,
        lat_longs = lat_longs,
        clades = rules.clades.output.clade_data,
        auspice_config = auspice_config
    output:
        auspice_json = rules.all.input.auspice_json,
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} {input.clades} \
            --lat-longs {input.lat_longs} \
            --colors {input.colors} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json}
        """

# would need to add the following line above             --clades {input.clades} {input.aa_muts} \

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
