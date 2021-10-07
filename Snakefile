"""
Calculate variant interface stability impact given a list of
"""
import os

localrules:
    all, complex_combine

COMPLEXES = [i for i in os.listdir('data/complex') if os.path.isdir(f'data/complex/{i}') and os.path.isfile(f'data/complex/{i}/model.pdb') and os.path.isfile(f'data/complex/{i}/individual_list')]

rule all:
    """
    Full pipeline
    """
    input:
        "data/complex.tsv"

rule repair:
    """
    Repair complex PDB file to use with FoldX
    """
    input:
        'data/complex/{complex}/model.pdb'

    output:
        'data/complex/{complex}/model_Repair.pdb',
        'data/complex/{complex}/model_Repair.fxout'

    resources:
        mem_mb = 8000

    log:
        "logs/complex_repair/{complex}.log"

    shell:
        "foldx --command=RepairPDB --pdb=model.pdb --pdb-dir=data/complex/{wildcards.complex} --clean-mode=3 --output-dir=data/complex/{wildcards.complex} &> {log}"

rule wt_analysis:
    """
    Analyse a WT complex interface to identify interface residues
    """
    input:
        pdb='data/complex/{complex}/model_Repair.pdb'

    output:
        'data/complex/{complex}/{interface}/wt/Indiv_energies_model_Repair_AC.fxout',
        'data/complex/{complex}/{interface}/wt/Interaction_model_Repair_AC.fxout',
        'data/complex/{complex}/{interface}/wt/Interface_Residues_model_Repair_AC.fxout',
        'data/complex/{complex}/{interface}/wt/Summary_model_Repair_AC.fxout'

    resources:
        mem_mb = 8000

    log:
        'logs/complex_wt_analysis/{complex}_{interface}.log'

    run:
        root = f'data/complex/{wildcards.complex}/{wildcards.interface}'
        shell(f'mkdir {root} &> {log} && echo "mkdir {root}" &> {log} || true')
        shell(f'mkdir {root}/wt &> {log} && echo "mkdir {root}/wt" &> {log} || true')
        shell(f"foldx --command=AnalyseComplex --pdb=model_Repair.pdb --pdb-dir=data/complex/{wildcards.complex} --clean-mode=3 --output-dir=data/complex/{wildcards.complex}/{wildcards.interface}/wt --analyseComplexChains={wildcards.interface.replace('_', ',')} &> {log}")

rule complex_mutant_models:
    """
    Generate PDBs with each mutation in the interface using FoldX BuildModel
    """
    input:
        muts='data/complex/{complex}/{interface}/individual_list',
        pdb='data/complex/{complex}/model_Repair.pdb'

    output:
        directory('data/complex/{complex}/{interface}/mutant_pdbs'),
        'data/complex/{complex}/{interface}/mutant_models_made'

    resources:
        mem_mb = 16000

    log:
        'logs/complex_mutant_models/{complex}_{interface}.log'

    shell:
        """
        mkdir data/complex/{wildcards.complex}/{wildcards.interface}/mutant_pdbs &> {log}
        foldx --command=BuildModel --pdb-dir=data/complex/{wildcards.complex} --pdb=model_Repair.pdb --mutant-file={input.muts} --output-dir=data/complex/{wildcards.complex}/{wildcards.interface}/mutant_pdbs --numberOfRuns=1 --clean-mode=3 --out-pdb=true &> {log}
        rm data/complex/{wildcards.complex}/{wildcards.interface}/mutant_pdbs/WT_* &> {log}
        touch 'data/complex/{wildcards.complex}/{wildcards.interface}/mutant_models_made' &> {log}
        """

def get_mut_complex_ram(wildcards):
    """
    Get RAM to use for each complex. Add large complexes here or modify to adapt RAM on file size.
    """
    ram = {'6x29': 64000, '6zoj': 64000}
    return ram.get(wildcards.complex, 16000)

rule complex_mut_analysis:
    """
    Analyse mutant interfaces using FoldX AnalyseComplex.
    """
    input:
        'data/complex/{complex}/{interface}/mutant_models_made'

    output:
        directory('data/complex/{complex}/{interface}/mutant'),
        'data/complex/{complex}/{interface}/mutant_analysis_done'

    resources:
        mem_mb = get_mut_complex_ram,

    threads: 8

    log:
        'logs/complex_mut_analysis/{complex}_{interface}.log'

    shell:
        """
        python bin/mut_analysis.py --processes 4 data/complex/{wildcards.complex}/{wildcards.interface}/mutant_pdbs {wildcards.interface} data/complex/{wildcards.complex}/{wildcards.interface}/mutant &> {log}
        touch data/complex/{wildcards.complex}/{wildcards.interface}/mutant_analysis_done &> {log}
        """

rule complex_combine:
    """
    Combine output files from running complex_mut_analysis
    """
    input:
        flag='data/complex/{complex}/{interface}/mutant_analysis_done',
        mutants='data/complex/{complex}/{interface}/individual_list',
        wt_interface='data/complex/{complex}/{interface}/wt/Interface_Residues_model_Repair_AC.fxout',
        wt_interaction='data/complex/{complex}/{interface}/wt/Interaction_model_Repair_AC.fxout',
        wt_indiv='data/complex/{complex}/{interface}/wt/Indiv_energies_model_Repair_AC.fxout',
        wt_summary='data/complex/{complex}/{interface}/wt/Summary_model_Repair_AC.fxout'

    output:
        indiv='data/complex/{complex}/{interface}/individual_energies.tsv',
        interaction='data/complex/{complex}/{interface}/interactions.tsv',
        interface='data/complex/{complex}/{interface}/interface_residues.tsv',
        summary='data/complex/{complex}/{interface}/summary.tsv'

    log:
        'logs/complex_combine/{complex}_{interface}.log'

    run:
        shell(f"python bin/combine.py {input.mutants} {input.wt_indiv} data/complex/{wildcards.complex}/{wildcards.interface}/mutant > {output.indiv} 2> {log}")
        shell(f"python bin/combine.py {input.mutants} {input.wt_interaction} data/complex/{wildcards.complex}/{wildcards.interface}/mutant > {output.interaction} 2> {log}")
        shell(f"python bin/combine.py {input.mutants} {input.wt_interface} data/complex/{wildcards.complex}/{wildcards.interface}/mutant > {output.interface} 2> {log}")
        shell(f"python bin/combine.py {input.mutants} {input.wt_summary} data/complex/{wildcards.complex}/{wildcards.interface}/mutant > {output.summary} 2> {log}")

def get_complex_tsv_files(complex):
    """
    Determine required input for complex_tsv.
    Modify this function to analyse different interfaces (currently just returns A/B)
    """
    return [f'data/complex/{complex}/A_B/interactions.tsv', f'data/complex/{complex}/A_B/interface_residues.tsv']

rule complex_tsv:
    """
    Combine complex FoldX results into one tsv
    """
    input:
        [get_complex_tsv_files(c) for c in COMPLEXES]

    output:
        'data/complex.tsv'

    log:
        'logs/complex_tsv.log'

    shell:
        f'python bin/tsv.py {" ".join(f"data/complex/{c}/model.yaml" for c in COMPLEXES)} > {{output}} 2> {{log}}'
